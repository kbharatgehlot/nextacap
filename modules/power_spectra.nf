#!/usr/bin/env nextflow

params.datapath = "/data/users/lofareor/chege/psp/test" //null
params.revision="rev001" //null
params.obsid="L254871" //null
params.msfiles = "/data/users/lofareor/chege/psp/ms_files.txt" //null
params.nodes = "node129" //"node[125-129]"
params.max_concurrent=15


workflow {
    init_ch = initPSDB(true, params.datapath)
    rev_ch = addRevision(init_ch.ps_dir, init_ch.default_toml_file, params.nodes, params.max_concurrent, params.revision )
    ps_ch = runPSPIPE(init_ch.ps_dir, rev_ch, params.obsid, params.msfiles)
    plotPS(ps_ch.ready, init_ch.ps_dir, rev_ch, params.obsid)
}


process initPSDB {
    // publishDir params.datapath

    input:
    val ready
    path datapath

    output:
    path "${datapath}/ps" , emit: ps_dir
    path "${datapath}/ps/default.toml", emit: default_toml_file
    // path "${datapath}/ps/config"

    script:

    """
    mkdir "${datapath}/ps"
    cd "${datapath}/ps"
    psdb init hba "${datapath}/ps" default
    """
}


process addRevision{
    publishDir "${params.datapath}/ps"

    input:
    path ps_dir
    path default_toml_file
    val nodes
    val max_concurrent
    val revname

    output:
    path "${revname}.toml"

    shell:
    '''
    #!/bin/bash
    psdb new_rev !{default_toml_file} !{revname}
    cat >"!{revname}.toml" <<EOL
    default_settings = "!{ps_dir}/!{default_toml_file}"
    [worker]
    nodes = \"!{nodes}\"
    max_concurrent = !{max_concurrent}
    run_on_file_host = true
    run_on_file_host_pattern = '\\/net/(node\\d{3})'
    [image]
    data_col = "DATA"
    [power_spectra]
    eor_bin_list = "!{ps_dir}/config/eor_bins_hba.parset"
    ps_config = "!{ps_dir}/config/ps_config_hba.parset"
    flagger = "!{ps_dir}/config/flagger.parset"
    [gpr]
    name = ""
    plot_results = true
    config_i = "!{ps_dir}/config/gpr_config_hba.parset"
    config_v = "!{ps_dir}/config/gpr_config_v.parset"
    EOL
    '''
    //same as writing this inside a bash script and calling it thisway
    // template '/home/users/chege/theleap/leap/templates/add_rev.sh'
}


process runPSPIPE {
    // publishDir "${params.datapath}/ps"

    input:
    path ps_dir

    path toml_file
    val obsid
    path msfiles

    output:
    val true, emit: ready
    // path "default"
    // path "merged_ms"
    // path "obs_ids"

    shell:
    // ,ssins seems to give an error when included. TODO:
    '''
    psdb add_obs !{toml_file} !{obsid} -m !{msfiles}
    pspipe merge_ms !{toml_file} !{obsid}
    pspipe image,gen_vis_cube !{toml_file} !{obsid}_flagged
    pspipe run_gpr !{toml_file} !{obsid}_flagged
    '''
}


process plotPS {
    publishDir "${params.datapath}/ps"

    input:
    val ready
    path ps_dir
    path toml_file
    val obsid

    output:
    path "${obsid}_flagged_IV_power_spectrum.png"

    script:
    """
    python3 /home/users/chege/theleap/leap/templates/plot_power_spctrum.py ${toml_file} --obsid ${obsid}_flagged
    """
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// I initilally wrote these tasks as separate processes but since they all take similar inputs I combined them into the runPSPIPE process
// I might have to add some if conditions incase say, the user wants to skip gpr
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// add_ch = addObs(rev_ch, params.obsid, params.msfiles)
// merge_ch = mergeMsets(add_ch, rev_ch, params.obsid)
// cube_ch = makeVisCube(merge_ch, rev_ch, params.obsid)
// gpr_ch = runGPR(cube_ch, rev_ch, params.obsid)


// process addObs {
//     input:
//     path toml_file
//     val obsid
//     path msfiles

//     output:
//     val true

//     script:
//     """
//     psdb add_obs ${toml_file} ${obsid} -m ${msfiles}
//     """
// }


// process mergeMsets {
//     input:
//     val ready
//     path toml_file
//     val obsid

//     output:
// process addObs {
//     val true
    
//     //,ssins
//     script:
//     """
//     pspipe merge_ms ${toml_file} ${obsid}
//     """
// }


// process makeVisCube {
//     input:
//     val ready
//     path toml_file
//     val obsid

//     output:
//     val true

//     script:
//     """
//     pspipe image,gen_vis_cube ${toml_file} ${obsid}_flagged
//     """
// }

// process runGPR {
//     input:
//     val ready
//     path toml_file
//     val obsid

//     output:
//     val true

//     script:
//     """
//     pspipe run_gpr ${toml_file} ${obsid}_flagged
//     """
// }