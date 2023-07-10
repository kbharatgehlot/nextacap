#!/usr/bin/env nextflow

params.datapath = "/data/users/lofareor/chege/psp/test" //null
params.revision="rev001" //null
params.obsid="L254871" //null
params.msfiles = "/data/users/lofareor/chege/psp/ms_files.txt" //null
params.nodes = "node129" //"node[125-129]"
params.max_concurrent=15
params.data_column="DATA"
params.pspipe_dir=null
params.merge_ms = false
params.delay_flag = false
params.ml_gpr = false


workflow {
    // init_ch = InitPSDB(true, params.datapath)
    rev_ch = AddRevision(true, params.obsid, params.data_column, params.datapath, params.pspipe_dir, params.nodes, params.max_concurrent, params.revision, params.merge_ms) //!check this
    ps_ch = RunPSPIPE(params.pspipe_dir, rev_ch.toml_file, params.obsid, params.msfiles, params.merge_ms, params.delay_flag, params.ml_gpr)
    // PlotPowerSpectrum(ps_ch.ready, params.pspipe_dir, rev_ch.toml_file, params.obsid)
}


process InitPSDB {
    // publishDir params.datapath
    input:
    val ready
    path datapath

    output:
    path "${datapath}/ps2" , emit: ps_dir
    // path "${datapath}/ps/default.toml", emit: default_toml_file
    // path "${datapath}/ps/config"

    script:

    """
    mkdir "${datapath}/ps2"
    cd "${datapath}/ps2"
    # psdb init hba "${datapath}/ps" default
    """
}


process AddRevision{
    publishDir "${ps_dir}"

    input:
    val ready
    val obsid
    val data_column
    path datapath
    val ps_dir
    val nodes
    val max_concurrent
    val revname
    val merge_ms

    output:
    val "${ps_dir}/${revname}.toml", emit: toml_file

    shell:
    // template 'add_rev.sh'
    '''
    #!/bin/bash

mkdir -p "!{ps_dir}"
cd "!{ps_dir}"
cp "!{projectDir}/configs/pspipe_toml_templates/default.toml" .
cp "!{projectDir}/configs/pspipe_toml_templates/eor_bins_hba.parset" .
cp "!{projectDir}/configs/pspipe_toml_templates/ps_config_hba.parset" .
cp "!{projectDir}/configs/pspipe_toml_templates/flagger.parset" .
cp "!{projectDir}/configs/pspipe_toml_templates/gpr_config_hba.parset" .
cp "!{projectDir}/configs/pspipe_toml_templates/gpr_config_v.parset" .
cp "!{projectDir}/configs/pspipe_toml_templates/gpr_ml_config_2023.parset" .
cp "!{projectDir}/configs/pspipe_toml_templates/flagger_pre_combine.parset" .

if !{merge_ms}; then
    image_data_col="DATA"
else
    image_data_col="!{data_column}"
fi

psdb new_rev !{ps_dir}/default.toml !{revname}
cat >"!{revname}.toml" <<EOL
default_settings = "!{ps_dir}/default.toml"
data_dir = "!{ps_dir}"
[worker]
nodes = \"!{nodes}\"
max_concurrent = !{max_concurrent}
run_on_file_host = true
run_on_file_host_pattern = '\\/net/(node\\d{3})'
env_file='/home/users/mertens/.activate_pspipe_dev.sh'
[merge_ms]
data_col = "!{data_column}"
aoflagger_strategy = '/home/users/mertens/projects/NCP/nights_np5_red1/lofar-sens2.lua'
[image]
data_col = "${image_data_col}"
channels_out = 'every3'
name="!{revname}"
[power_spectra]
eor_bin_list = "!{ps_dir}/eor_bins_hba.parset"
ps_config = "!{ps_dir}/ps_config_hba.parset"
flagger = "!{ps_dir}/flagger.parset"
[gpr]
config_i = "!{ps_dir}/gpr_config_hba.parset"
config_v = "!{ps_dir}/gpr_config_v.parset"
[ml_gpr]
name = 'eor_vae_2023'
config = "!{ps_dir}/gpr_ml_config_2023.parset"
[combine]
pre_flag = "!{ps_dir}/flagger_pre_combine.parset"
EOL
    '''
}


process RunPSPIPE {
    // debug true

    input:
    path ps_dir
    path toml_file
    val obsid
    path msfiles
    val merge_ms
    val delay_flag
    val vis_flag
    val gpr
    val ml_gpr

    output:
    val true, emit: ready

    shell:
    // template 'run_pspipe.sh'
    // ,ssins seems to give an error when included. TODO:
    // '''
    // psdb add_obs !{toml_file} !{obsid} -m !{msfiles}
    // pspipe merge_ms,delay_flagger !{toml_file} !{obsid}
    // pspipe image,gen_vis_cube !{toml_file} !{obsid}_flagged
    // pspipe run_ml_gpr !{toml_file} !{obsid}_flagged
    // '''
    '''
    psdb add_obs !{toml_file} !{obsid} -m !{msfiles}
    obs="!{obsid}"

    mkdir -p !{launchDir}/logs

    if !{merge_ms}; then
        echo "Merging Ms files"
        if !{delay_flag}; then
            echo "Using delay flagger"
            pspipe merge_ms,delay_flagger !{toml_file} !{obsid} > !{launchDir}/logs/ps_ms_merging_with_aoflagger_and_delay_flagger.log 2>&1
        elif !{vis_flag}; then
            pspipe merge_ms,vis_flagger !{toml_file} !{obsid} > !{launchDir}/logs/ps_ms_merging_with_aoflagger_and_vis_flagger.log 2>&1
        else
            pspipe merge_ms !{toml_file} !{obsid} > !{launchDir}/logs/ps_ms_merging_aoflagger_only.log 2>&1
        fi

        obs="!{obsid}_flagged"
    fi
    echo "making image cube"
    pspipe image,gen_vis_cube !{toml_file} ${obs} > !{launchDir}/logs/ps_image_gen_vis_cube.log 2>&1

    if !{ml_gpr}; then
        echo "Running foreground subtraction with ML_GPR"
        pspipe run_ml_gpr !{toml_file} ${obs} > !{launchDir}/logs/ps_ml_gpr.log 2>&1
    elif !{gpr}; then
        echo "Running foreground subtraction with GPR"
        pspipe run_gpr !{toml_file} ${obs} > !{launchDir}/logs/ps_gpr.log 2>&1
    else
        echo "GPR foreground subtraction NOT applied"
    fi

    # python3 !{projectDir}/templates/plot_power_spctrum.py !{toml_file} --obsid ${obs} --outdir !{ps_dir} > !{launchDir}/logs/plot_ps.log 2>&1

    '''

}


process PlotPowerSpectrum {
    publishDir "${ps_dir}", pattern: "*.png", mode: "move", overwrite: true

    input:
    val ready
    path ps_dir
    path toml_file
    val obsid

    output:
    path "${obsid}_*.png"

    shell:
    """
    python3 ${projectDir}/templates/plot_power_spctrum.py !{toml_file} --obsid !{obsid}
    """
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// I initilally wrote these tasks as separate processes but since they all take similar inputs I combined them into the RunPSPIPE process
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
