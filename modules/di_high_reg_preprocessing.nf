#!/usr/bin/env nextflow

// params.outdir = "${launchDir}"
params.ms_files = null


process copyFlags {
    debug true
    //cpus 12

    input:
    path ms

    output:
    path "${ms}"

    // sleep 60s
    script:
    """
    python3  /home/users/chege/theleap/leap/templates/cp_flag_col.py -i ${ms} -t table.f4_TSM0 -o table.f4_TSM0_copy
    """
}


process copyFlagsBack {
    debug true
    //cpus 12

    input:
    path ms

    output:
    path "${ms}"

    // sleep 60s
    script:
    """
    python3  /home/users/chege/theleap/leap/templates/cp_flag_col.py -i ${ms} -t table.f4_TSM0_copy -o table.f4_TSM0
    """
}


process scaleData {
    debug true
    //cpus 12

    input:
    path ms

    output: 
    path "${ms}"

    script:
    """
    python3 /home/users/chege/theleap/leap/templates/scaledata.py -i ${ms} -f 0
    """

}


process clipData {
    debug true
    //cpus 12

    input:
    path ms

    output:
    path "${ms}"


    """
    python3 /home/users/chege/theleap/leap/templates/clip_data.py -i ${ms} --flag_intrastations --flag_badbaselines -o DATA 
    """
}

workflow {
    _ = file(params.ms_files, glob: false, checkIfExists: true) //check it exists
    def msetsList = new File(params.ms_files).collect {it}
    msets_ch = Channel.fromList(msetsList)
    copy_flags_ch = copyFlags(msets_ch)
    copy_flags_back_ch = copyFlagsBack(copy_flags_ch)
    scale_data_ch = scaleData(copy_flags_back_ch)
    clip_ch = clipData(scale_data_ch).collect()
}



// process fixBeamInfo {
//     input:
//     path ms

//     script:
//     """
//     python3 /home/users/chege/pipeline/fixbeaminfo.py ${ms}
//     """
// }