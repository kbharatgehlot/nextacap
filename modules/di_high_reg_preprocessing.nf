#!/usr/bin/env nextflow

params.outdir = "${launchDir}"

process copyFlagsBack {
    input:
    path ms

    script:
    """
    python3  /home/users/chege/theleap/leap/templates/cp_flag_col.py -i ${ms} -t table.f4_TSM0_copy -o table.f4_TSM0
    sleep 100s
    """
}

process clipData {
    debug true
    cpus 12

    input:
    path ms

    output:
    val ms

    """
    python3 /home/users/chege/theleap/leap/templates/clip_data.py -i ${ms} --flag_intrastations --flag_badbaselines --flag_longbaselines -o DATA
    """
}
    
workflow {
    _ = file(params.msfiles, glob: false, checkIfExists: true) //check it exists
    def msetsList = new File(params.msfiles).collect {it}
    msets_ch = Channel.fromList(msetsList)
    copyFlagsBack(msets_ch)
    clip_ch = clipData(msets_ch).collect()
}



// process copyFlags {
//     input:
//     path ms

//     script:
//     """
//     python3  /home/users/chege/cp_flag_col.py -i ${ms}
//     """
// }

// process scaleData {
//     input:
//     path ms

//     script:
//     """
//     python3 /home/users/chege/pipeline/scaledata.py ${ms}
//     """

// }

// process fixBeamInfo {
//     input:
//     path ms

//     script:
//     """
//     python3 /home/users/chege/pipeline/fixbeaminfo.py ${ms}
//     """
// }