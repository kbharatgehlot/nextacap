#!/usr/bin/env nextflow

// params.outdir = "${launchDir}"
params.ms_files = null

process clipData {
    debug true
    cpus 12

    input:
    path ms

    output:
    val ms
    // This is only done because sagecal does not like these outlying visibilities
    // The flag_intrastations and flag_badbaselines here is not needed
    // Always check the clip data statistics
    """
    python3 /home/users/chege/theleap/leap/templates/clip_data.py -i ${ms} -o DATA
    """
}

workflow {
    _ = file(params.ms_files, glob: false, checkIfExists: true) //check it exists
    def msetsList = new File(params.ms_files).collect {it}
    msets_ch = Channel.fromList(msetsList)
    clip_ch = clipData(msets_ch).collect()
}