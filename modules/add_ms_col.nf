#!/usr/bin/env nextflow

params.ms_files = null
params.col = null


process addCol {
    debug true
    cpus 12

    input:
    path ms
    val col

    output:
    val ms

    script:
    """
    python3 /home/users/chege/theleap/leap/templates/add_col.py -i ${ms} -c ${col}
    """

    // alternative method
    // DP3 msin={ms} msin.datacolumn={incol} msout=. msout.datacolumn={outcol} steps=[] 
}

workflow {
    _ = file(params.ms_files, glob: false, checkIfExists: true) //check it exists
    def msetsList = new File(params.ms_files).collect {it}
    msets_ch = Channel.fromList(msetsList)
    clip_ch = addCol(msets_ch, params.col).collect()
}