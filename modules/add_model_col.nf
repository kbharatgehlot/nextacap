#!/usr/bin/env nextflow

params.ms_files = null
params.datacol = null
params.correctedcol = null
params.modelcol = null

process addModelCol {
    debug true
    cpus 12

    input:
    path ms
    val datacol
    val correctedcol
    val modelcol

    output:
    val ms

    script:
    """
    python3 /home/users/chege/theleap/leap/templates/add_model_col.py -i ${ms} -d ${datacol} -c ${correctedcol} -m ${modelcol}
    """
}

workflow {
    _ = file(params.ms_files, glob: false, checkIfExists: true) //check it exists
    def msetsList = new File(params.ms_files).collect {it}
    msets_ch = Channel.fromList(msetsList)
    addmodel_ch = addModelCol(msets_ch, params.datacol, params.correctedcol, params.modelcol).collect()
}