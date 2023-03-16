#!/usr/bin/env nextflow

params.datapath="${launchDir}"
params.outdir = "${launchDir}"
params.msfiles= 'ms_files.txt'

params.number_of_threads = 4

params.command = null

workflow {
    // Check that we have the 'params.msfiles' list and load it into a channel
    _ = file(params.msfiles, glob: false, checkIfExists: true) //check it exists
    def msetsList = new File(params.msfiles).collect {it}
    msets_ch = Channel.fromList(msetsList) // .buffer( size: 4, remainder: true )

    sage_std_ch = runSagecalStandalone(msets_ch, params.shapelets_modes).collect()
}

process runSagecalStandalone {
    debug true
    cpus params.number_of_threads
    publishDir params.outdir

    input:
    path ms
    val modes

    output:

    tuple path("${ms.Name}.solutions"), path("${ms.Name}_bandpass.log"), val(true)

    script:
    standalone_sagecal_command = make_standalone_sagecal_command() + " -p ${ms.Name}.solutions -d ${ms}"

    """
    cp ${modes} \$PWD
    ${standalone_sagecal_command} > "${ms.Name}_bandpass.log"  2>&1
    """
}

def make_standalone_sagecal_command() {
    return params.command
}
