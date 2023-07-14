#!/usr/bin/env nextflow

params.datapath="${launchDir}"
params.outdir = "${launchDir}"
params.ms_files= null
params.solsdir="${launchDir}"

// params.number_of_threads = 4

params.command = null

workflow {
    // Check that we have the 'params.ms_files' list and load it into a channel
    _ = file(params.ms_files, glob: false, checkIfExists: true) //check it exists
    def msetsList = new File(params.ms_files).collect {it}
    msets_ch = Channel.fromList(msetsList) // .buffer( size: 4, remainder: true )

    sage_std_ch = runSagecalStandalone(msets_ch, params.shapelets.modes, params.solsdir).collect()
}

process runSagecalStandalone {
    debug true
    // cpus params.number_of_threads
    // errorStrategy { task.exitStatus == 139 ? 'retry' : 'terminate' } //134 is core dumping error
    time '2h'
    errorStrategy 'retry'
    maxRetries 4
    
    publishDir "${params.outdir}/bandpass_logs", pattern: "*bandpass.log", mode: "move", overwrite: true
    publishDir "${solsdir}", pattern: "*.solutions", mode: "move", overwrite: true
    label 'parallel_jobs'

    input:
    path ms
    val modes
    val solsdir

    output:

    tuple path("${ms.Name}.solutions"), path("${ms.Name}_bandpass.log"), val(true)

    script:
    standalone_sagecal_command = make_standalone_sagecal_command() + " -p ${ms.Name}.solutions -d ${ms}"

    """
    cp ${modes} \$PWD
    ${standalone_sagecal_command} > "${ms.Name}_bandpass.log"  2>&1
    rm *.modes
    """
}


workflow Predict {
    // Check that we have the 'params.ms_files' list and load it into a channel
    _ = file(params.ms_files, glob: false, checkIfExists: true) //check it exists
    def msetsList = new File(params.ms_files).collect {it}
    msets_channel = Channel.fromList(msetsList)

    predict_channel = runSagecalStandalonePredict(msets_channel, params.shapelets.modes).collect()
}


process runSagecalStandalonePredict {
    debug true
    errorStrategy 'retry'
    maxRetries 2
    cpus 12
    memory '100 GB'
    // label 'parallel_jobs'
    // publishDir "${params.outdir}/logs", pattern: "*predict.log", mode: "move", overwrite: true

    input:
    path ms
    val modes

    output:
    path "*predict.log" // tuple  val(true)

    script:
    standalone_sagecal_command = make_standalone_sagecal_command() + " -d ${ms}"

    // cp ${modes} \$PWD
    """
    ${standalone_sagecal_command} > "${ms.Name}_predict.log"  2>&1
    """
}


def make_standalone_sagecal_command() {
    return params.command
}

