#!/usr/bin/env nextflow

params.path = null
params.ms_files="ms_files_002.txt"
params.data_column=null

process gen003Vis {
    publishDir "${params.path}", mode: 'move'

    input:
    path ms002
    val data_column


    output:
    tuple path('*_SAP*_SB*_uv_003*.MS'), val(true)

    shell:
    '''
    #!/bin/bash
    ms003=$(echo "!{ms002}" | sed "s/002/003/")
    cat >"!{ms002.Name}_avg_to_003.parset" <<EOL
    msin = !{ms002}
    msout = ${ms003}
    msin.startchan = 0
    msin.nchan = 3
    msin.datacolumn = "!{data_column}"
    msout.datacolumn = DATA
    steps = [avg1]
    avg1.type = average
    avg1.freqstep = 1
    avg1.timestep = 5
    EOL
    DP3 "!{ms002.Name}_avg_to_003.parset"
    '''
    //same as writing this inside a bash script and calling it thisway
    // template '/home/users/chege/theleap/leap/templates/avgto003.sh'
}

workflow {
    _ = file(params.ms_files, glob: false, checkIfExists: true) //check it exists
    def msetsList = new File(params.ms_files).collect {it}
    msets_ch = Channel.fromList(msetsList)
    gen3vis_ch = gen003Vis(msets_ch, params.data_column).collect()
}