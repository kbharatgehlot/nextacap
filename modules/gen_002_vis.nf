
#!/usr/bin/env nextflow

params.path = null
params.msfiles=null

params.number_of_threads = 12

process flagGen002Vis {
    publishDir "${params.path}" mode: 'move'

    input:
    path ms001

    output:
    val true

    shell:
    '''
    #!/bin/bash
    ms002==$(echo "!{ms001}" | sed "s/001/002/")
    cat >"!{ms001.Name}_avg_to_003.parset" <<EOL
    msin = !{ms001}
    msout = ${ms002}
    msin.startchan = 0
    msin.nchan = 15
    msin.datacolumn = DATA
    msout.datacolumn = DATA
    steps = [flag1,interpolate,avg1]
    flag1.type=aoflagger
    flag1.memoryperc=20
    avg1.type =average
    avg1.freqstep = 5
    avg1.timestep = 1
    EOL
    DP3 "!{ms001.Name}_flag_avg_to_002.parset"
    '''
    //same as writing this inside a bash script and calling it thisway
    // template '/home/users/chege/theleap/leap/templates/avgto003.sh'
}

workflow {
    _ = file(params.msfiles, glob: false, checkIfExists: true) //check it exists
    def msetsList = new File(params.msfiles).collect {it}
    msets_ch = Channel.fromList(msetsList)
    gen3vis_ch = gen003Vis(msets_ch).collect()
}





