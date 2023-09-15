#!/usr/bin/env nextflow

/*
* README
* copys data from similarly named directories from one list of nodes to another list of nodes  based on a wildcard pattern
* USAGE
* ~/mysoftware/nextflow run ~/theleap/distribute_data.nf
*/

//inputs
params.obsid =null   //e.g "L254871"
params.dirA =null //Should be availble across all the nodes from which we are moving the data e.g. "/data/users/lofareor/NCP/pipeline/redshift1/${params.obsid}/"
params.dirB = null //The data will be moved to this path on each of the nodes we are moving the data to. if it does not exist it will be created. if it exists it will be overwritten. e.g "/data/users/lofareor/chege/leap_test_data/"
params.prefix = null // This is the glob pattern of the datato be searched for in each node e.g. "*_uv_001_simulation.MS" e.g. "*_uv_002_simulation_${params.label}.MS"
params.from_nodes = null //NOdes to get the data from. comma separted string with no spaces or just an integer if single node e.g.  '120,121,122,123,125' or 120 
params.to_nodes = null // Nodes to copy the data to e.g.  '120,121,122,123,125' or 120 
params.files_per_node = null //Number of files to be placed per node. The remainder will go to the master node. Make sure that total_ms_files/files_per_node < number_of_to_nodes but also have a remainder for the master node.
params.dry_run = false //Use this to test what happens before actully starting to copy the data

workflow {
    from_nodes = parseNodes(params.from_nodes)
    to_nodes = parseNodes(params.to_nodes)
    println(from_nodes)
    println(to_nodes)

    if (from_nodes.size() == to_nodes.size()){
        nodesA_ch = Channel.fromList( from_nodes )
        nodesB_ch = Channel.fromList( to_nodes )
    }

    else if (from_nodes.size() == 1) {
        nodesA_ch = Channel.fromList( from_nodes )
        nodesB_ch = Channel.fromList( to_nodes ).collect()
    }

    else {
        nodesA_ch = Channel.fromList( from_nodes ).collect()
        nodesB_ch = Channel.fromList( to_nodes ).collect()
    }

    copy = copy_from_nodeA2B(nodesA_ch, nodesB_ch)
}

process copy_from_nodeA2B {
    debug true

    input:
    val nodeA
    val nodeB

    output:
    stdout

    script:
    if (params.dry_run)

    """
    python3 ${projectDir}/templates/migrate_data.py -a ${nodeA} -b ${nodeB} -x ${params.dirA} -y ${params.dirB} -d ${params.prefix} -n ${params.files_per_node} -t
    """

    else

    """
    python3 ${projectDir}/templates/migrate_data.py -a ${nodeA} -b ${nodeB} -x ${params.dirA} -y ${params.dirB} -d ${params.prefix} -n ${params.files_per_node}
    """

}

//================== Closures ============================================
//without appending `node` prefix
def parseNodes( nodes ) {
    if (nodes instanceof String){
        nodes_list = nodes.split(',').collect{"${it}"} as List
    }
    //when a single node is given..
    else if (nodes instanceof Integer) {
        nodes_list = nodes.collect{"${it}"} as List
    }
    return nodes_list
}





// include { PsshNodesList  } from './main.nf'

// params.pssh_hosts_txt_file="${launchDir}/psshhosts.txt"
// params.htmltxt = null //needed when downloading staged data


// def parseNodesWithPrefix( nodes ) {
//     if (nodes instanceof String){
//         nodes_list = nodes.split(',').collect{"node${it}"} as List
//     }
//     //when a single node is given..
//     else if (nodes instanceof Integer) {
//         nodes_list = nodes.collect{"node${it}"} as List
//     }
//     return nodes_list
// }

// process DownloadStagedMS {
//     debug true

//     input:
//     val htmltxt
//     val pssh_hosts_txt_file
//     val datapath

//     output:
//     val true

//     script:
//     downloading_file="/home/users/chege/theleap/leap/download_staged_data.nf"
//     """
//     pssh -v -i -h ${pssh_hosts_txt_file} -t 0 -x "cd ${datapath}; bash" ${params.nextflow_executable} run ${downloading_file} --htmltxt ${htmltxt} -entry Download > ${launchDir}/downloading_staged_data.log 2>&1
//     """
// }


// workflow SplitHtmlTxtFile {
//     chunks_ch = Channel.fromPath(params.htmltxt).splitText( by: params.files_per_node, file: true )
// }


// workflow Download {
//     to_nodes = parseNodesWithPrefix(params.to_nodes)
//     pssh_ch = PsshNodesList(true, to_nodes, params.pssh_hosts_txt_file)

//     chunks_ch = Channel.fromPath(params.htmltxt).splitText( by: params.files_per_node, file: true )
//     // chunks_ch.view()
//     DownloadStagedMS(chunks_ch, pssh_ch[0], params.dirB)

//     // .map { file -> tuple(file.name, file) }
    
//     //    .splitText(by: 2, file: true)
//     // chunks_ch.view()
//     //    .set { chunks_ch }
    
// }
