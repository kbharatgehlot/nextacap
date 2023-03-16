#!/usr/bin/env nextflow
import groovy.json.JsonOutput

include {
    makeEffNr;
    convertSageSolutions;
    convertSageZSol;
    plotDDSageSols;
    plotDISageSols;
} from './modules/gains_analysis'

include {
    initPSDB;
    addRevision;
    runPSPIPE;
    plotPS;
} from './modules/power_spectra.nf'

params.obsid = "L254871"
params.nodes = "125,126,127,128,129" //"125,126,117,128,129,130,131"
params.slots = 4 //TODO: determine this automatically from given nodes.
params.datapath = "/data/users/lofareor/chege/leap_tests/leaptest_mpi"  //"/data/users/lofareor/chege/spatial_reg_test/1hour/spectral/L254871" ////"/home/users/chege/theleap/pipe/testdir" //"/data/users/lofareor/chege/leap_tests/bandpass" //
params.label="R" //"2022" //"R" //
params.redshift=1
params.timechunks = 1
params.msfiles= "ms_files.txt"
params.all_msets_txt_list = "${params.datapath}/all_ms_files.txt"
params.mpirun_hosts_txt_file = "${launchDir}/mpirun_hosts_list.txt" // nodes host filename for mpirun #${params.datapath}
params.pssh_hosts_txt_file = "${launchDir}/pssh_hosts_list.txt" // nodes host filename for pssh #${params.datapath}/

//pspipe
params.revision= "rev001"
params.max_concurrent = 15


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Help message
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

def helpMessage() {
log.info """
Typical usage command:
    nextflow run ../leap/main.nf --mode mpi --stage DD --obsid L254871 --nodes "125,126,127,128,129"

Mandatory arguments:
--obsid                        Observation ID. e.g L254874
--datapath                     The path where data is stored. Should be valid for all the nodeswhere data is distributed.
--nodes                        Node numbers where data is distributed as a comma separated string. e.g. "125,126,127,128,129" or a single number if using a single node. The last node in the list is used as the masternode.
--slots                        Number of slots per node for mpi run. For now, all given nodes are assigned the same slots number so maybe don't use nodes 119 and 124 on dawn.

--mode                         Either of the 2 sagecal modes "mpi" or "standalone".
--stage                        "DI" or "DD".
--stage_number                 "002" or "003", data labeling depending the analysis stage
--label                        An extra data string label, like "_R" or '_2022' previously used for LOFAR NCP data to refer to a specific processing season.
--redshift                     One of the 3 LOFAR redshift bins; 1, 2 or 3.

Optional arguments:
--number_of_processes          Total processes needed for sagecal mpi mode.
--params.ms_pattern            A pattern matching the naming of the measurement sets to be processed.
--timechunks                   Stop after outputting this number of solutions. Otherwise calibrate full observation. Only available in mpi mode.
--help                         This usage statement, and exit.
"""
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Workflows
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

workflow {
    //Check that all params are valid
    //save them in a json file for modification if needed and for posterity
    VET_PARAMS()

    //Generate mpirun and pssh hosts files
    SET_NODES(VET_PARAMS.out)

    //Call sagecall with an 'all_nodes_standby' signal and the pssh_hosts.txt file
    SAGECAL(SET_NODES.out[0], SET_NODES.out[1])

    //Convert sagecal gains outputs and make diagnostic plots
    ANALYSE_GAINS(SAGECAL.out)

    //Make power spectra using PSPIPE
    POWER_SPECTRUM(SAGECAL.out)

}


workflow VET_PARAMS {
    main:
        nodes_list  = parseNodes()
        params_check_ch = checkingParams()
        checkMSNaming()
        checkDataDistribution(nodes_list)

        //we can now generate the sagecal command
        params.sagecal_command=make_sagecal_command()

        //All params are available now
        savingParams(params_check_ch.all_params_valid)

        //only show params info after validating all params
        if (params.verbose){
            paramsInfo()
        }

    emit:
        savingParams.out.params_json_written
}


workflow SET_NODES {

    take:
        params_json_written

    main:
        if (params.mode == "mpi") {
            mpirun_ch = mpirunNodesList(params_json_written, nodes_list, params.slots, params.mpirun_hosts_txt_file)
            psshNodesList(mpirun_ch.mpi_standby, nodes_list, params.pssh_hosts_txt_file)
        }
        else{
            psshNodesList(params_json_written, nodes_list, params.pssh_hosts_txt_file)
        }

    emit:
        psshNodesList.out.all_nodes_standby
        psshNodesList.out.pssh_hosts_txt
}



workflow SAGECAL {
    take:
        all_nodes_standby
        pssh_hosts_txt

    main:
        //Run the necessary sagecal mode
        if (params.mode == "standalone") {

            sagecalStandalone(all_nodes_standby, params.sagecal_command, pssh_hosts_txt, params.datapath, params.standalone_sage_nf_file, params.msfiles)

            sage_ch = sagecalStandalone.out
        }

        else if (params.mode == "mpi") {

            cp_ch = getModels(all_nodes_standby, pssh_hosts_txt, params.sky_model, params.clusters_file, params.admm_rho_file, params.shapelets_modes, params.datapath)

            preprocess_ch = preProcess(cp_ch, pssh_hosts_txt, params.preprocessing_file, params.datapath, params.msfiles)

            sagecalMPI(preprocess_ch.collect(), params.sagecal_command, params.datapath)

            sage_ch = sagecalMPI.out
        }

        //I needed a way to emit this workflows output regardless of the sagecal mode. This dummy process helps with that.
        sageDone(sage_ch)

    emit:
        sageDone.out

}

workflow ANALYSE_GAINS {
    
    take:
        sagecal_complete

    main:
    
        eff_ch = makeEffNr(sagecal_complete, params.clusters_file)

        conv_ch = convertSageSolutions(eff_ch, params.obsid, params.ms_pattern, params.nodes, params.pssh_hosts_txt_file, params.solsdir, params.datapath, params.stage_number)

        if (params.mode=="mpi"){
            //convert global solutions
            convertSageZSol(conv_ch.ready, eff_ch, params.obsid, params.solsdir, params.datapath)

            if (params.stage == "DI" ) {
                plotDISageSols(conv_ch.npy, conv_ch.npz, eff_ch, params.obsid, params.solsdir, params.plotgains.fmin, params.plotgains.fmax)
            }

            else{
                plotDDSageSols(conv_ch.npy, conv_ch.npz, eff_ch, params.obsid, params.solsdir, params.plotgains.fmin, params.plotgains.fmax)
            }
        }
        else {
            plotDISageSols(conv_ch.npy, conv_ch.npz, eff_ch, params.obsid, params.solsdir, params.plotgains.fmin, params.plotgains.fmax)
        }
}


workflow POWER_SPECTRUM {
    take:
        sagecal_complete

    main:
        init_ch = initPSDB(sagecal_complete, params.datapath)
        rev_ch = addRevision(init_ch.ps_dir, init_ch.default_toml_file, nodes_list[-1], params.max_concurrent, params.revision ) // TODO: stop using only the final node 
        ps_ch = runPSPIPE(init_ch.ps_dir, rev_ch, params.obsid, params.all_msets_txt_list)
        plotPS(ps_ch.ready, init_ch.ps_dir, rev_ch, params.obsid)
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Processes
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

process checkingParams {

    output:
    val true, emit: all_params_valid

    exec:
    checkParams()
}


process savingParams {
    debug true
    publishDir params.outdir

    input:
    val all_params_valid

    output:
        path 'params.json'
        val true, emit: params_json_written

    script:
    """
        echo '${JsonOutput.prettyPrint(JsonOutput.toJson(params))}' > params.json
    """
}


process mpirunNodesList {
    input:
    val ready
    val nodes
    val slots
    val mpirun_hosts_txt_file

    output:
    val "${mpirun_hosts_txt_file}", emit: mpirun_hosts_txt
    val true, emit: mpi_standby

    exec:

    write_nodes_for_mpirun(nodes, slots, mpirun_hosts_txt_file)
}

process psshNodesList {
    input:
    val ready
    val nodes
    val pssh_hosts_txt_file

    output:
    val "${pssh_hosts_txt_file}", emit:  pssh_hosts_txt
    val true, emit: all_nodes_standby

    // seq  -f "node%g" ${first_node} ${last_node} > pssh_hosts_txt_file.txt
    exec:

    write_nodes_for_pssh(nodes, pssh_hosts_txt_file)
}


//cp_ch = getModels(params.pssh_hosts_txt_file, params.sky_model, params.clusters_file, params.admm_rho_file, params.shapelets_modes, params.datapath)
// copying the models fail when i provide them as path
process getModels {
    input:
    val ready
    path pssh_hosts_txt_file
    val sky_model
    val clusters_file
    val admm_rho_file
    val shapelets_modes
    val datapath

    output:
    val true

    script:

    """
    pssh -v -i -h ${pssh_hosts_txt_file} -t 0 cp ${sky_model} ${clusters_file} ${admm_rho_file} ${shapelets_modes} ${datapath} > ${launchDir}/models_copying.log 2>&1
    """
}

process preProcess {
    input:
    val ready
    val pssh_hosts_txt_file
    val preprocessing_file
    val datapath
    val msfiles

    output:
    val true

    script:
    """
    pssh -v -i -h ${pssh_hosts_txt_file} -t 0 -x "cd ${datapath}; bash" ~/mysoftware/nextflow run ${preprocessing_file} --msfiles ${msfiles} > ${launchDir}/preprocessing.log 2>&1
    """
}

// sagecalStandalone(true, params.sagecal_command, params.pssh_hosts_txt_file, params.datapath, params.standalone_sage_nf_file, params.msfiles)
process sagecalStandalone {
    input:
    val ready
    val command
    path pssh_hosts_txt_file
    val datapath
    val standalone_sage_nf_file
    val msfiles

    output:
    val true , emit: sagecal_complete

    script:
    """
    pssh -v -i -h ${pssh_hosts_txt_file} -t 0 -x "cd ${datapath}; bash" ~/mysoftware/nextflow run ${standalone_sage_nf_file} --msfiles ${msfiles} --command "'${command}'"  > ${launchDir}/sagecal_standalone.log 2>&1
    """
}

// sagecalMPI(preprocess_ch, params.sagecal_command, params.datapath)
process sagecalMPI {
    debug true
    cpus params.number_of_threads
    publishDir params.outdir

    input:
    val ready
    val sagecal_command
    val datapath

    output:
    val true , emit: sagecal_complete
    //path(consensus_solutions_file), path(logfile), 

    script:

    """
    cd ${datapath}
    ${sagecal_command}
    """

}

process sageDone {
    input:
    val ready

    output:
    val true , emit: sagecal_complete

    script:

    """
    echo "this is a hack to help the 'SAGECAL' workflow emit a 'sagecal_complete 'signal to the following wokflows"
    """
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Closures
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def checkParams() {
    //check obsid and nodes
    assert params.datapath : log.error("Error:'--datapath' needed")
    assert params.nodes : log.error("Error: Please provide a comma separated string of nodes. Try --nodes '100,101,102,103,104'")
    assert params.obsid : log.error("Error: Please provide an obsid for proper namings. Try --obsid 'your_obsid'")
    assert ["standalone", "mpi"].contains(params.mode) : log.error("Error: Unknown sagecal mode: ${params.mode}. Can only be either 'mpi' or 'standalone'")

    if (params.mode == "mpi"){
        assert params.stage : log.error("Error: In MPI mode, a calibration stage must be chosen. Can be either --stage 'DI' or --stage 'DD'")
        assert ["DI", "DD"].contains(params.stage) : log.error("Error: Unknown lofar calibration stage: ${params.stage}. Can only be either 'DI' or 'DD'")
    }
}


//return a list of the nodes with the 'node' prefix given an input string e.g. [node129]
def parseNodes(){
    if (params.nodes instanceof String){
        nodes_list = params.nodes.split(',').collect{"node${it}"} as List
    }
    //when a single node is given..
    else if (params.nodes instanceof Integer) {
        nodes_list = params.nodes.collect{"node${it}"} as List
    }
    //Define 2 more required params
    params.masternode = nodes_list[-1] //The last node in the list is used as the masternode
    params.number_of_processes = params.slots * nodes_list.size

    return nodes_list
}


//Default LOFAR data has a naming like 'OBSID_SAPNumber_SBNumber_uv_label.MS'
//we check whether this is the convention otherwise use any ms_pattern provided by the user
def checkMSNaming() {

    if (params.obsid && params.stage_number && params.label) {
        //when given through the CLI, all string numbers tend to be converted to a single integer e.g. "003" becomes 3 causing a path error.
        //here I revert it to string and pad it with 2 zeros again
        if (params.stage_number instanceof Integer) {
            stage_number=String.format("%03d", params.stage_number)
        }
        else {
            stage_number = params.stage_number
        }

        params.ms_pattern = "${params.obsid}_SAP000_SB???_uv_${stage_number}_${params.label}.MS" //default LOFAR data MS naming

        calibration_outputs_label = "${params.obsid}_${params.stage}_${stage_number}_${params.label}"
        params.solsdir = "${params.datapath}/solutions_${calibration_outputs_label}"
        params.logfile = "sagecal_${calibration_outputs_label}.log" // sagecal output log file
    }
    
    else {
        println("${params.obsid}, ${params.stage_number}, ${params.label}. Not following default naming convention")
        assert params.ms_pattern : log.error("'--ms_pattern' needed") //MS pattern to search for
    }
    return params.ms_pattern
}

//Given a list_of_strings, write out a .txt file with each string per line 
def writeListToTxt (txt_file_name, list_of_strings) {
        File txt_file = new File(txt_file_name)
        txt_file.withWriter{ out ->
        list_of_strings.each {out.println it}
    }
}

//Given a directory path and a glob pattern, check that the directory exists and has at least one file corresonding to the glob pattern
def checkFilesExistence (files_path, files_glob_pattern) {
    //First check that the directory exists
    dir = file(files_path, type: 'dir', glob: true, checkIfExists: true)
    assert dir.isDirectory() : log.error("Error: datapath param is not a valid path on ${node}. Got '--datapath=${datapath}'")

    //check that we have a least one file with the given glob pattern
    ms_files_path = "${files_path}/${files_glob_pattern}"
    list_of_msets = file(ms_files_path, type: 'dir', glob: true, checkIfExists: true)
    assert list_of_msets.size()>0 : log.error("Error: No measurement sets matching pattern '${ms_files_path}' found.")

    return list_of_msets
}


def checkDataDistribution(nodes_list) {
    // make sure the datapath directory and MS files exist on each node
    // write a list of ms files in each node's datapath

    params.ms_files_per_node  = [:] // A map to populate the ms files in each node
    params.sub_bands_per_node = [:] //A map to populate the subband numbers of the files in each node

    def isSubBand = { it.toString().split("_")[-4].split("SB")[1].toInteger()}
    def isMS = { it.toString().split("/")[-1] }

    all_msets = []
    for (node in nodes_list) {
        fullpath="/net/${node}${params.datapath}"
        list_of_msets = checkFilesExistence(fullpath, params.ms_pattern)
        all_msets.addAll(list_of_msets)

        params.ms_files_per_node."${node}" = list_of_msets.collect(isMS)
        params.sub_bands_per_node."${node}" = list_of_msets.collect(isSubBand)

        //write out a txt file per node with the msets in that node 
        msets_txt_file = "${fullpath}/${params.msfiles}"
        writeListToTxt (msets_txt_file, list_of_msets)

    }

    //write a single file listing all msets from all nodes
    writeListToTxt (params.all_msets_txt_list, all_msets)

    if (params.verbose){
        log.info ("""
                ________________________________________________________________________________
                Info:  MS files
            """.stripIndent()
        )

        params.ms_files_per_node.each { entry -> println "$entry.key : $entry.value.size MS files : $entry.value" }
    }
    return
}


//print out the parameters for the user
def paramsInfo(){
    params_string = """
        ________________________________________________________________________________
        Info: Modes, Nodes and data

        Sagecal Mode                ${params.mode}
        Cluster nodes               ${nodes_list}
        Master node                 ${params.masternode}
        Slots per node              ${params.slots}
        MS pattern                  ${params.ms_pattern}
        NProcesses                  ${params.number_of_processes}
        Datapath                    ${params.datapath}
        Redshift                    ${params.redshift}
        Stage_number                ${params.stage_number}
        ________________________________________________________________________________
        Info:  Models*
        
        Sky Model:                  ${params.sky_model}
        Clusters file:              ${params.clusters_file}
        Shapelets modes files:      ${params.shapelets_modes}

        *copied to ${params.datapath} on each node.

        ________________________________________________________________________________
        Info: Sagecal parameters
        
        MS pattern              -f      ${params.ms_pattern}
        Sky model               -s      ${params.sky_model}                        
        Sky model format        -F      ${params.sky_model_format}
        Clusters file           -c      ${params.clusters_file}
        Beam_model              -B      ${params.beam_model} 
        Input_column            -I      ${params.input_column}
        Output column           -O      ${params.output_column}
        Global solutions file   -p      ${params.solutions_file} (if in MPI mode else ${params.ms_pattern}.solutions)
        Time sampers per sol    -t      ${params.time_samples_per_sol}
        Min lambda cut          -x      ${params.min_lambda}
        Max lambda cut          -y      ${params.max_lambda}
        # of threads            -n      ${params.number_of_threads}
        Max EM iter             -e      ${params.max_em_iterations}
        Single max EM iter      -g      ${params.single_max_em_iterations}
        Max LBFGS iter          -l      ${params.max_lbfgs_iterations}
        LBFGS memory size       -m      ${params.lbfgs_memory_size}
        Solver                  -j      ${params.solver}
        Robust nu lower         -L      ${params.robust_nu_lower}
        Robust nu upper         -H      ${params.robust_nu_upper}
        Pre-whiten              -W      ${params.pre_whiten}
        Apply global solutions  -U      ${params.apply_global_solution}
        Use GPU                 -E      ${params.use_gpu}
        Verbose                 -V      ${params.verbose}
    """.stripIndent()

    if (params.mode == "mpi") {
        log.info("""
                ________________________________________________________________________________
                Info: Running default LOFAR ${params.stage} calibration stage i.e. sagecal mpi with consensus optimization, high regularisation.""".stripIndent())
        consensus_params = """
                ________________________________________________________________________________
                Info: Spectral regularisation params
                
                ADMM iterations         -A      ${params.admm_iterations}
                Consensus Bpol terms    -P      ${params.consensus_polynomial_terms}
                ADMM rho file           -G      ${params.admm_rho_file}
                # of sols               -T      ${params.timechunks}
        """.stripIndent()

        full_params_string = params_string + consensus_params
    }

    else if (params.mode == "standalone") {
        log.info("""
            ________________________________________________________________________________
            Info: Running default LOFAR bandpass calibration (No regularisation)""".stripIndent())
        full_params_string = params_string
    }

    sagecal_string = """
                ________________________________________________________________________________
                Info: Sagecal command to be run
                
                ${params.sagecal_command}
                ________________________________________________________________________________
                """.stripIndent()

    log.info(full_params_string + sagecal_string)

    return
}


def make_sagecal_command() {

            String sage_command =  """
-s ${params.sky_model} \
-F ${params.sky_model_format} \
-c ${params.clusters_file} \
-B ${params.beam_model} \
-I ${params.input_column} \
-O ${params.output_column} \
-t ${params.time_samples_per_sol} \
-x ${params.min_lambda} \
-y ${params.max_lambda} \
-n ${params.number_of_threads} \
-e ${params.max_em_iterations} \
-g ${params.single_max_em_iterations} \
-l ${params.max_lbfgs_iterations} \
-m ${params.lbfgs_memory_size} \
-j ${params.solver} \
-L ${params.robust_nu_lower} \
-H ${params.robust_nu_upper} \
-W ${params.pre_whiten} \
-U ${params.apply_global_solution} \
-E ${params.use_gpu}
""".stripIndent()

    if (params.mode == "mpi") {
        if (params.timechunks) {
            sage_command =  sage_command.strip() + " -T ${params.timechunks}"
        }

        sage_command =  sage_command.strip() + " -A ${params.admm_iterations} -P ${params.consensus_polynomial_terms} -G ${params.admm_rho_file} -p ${params.solutions_file} -f \"${params.ms_pattern}\" -V > ${params.logfile} 2>&1"

        mpi_command = make_mpi_command()

        command = mpi_command + " " +  sage_command

    }

    else if (params.mode == "standalone"){
        command = "sagecal_gpu " + sage_command.strip() + " -k ${params.cluster_id} -a ${params.action} "
    }
    return command.stripIndent()
}


//write the mpi run command
def make_mpi_command() {
    return "${params.mpi} -np ${params.number_of_processes} -hostfile ${params.mpirun_hosts_txt_file} --map-by slot sagecal-mpi_gpu"
}


//mpi machines.txt file
def write_nodes_for_mpirun(nodes_list, slots, mpirun_hosts_txt_file) {
    mpirun_hosts_txt_file = new File(mpirun_hosts_txt_file)
    if (mpirun_hosts_txt_file.exists()){
        mpirun_hosts_txt_file.delete()
    }
    mpirun_hosts_txt_file.createNewFile()
    mpirun_hosts_txt_file.withWriter { out ->
    nodes_list.each {
      out.println("${it}.fast.dawn.rug.nl slots=$slots max_slots=$slots")
    }
  }
}

return

//pssh hosts_list.txt file   
def write_nodes_for_pssh(nodes_list, pssh_hosts_txt_file) {
    pssh_hosts_txt_file = new File(pssh_hosts_txt_file)
    if (pssh_hosts_txt_file.exists()){
        pssh_hosts_txt_file.delete()
    }
    pssh_hosts_txt_file.createNewFile()
    pssh_hosts_txt_file.withWriter { out ->
    nodes_list.each {
      out.println("${it}") // node129
    }
  }
}


//error handler
workflow.onError {
    log.info("Error. Pipeline execution stopped with the following message: ${workflow.errorMessage}")
}

//completion handler
workflow.onComplete {
	log.info ( workflow.success ? '\nDone!' : '\nOops .. something went wrong' )
}



//Earlier I had the separate sagecal modes in these separate workflows.
//I then combined them into the one SAGECAL workflow above h=which is neater to call in the overall workflow.
//Might delete later
// //Run the necessary sagecal mode
// if (params.mode == "standalone") {

//    STANDALONE_SAGECAL(SET_NODES.out[0], SET_NODES.out[1])

// }

// else if (params.mode == "mpi") {

//   MPI_SAGECAL(SET_NODES.out[0], SET_NODES.out[1])

// }

// workflow STANDALONE_SAGECAL {
//     take:
//         all_nodes_standby
//         pssh_hosts_txt

//     main:
//         sagecalStandalone(all_nodes_standby, params.sagecal_command, pssh_hosts_txt, params.datapath, params.standalone_sage_nf_file, params.msfiles)
// }


// workflow MPI_SAGECAL {
//     take:
//         all_nodes_standby
//         pssh_hosts_txt

//     main:
//         cp_ch = getModels(all_nodes_standby, pssh_hosts_txt, params.sky_model, params.clusters_file, params.admm_rho_file, params.shapelets_modes, params.datapath)

//         preprocess_ch = preProcess(cp_ch, pssh_hosts_txt, params.preprocessing_file, params.datapath, params.msfiles)

//         sagecalMPI(preprocess_ch.collect(), params.sagecal_command, params.datapath)
    

//     emit:
//         sagecalMPI.out.sagecal_complete

// }