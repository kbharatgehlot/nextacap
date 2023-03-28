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

params.tasks = "${workflow.profile}".split(",")

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
    
    //flag and average to 002
    // FLAG_GEN_002_VIS(SET_NODES.out[0], SET_NODES.out[1])


    SAGECAL_MPI_DI(SET_NODES.out[0])


    // ANALYSE_GAINS(SAGECAL_MPI_DI.out, "mpi_di")


    SAGECAL_BANDPASS(SAGECAL_MPI_DI.out)


    // ANALYSE_GAINS(SAGECAL_BANDPASS.out, "bandpass")

    //average to 003
    GEN_003_VIS(SAGECAL_BANDPASS.out)

    
    SAGECAL_MPI_DD(GEN_003_VIS.out)


    ANALYSE_GAINS(SAGECAL_MPI_DD.out, "mpi_dd")


    POWER_SPECTRUM(SAGECAL_MPI_DD.out)

}


workflow init {
    VET_PARAMS()
}


workflow VET_PARAMS {
    main:
        nodes_list  = parseNodes()
        params_check_ch = checkingParams()

        checkInitialDataDistribution(nodes_list)

        writeMSListPerNodeStageNumber()

        //we can now generate the sagecal command for each mode
        generateSagecalCommand()
        

        //All params are available now
        savingParams(params_check_ch.all_params_valid)

        //only show params info after validating all params //TODO:
        // if (params.verbose){
        //     paramsInfo()
        // }

    emit:
        savingParams.out.params_json_written
}


workflow SET_NODES {

    take:
        params_json_written

    main:
        if (params.tasks.contains("mpi_di") || params.tasks.contains("mpi_dd")) {
            mpirun_ch = mpirunNodesList(params_json_written, nodes_list, params.cluster.slots, params.cluster.mpirun_hosts_txt_file)
            psshNodesList(mpirun_ch.mpi_standby, nodes_list, params.cluster.pssh_hosts_txt_file)
        }
        else if (params.tasks.contains("bandpass")) {
            psshNodesList(params_json_written, nodes_list, params.cluster.pssh_hosts_txt_file)
        }

    emit:
        psshNodesList.out.all_nodes_standby
        psshNodesList.out.pssh_hosts_txt
}


workflow FLAG_GEN_002_VIS  {
    take:
        all_nodes_standby
        pssh_hosts_txt

    main:
        flagAverageVisTo002(all_nodes_standby, pssh_hosts_txt, params.gen_002_vis.nf_module, params.data.path, params.data.ms_files_001)
    
    emit:
        flagAverageVisTo002.out

}


workflow SAGECAL_MPI_DI {
    take:
        gen_002_ready

     main:
        cp_ch = getModels(gen_002_ready, params.cluster.pssh_hosts_txt_file, params.shapelets.modes, params.data.path)
        di_preprocess_ch = preProcess_di(cp_ch.collect(), params.cluster.pssh_hosts_txt_file, params.mpi_di.preprocessing_file, params.data.path, params.data.ms_files_002)
        sage_mpi_di_ch = sagecalMPI(di_preprocess_ch.collect(), params.mpi_di.sagecal_command, params.data.path)

        eff_ch = makeEffNr(sage_mpi_di_ch.sagecal_complete, params.mpi_di.clusters_file)

        conv_ch = convertSageSolutions(eff_ch, params.data.obsid, params.mpi_di.ms_pattern, params.cluster.nodes, params.cluster.pssh_hosts_txt_file, params.mpi_di.solsdir, params.data.path, params.mpi_di.stage_number)

        convertSageZSol(conv_ch.ready, eff_ch, params.data.obsid, params.mpi_di.solsdir, params.data.path)
        plotDISageSols(conv_ch.npy, conv_ch.npz, eff_ch, params.data.obsid, params.mpi_di.solsdir, params.gains.fmin, params.gains.fmax)
    emit:
        // plotDISageSols.out
        true
}


workflow SAGECAL_BANDPASS {
    take:
        sage_mpi_di_done

    main:
        bandpass_ch = sagecalStandalone(sage_mpi_di_done, params.bandpass.sagecal_command, params.cluster.pssh_hosts_txt_file, params.data.path, params.bandpass.nf_module, params.data.ms_files_002)

        eff_ch = makeEffNr(bandpass_ch.sagecal_complete, params.bandpass.clusters_file)

        conv_ch = convertSageSolutions(eff_ch, params.data.obsid, params.bandpass.ms_pattern, params.cluster.nodes, params.cluster.pssh_hosts_txt_file, params.bandpass.solsdir, params.data.path, params.bandpass.stage_number)

        plotDISageSols(conv_ch.npy, conv_ch.npz, eff_ch, params.data.obsid, params.bandpass.solsdir, params.gains.fmin, params.gains.fmax)

    emit:
        // plotDISageSols.out
        true
}


workflow GEN_003_VIS  {
    take:
        sagecal_bandpass_done

    main:
        averageVisTo003(sagecal_bandpass_done, params.cluster.pssh_hosts_txt_file, params.gen_003_vis.nf_module, params.data.path, params.data.ms_files_002)

    emit:
        averageVisTo003.out
}


workflow SAGECAL_MPI_DD {
    take:
        gen_003_complete

     main:
        // cp_ch = getModels(all_nodes_standby, params.cluster.pssh_hosts_txt_file, params.mpi_dd.sky_model, params.mpi_dd.clusters_file, params.mpi_dd.admm_rho_file, params.shapelets.modes, params.data.path)
        dd_preprocess_ch = preProcess_dd(gen_003_complete, params.cluster.pssh_hosts_txt_file, params.mpi_dd.preprocessing_file, params.data.path, params.data.ms_files_003)
        sagecalMPI(dd_preprocess_ch.collect(), params.mpi_dd.sagecal_command, params.data.path)
    emit:
        sagecalMPI.out
}


workflow ANALYSE_GAINS {
    
    take:
        sagecal_mpi_dd_complete
        gains_type

    main:
        eff_ch = makeEffNr(sagecal_mpi_dd_complete, params."${gains_type}".clusters_file)

        conv_ch = convertSageSolutions(eff_ch, params.data.obsid, params."${gains_type}".ms_pattern, params.cluster.nodes, params.cluster.pssh_hosts_txt_file, params."${gains_type}".solsdir, params.data.path, params."${gains_type}".stage_number)


        if (gains_type=="mpi_di") {
            convertSageZSol(conv_ch.ready, eff_ch, params.data.obsid, params.mpi_di.solsdir, params.data.path)
            plotDISageSols(conv_ch.npy, conv_ch.npz, eff_ch, params.data.obsid, params.mpi_di.solsdir, params.gains.fmin, params.gains.fmax)
        }
        else if (gains_type=="mpi_dd") {
                convertSageZSol(conv_ch.ready, eff_ch, params.data.obsid, params.mpi_dd.solsdir, params.data.path)
                plotDDSageSols(conv_ch.npy, conv_ch.npz, eff_ch, params.data.obsid, params.mpi_dd.solsdir, params.gains.fmin, params.gains.fmax)
            }

        else if (gains_type=="bandpass") {
                plotDISageSols(conv_ch.npy, conv_ch.npz, eff_ch, params.data.obsid, params.bandpass.solsdir, params.gains.fmin, params.gains.fmax)
            }
}


workflow POWER_SPECTRUM {
    take:
        sagecal_complete

    main:
        allMsetsPerStageNumber("003", params.data.sub_bands_per_node)
        init_ch = initPSDB(sagecal_complete, params.data.path)
        rev_ch = addRevision(init_ch.ps_dir, init_ch.default_toml_file, nodes_list[-1], params.pspipe.max_concurrent, params.pspipe.revision ) // TODO: stop using only the final node 
        ps_ch = runPSPIPE(init_ch.ps_dir, rev_ch, params.data.obsid, params.data.all_ms_files_003)
        plotPS(ps_ch.ready, init_ch.ps_dir, rev_ch, params.data.obsid)
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
    publishDir params.logs_dir, mode: 'move'

    input:
    val all_params_valid

    output:
        file 'params.json'
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
    // debug true

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

//cp_ch = getModels(params.cluster.pssh_hosts_txt_file, params.sky_model, params.clusters_file, params.admm_rho_file, params.shapelets_modes, params.data.path)
// copying the models fail when i provide them as path
process getModels {
    debug true

    input:
    val ready
    path pssh_hosts_txt_file
    val shapelets_modes
    val datapath

    output:
    val true

    script:

    """
    pssh -v -i -h ${pssh_hosts_txt_file} -t 0 cp ${shapelets_modes} ${datapath} > ${params.logs_dir}/models_copying.log 2>&1
    """
}


process preProcess_di {
    debug true

    input:
    val ready
    val pssh_hosts_txt_file
    val preprocessing_file
    val datapath
    val ms_files

    output:
    val true

    script:
    """
    pssh -v -i -h ${pssh_hosts_txt_file} -t 0 -x "cd ${datapath}; bash" ~/mysoftware/nextflow run ${preprocessing_file} --ms_files ${ms_files} > ${params.logs_dir}/preprocessing_di.log 2>&1
    """
}


process preProcess_dd {
    debug true

    input:
    val ready
    val pssh_hosts_txt_file
    val preprocessing_file
    val datapath
    val ms_files

    output:
    val true

    script:
    """
    pssh -v -i -h ${pssh_hosts_txt_file} -t 0 -x "cd ${datapath}; bash" ~/mysoftware/nextflow run ${preprocessing_file} --ms_files ${ms_files} > ${params.logs_dir}/preprocessing_dd.log 2>&1
    """
}

process flagAverageVisTo002 {
    debug true

    input:
    val ready
    val pssh_hosts_txt_file
    val dp3_average_to_002_file
    val datapath
    val ms_files

    output:
    val true

    script:
    """
    pssh -v -i -h ${pssh_hosts_txt_file} -t 0 -x "cd ${datapath}; bash" ~/mysoftware/nextflow run ${dp3_average_to_002_file} --ms_files ${ms_files} > ${params.logs_dir}/flag_average_to_002.log 2>&1
    """
}


process averageVisTo003 {
    debug true

    input:
    val ready
    val pssh_hosts_txt_file
    val dp3_average_to_003_file
    val datapath
    val ms_files

    output:
    val true

    script:
    //--ms_files ${ms_files}
    """
    pssh -v -i -h ${pssh_hosts_txt_file} -t 0 -x "cd ${datapath}; bash" ~/mysoftware/nextflow run ${dp3_average_to_003_file} --path ${datapath} > ${params.logs_dir}/average_to_003.log 2>&1
    """
}


// sagecalStandalone(true, params.sagecal_command, params.cluster.pssh_hosts_txt_file, params.data.path, params.standalone_sage_nf_file, params.ms_files)
process sagecalStandalone {
    debug true
    // errorStrategy 'ignore'

    input:
    val ready
    val command
    path pssh_hosts_txt_file
    val datapath
    val standalone_sage_nf_file
    val ms_files

    output:
    val true , emit: sagecal_complete

    script:
    """
    #echo "done"
    pssh -v -i -h ${pssh_hosts_txt_file} -t 0 -x "cd ${datapath}; bash" ~/mysoftware/nextflow run ${standalone_sage_nf_file} --ms_files ${ms_files} --command "'${command}'"  > ${params.logs_dir}/sagecal_standalone.log 2>&1
    """
}


process sagecalMPI {
    debug true
    cpus params.mpi_di.number_of_threads

    input:
    val ready
    val sagecal_command
    val datapath

    output:
    val true , emit: sagecal_complete
    //path(consensus_solutions_file), path(logfile),

    script:
    //when done we remove the many shapelets .modes files
    """
    echo "done"
    #cd ${datapath}
    #${sagecal_command}
    #rm *.modes 
    """
    //the rm modes should give the absolute path
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Closures
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def checkParams() {
    //check obsid and nodes
    assert params.data.path : log.error("Error:'--data.path' needed")
    assert params.cluster.nodes : log.error("Error: Please provide a comma separated string of nodes. Try --cluster.nodes '100,101,102,103,104'")
}

/*
Generate a sagecal command fo rthe different LOFAr EoR calibration stages
*/
def generateSagecalCommand() {
    if (params.tasks.contains("mpi_di")){
        params.mpi_di.ms_pattern = lofarDefaultMsnamePattern(params.mpi_di.stage_number)
        params.mpi_di.sagecal_command =  sagecalMPIDICommand()
    }
    if (params.tasks.contains("bandpass")){
        params.bandpass.ms_pattern = lofarDefaultMsnamePattern(params.bandpass.stage_number)
        params.bandpass.sagecal_command =  sagecalBandpassCommand()
    }
    if (params.tasks.contains("mpi_dd")){
        params.mpi_dd.ms_pattern = lofarDefaultMsnamePattern(params.mpi_dd.stage_number)
        params.mpi_dd.sagecal_command =  sagecalMPIDDCommand()
    }
}


/*
return a list of the nodes with the 'node' prefix given an input string e.g. [node129]
*/
def parseNodes(){
    if (params.cluster.nodes instanceof String){
        nodes_list = params.cluster.nodes.split(',').collect{"node${it}"} as List
    }
    //when a single node is given..
    else if (params.cluster.nodes instanceof Integer) {
        nodes_list = params.cluster.nodes.collect{"node${it}"} as List
    }
    //Define 2 more required params
    params.cluster.masternode = nodes_list[-1] //The last node in the list is used as the masternode
    params.cluster.number_of_processes = params.cluster.slots * nodes_list.size

    return nodes_list
}


/*
mpirun needs a file listing the nodes to run commands on and the number of slots to use per node.
This command makes sucha a file given a nodeslist, slots and a file path.
input:      list of strings e.g. ["node100", "node101"]
            slots [int]
            file name string e.g "hosts_list.txt"
output:     The written file
*/
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


/*
pssh needs a file listing the nodes to run commands on
This command makes sucha a file given a nodeslist and a file path.
input:      list of strings e.g. ["node100", "node101"]
            file name string e.g "hosts_list.txt"
output:     The written file
*/
def write_nodes_for_pssh(nodes_list, pssh_hosts_txt_file) {
    pssh_hosts_txt_file = new File(pssh_hosts_txt_file)
    if (pssh_hosts_txt_file.exists()){
        pssh_hosts_txt_file.delete()
    }
    pssh_hosts_txt_file.createNewFile()
    pssh_hosts_txt_file.withWriter { out ->
    nodes_list.each {
      out.println("${it}")
    }
  }
}


/*
Given a directory path and a glob pattern,
Check that the directory exists and has at least one file corresonding to the glob pattern
*/
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

/*
Given the tasks requested using 'params.task' return a map of each task with its respective 'stage_number'
This is forthe conventional MS naming
*/
def getRequiredtasksStageNumbers(){
    def tasks_number_map = ["gen002vis": 1, "gen003vis": 2, "mpi_di": 2, "bandpass": 2, "mpi_dd": 3  ]
    def requested_tasks = [:]
    for (task in params.tasks){
        if (tasks_number_map.containsKey(task)){
            requested_tasks."${task}" = tasks_number_map."${task}"
        }
    }
    return requested_tasks
}

/*
We want to know the inital stage number to expect for the first set of measurement sets
Can either be 001, 003, 0r 003.
We use this to check that the MS files exist and the data subbands distribution per node.
*/
def getInitialStageNumber() {
    requested_tasks = getRequiredtasksStageNumbers()
    initial_stage_number = requested_tasks.min { it.value }.value
    initial_stage_number = String.format("%03d", initial_stage_number)
    log.info("Info: initial stage number: ${initial_stage_number}")
    return initial_stage_number 
}


/*
Some proesses may need a list of all measurement sets across all nodes. e.g pspipe
Make such a file given a map containing the subbands per node and the stage number of the data needed
*/
def allMsetsPerStageNumber(stage_number, sub_bands_per_node){
    mses = []
    for (node in sub_bands_per_node) {
        for (sb in node.value) {
            fyl = lofarDefaultMsnamePattern(stage_number).replace("???", "${sb}")
            full_fyl_path = "/net/${node.key}${params.data.path}/${fyl}"
            mses.add(full_fyl_path)
        }
    }

    txt = "/net/${params.cluster.masternode}/${params.data.path}/all_ms_files_${stage_number}.txt"
    writeListToTxt(txt, mses)
}


def _writeMSListPerNodeStageNumber(stage_number, sub_bands_per_node){
    stage_number = String.format("%03d", stage_number)
    for ( node in sub_bands_per_node ) {
        mses = []
        txt = "/net/${node.key}/${params.data.path}/ms_files_${stage_number}.txt"
        for (sb in node.value) {
            prefix = lofarDefaultMsnamePattern(stage_number).replace("???", "${sb}")
            ms = "/net/${node.key}/${params.data.path}/${prefix}"
            mses.add(ms)
        }
        writeListToTxt(txt, mses)
    }
}


def writeMSListPerNodeStageNumber() {
    requested_tasks = getRequiredtasksStageNumbers()
    stage_num_per_task = requested_tasks.values() as List

    for (n in stage_num_per_task.unique()) {
        _writeMSListPerNodeStageNumber(n, params.data.sub_bands_per_node)
    }
}


/*
based on the tasks, get the minimum stage number
Generate the corresponding LOFAR MS pattern
Check if we have data in the provided nodes following that pattern
save the subbands of the data found in each node
save the msfiles found in each node
write a single file listing all msets from all nodes
*/
def checkInitialDataDistribution(nodes_list) {
    initial_stage_number = getInitialStageNumber()
    initial_ms_pattern = lofarDefaultMsnamePattern(initial_stage_number)
    log.info("Info: initial ms pattern: ${initial_ms_pattern}")

    params.data.sub_bands_per_node = [:]

    def isSubBand = { it.toString().split("_")[-4].split("SB")[1].toInteger()}

    for (node in nodes_list) {
        fullpath="/net/${node}${params.data.path}"
        list_of_msets = checkFilesExistence(fullpath, initial_ms_pattern)
        params.data.sub_bands_per_node."${node}" = list_of_msets.collect(isSubBand)

    }

}

    
/*
Using the LOFAR EoR naming convention,
return the expected MS name pattern at a given analysis stage level
*/
def lofarDefaultMsnamePattern(stage_number) {
    if (stage_number instanceof Integer) {
        stage_number=String.format("%03d", stage_number)
    }
    else{
        assert stage_number instanceof String
        assert stage_number.length() == 3
    }

    ms_pattern = "${params.data.obsid}_SAP${params.data.sub_array_pointing}_SB???_uv_${stage_number}_${params.data.label}.MS"

    return ms_pattern
}


def sagecalBandpassCommand() {
    if (params.bandpass.sagecal_command) {
        log.info("Info: using user provided sagecal bandpass command")
        return params.bandpass.sagecal_command
    }
    else {
        String sage_bandpass_command =  """
-s ${params.bandpass.sky_model} \
-F ${params.bandpass.sky_model_format} \
-c ${params.bandpass.clusters_file} \
-B ${params.bandpass.beam_model} \
-I ${params.bandpass.input_column} \
-O ${params.bandpass.output_column} \
-t ${params.bandpass.time_samples_per_sol} \
-x ${params.bandpass.min_lambda} \
-y ${params.bandpass.max_lambda} \
-n ${params.bandpass.number_of_threads} \
-e ${params.bandpass.max_em_iterations} \
-g ${params.bandpass.single_max_em_iterations} \
-l ${params.bandpass.max_lbfgs_iterations} \
-m ${params.bandpass.lbfgs_memory_size} \
-j ${params.bandpass.solver} \
-L ${params.bandpass.robust_nu_lower} \
-H ${params.bandpass.robust_nu_upper} \
-W ${params.bandpass.pre_whiten} \
-U ${params.bandpass.apply_global_solution} \
-k ${params.bandpass.gains_correction_cluster_id} \
-a ${params.bandpass.action} \
-D ${params.bandpass.enable_diagnostics} \
-E ${params.bandpass.use_gpu}
""".stripIndent()
    
        bandpass_command = "sagecal_gpu " + sage_bandpass_command.strip()
        bandpass_command= bandpass_command.stripIndent()

        return bandpass_command
    }
}

def sagecalMPIDICommand(){
    if (params.mpi_di.sagecal_command) {
        log.info("using user provided sagecal MPI DI command")
        return params.mpi_di.sagecal_command
    }
    else{
            String sage_mpi_di_command =  """ \
-s ${params.mpi_di.sky_model} \
-F ${params.mpi_di.sky_model_format} \
-c ${params.mpi_di.clusters_file} \
-B ${params.mpi_di.beam_model} \
-I ${params.mpi_di.input_column} \
-O ${params.mpi_di.output_column} \
-t ${params.mpi_di.time_samples_per_sol} \
-x ${params.mpi_di.min_lambda} \
-y ${params.mpi_di.max_lambda} \
-n ${params.mpi_di.number_of_threads} \
-e ${params.mpi_di.max_em_iterations} \
-g ${params.mpi_di.single_max_em_iterations} \
-l ${params.mpi_di.max_lbfgs_iterations} \
-m ${params.mpi_di.lbfgs_memory_size} \
-j ${params.mpi_di.solver} \
-L ${params.mpi_di.robust_nu_lower} \
-H ${params.mpi_di.robust_nu_upper} \
-W ${params.mpi_di.pre_whiten} \
-U ${params.mpi_di.apply_global_solution} \
-A ${params.mpi_di.admm_iterations} \
-P ${params.mpi_di.consensus_polynomial_terms} \
-G ${params.mpi_di.admm_rho_file} \
-p ${params.mpi_di.solutions_file} \
-f \"${params.mpi_di.ms_pattern}\" \
-E ${params.mpi_di.use_gpu} \
-V > ${params.mpi_di.logfile}  2>&1 
""".stripIndent()


        mpi_command = make_mpi_command()

        mpi_di_command = mpi_command + " " +  sage_mpi_di_command.strip()
        mpi_di_command= mpi_di_command.stripIndent()

        return mpi_di_command
    }
}


def sagecalMPIDDCommand(){
    if (params.mpi_dd.sagecal_command) {
        log.info("using user provided sagecal MPI DD command")
        return params.mpi_dd.sagecal_command
    }
    else {
            String sage_mpi_dd_command =  """ \
-s ${params.mpi_dd.sky_model} \
-F ${params.mpi_dd.sky_model_format} \
-c ${params.mpi_dd.clusters_file} \
-B ${params.mpi_dd.beam_model} \
-I ${params.mpi_dd.input_column} \
-O ${params.mpi_dd.output_column} \
-t ${params.mpi_dd.time_samples_per_sol} \
-x ${params.mpi_dd.min_lambda} \
-y ${params.mpi_dd.max_lambda} \
-n ${params.mpi_dd.number_of_threads} \
-e ${params.mpi_dd.max_em_iterations} \
-g ${params.mpi_dd.single_max_em_iterations} \
-l ${params.mpi_dd.max_lbfgs_iterations} \
-m ${params.mpi_dd.lbfgs_memory_size} \
-j ${params.mpi_dd.solver} \
-L ${params.mpi_dd.robust_nu_lower} \
-H ${params.mpi_dd.robust_nu_upper} \
-W ${params.mpi_dd.pre_whiten} \
-U ${params.mpi_dd.apply_global_solution} \
-A ${params.mpi_dd.admm_iterations} \
-P ${params.mpi_dd.consensus_polynomial_terms} \
-G ${params.mpi_dd.admm_rho_file} \
-p ${params.mpi_dd.solutions_file} \
-f \"${params.mpi_dd.ms_pattern}\" \
-E ${params.mpi_dd.use_gpu} \
-V > ${params.mpi_dd.logfile}  2>&1 
""".stripIndent()

        mpi_command = make_mpi_command()

        mpi_dd_command = mpi_command.strip() + " " +  sage_mpi_dd_command.strip()
        mpi_dd_command = mpi_dd_command.stripIndent()

        return mpi_dd_command
    }
}

//write the mpi run command
def make_mpi_command() {
    return "${params.cluster.mpi} -np ${params.cluster.number_of_processes} -hostfile ${params.cluster.mpirun_hosts_txt_file} --map-by slot sagecal-mpi_gpu"
}


//Given a list_of_strings, write out a .txt file with each string per line 
def writeListToTxt (txt_file_name, list_of_strings) {
        File txt_file = new File(txt_file_name)
        txt_file.withWriter{ out ->
        list_of_strings.each {out.println it}
    }
}


