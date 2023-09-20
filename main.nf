#!/usr/bin/env nextflow

// Show help message and exit
if (params.help) {
    helpMessage()
    exit 0
}

import groovy.json.JsonOutput
//import workflows
include {
    MakeEffectiveClustersNumberFile;
    ConvertSagecalSolutions;
    ConvertSagecalGlobalSolutions;
    PlotSagecalDISolutions;
    PlotSagecalDDSolutions;
} from './modules/gains_analysis.nf'

include {
    InitPSDB;
    AddRevision;
    RunPSPIPE;
    PlotPowerSpectrum;
} from './modules/power_spectra.nf'

//The tasks to be run are declared using the '-profile' command line option as a comma-separated string
//we split this string into a list.
 //"${workflow.profile}".split(",")
PRODUCTION_PIPELINE="standard,gen002vis,mpi_di,bandpass,gen003vis,mpi_dd,gains,ps,wsclean"
PRODUCTION_PIPELINE_FROM_FIRST_DI="standard,mpi_di,bandpass,gen003vis,mpi_dd,gains,ps,wsclean" //No gen002vis
PRODUCTION_PIPELINE_FROM_BANDPASS="standard,bandpass,gen003vis,mpi_dd,gains,ps,wsclean"
PRODUCTION_PIPELINE_FROM_POST_DI_AVERAGING="standard,gen003vis,mpi_dd,gains,ps,wsclean"
PRODUCTION_PIPELINE_FROM_DD="standard,mpi_dd,gains,ps,wsclean"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                        // INITIALIZATION Workflow //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

workflow Initialise {
    VET_PARAMS()
    SET_NODES(VET_PARAMS.out)
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                        //  MAIN WORKFLOW SELECTOR //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
workflow {

    if ("${workflow.profile}"=="standard") {
        defaultLOFAREoRPipelineMessage()
        params.tasks = PRODUCTION_PIPELINE.split(",")
        Production()
    }
    else {
        requested_tasks_set = "${workflow.profile}".split(",") as Set
        pipeline_set =  PRODUCTION_PIPELINE.split(",") as Set
        pipeline_from_first_di_set =  PRODUCTION_PIPELINE_FROM_FIRST_DI.split(",") as Set
        pipeline_from_bandpass_set =  PRODUCTION_PIPELINE_FROM_BANDPASS.split(",") as Set
        pipeline_from_post_di_averaging_set =  PRODUCTION_PIPELINE_FROM_POST_DI_AVERAGING.split(",") as Set
        pipeline_from_dd_set =  PRODUCTION_PIPELINE_FROM_DD.split(",") as Set

        if (requested_tasks_set == pipeline_set) {
            defaultLOFAREoRPipelineMessage()
            params.tasks = PRODUCTION_PIPELINE.split(",")
            if (params.init) {
                Initialise()
            }
            else{
                Production()
            }
        }
        else if (requested_tasks_set == pipeline_from_first_di_set) {
            defaultLOFArEoRPipelineFrom002Message()
            params.tasks = PRODUCTION_PIPELINE_FROM_FIRST_DI.split(",")
            if (params.init) {
                Initialise()
            }
            else {
                ProductionPipelineFromFirstDI()
            }
        }
        else if (requested_tasks_set == pipeline_from_bandpass_set) {
            defaultLOFArEoRPipelineFromBandpassMessage()
            params.tasks = PRODUCTION_PIPELINE_FROM_BANDPASS.split(",")
            if (params.init) {
                Initialise()
            }
            else {
                ProductionPipelineFromBandpass()
            }
        }
        else if (requested_tasks_set == pipeline_from_post_di_averaging_set) {
            defaultLOFArEoRPipelineFromPostDIAveraging()
            params.tasks = PRODUCTION_PIPELINE_FROM_POST_DI_AVERAGING.split(",")
            if (params.init) {
                Initialise()
            }
            else {
                ProductionPipelineFromPostDIAveraging()
            }
        }
        else if (requested_tasks_set == pipeline_from_dd_set) {
            defaultLOFArEoRPipelineFromDDMessage()
            params.tasks = PRODUCTION_PIPELINE_FROM_DD.split(",")
            if (params.init) {
                Initialise()
            }
            else {
                ProductionPipelineFromDD()
            }
        }
        else {
            unsupportedTasksCombinationError()
            exit 0
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                        //  RUNNING WORKFLOW INFORMATION    //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def unsupportedTasksCombinationError() {
    log.error """
Error: The requested tasks combination is currrently unsupported (`${workflow.profile}`)
Error: Supported task combinations include:\n
    --> ${PRODUCTION_PIPELINE} (full production pipeline, default)\n 
    --> ${PRODUCTION_PIPELINE_FROM_FIRST_DI} (production pipeline from the first DI step)\n
    --> ${PRODUCTION_PIPELINE_FROM_BANDPASS} (production pipeline from the second DI step)\n
    --> ${PRODUCTION_PIPELINE_FROM_DD} (production pipeline from the second DI step)\n
"""
}


def defaultLOFAREoRPipelineMessage() {
    log.info """
[nextleap]: Running default LOFAR EoR KSP production pipeline:
[nextleap]: ${PRODUCTION_PIPELINE}
"""
}

def defaultLOFArEoRPipelineFrom002Message() {
    log.info """
[nextleap]: Running default LOFAR EoR KSP production pipeline from the the first Direction Independent calibration step (fast time-varying effects with spectral regularisation):
[nextleap]: ${PRODUCTION_PIPELINE_FROM_FIRST_DI}
"""
}

def defaultLOFArEoRPipelineFromBandpassMessage() {
    log.info """
[nextleap]: Running default LOFAR EoR KSP production pipeline from the second Direction Independent calibration step (bandpass without spectral regularisation):
[nextleap]: ${PRODUCTION_PIPELINE_FROM_FIRST_DI}
"""
}

def defaultLOFArEoRPipelineFromPostDIAveraging() {
    log.info """
[nextleap]: Running default LOFAR EoR KSP production pipeline from an averaging step towards Direction Dependent Calibration:
[nextleap]: ${PRODUCTION_PIPELINE_FROM_POST_DI_AVERAGING}
"""
}

def defaultLOFArEoRPipelineFromDDMessage() {
    log.info """
[nextleap]: Running default LOFAR EoR KSP production pipeline from the Direction Dependent calibration step:
[nextleap]: ${PRODUCTION_PIPELINE_FROM_DD}
"""
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                        // PRODUCTION Workflow //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

workflow Production {
    //Check that all params are valid
    //save them in a json file for modification if needed and for posterity
    VET_PARAMS()

    //Generate mpirun and pssh hosts files
    SET_NODES(VET_PARAMS.out)
    
    // flag and average to 002
    FLAG_GEN_002_VIS(SET_NODES.out[0], SET_NODES.out[1])

    //write out a txt file with all 002 MS files (A single one located in the master node data path)
    allMsetsPerStageNumber("002", params.data.sub_bands_per_node)

    // run mpi DI
    SAGECAL_MPI_DI(FLAG_GEN_002_VIS.out)

    //run bandpass
    SAGECAL_BANDPASS(SAGECAL_MPI_DI.out)

    // average to 003
    GEN_003_VIS(SAGECAL_BANDPASS.out)

    //write out a txt file with all 003 MS files (A single one located in the master node data path)
    allMsetsPerStageNumber("003", params.data.sub_bands_per_node)

    //Run sagecal mpi DD
    SAGECAL_MPI_DD(GEN_003_VIS.out)

    //Plot sagecal DD gains
    ANALYSE_GAINS(SAGECAL_MPI_DD.out, "mpi_dd")

    //run pspipe
    POWER_SPECTRUM(SAGECAL_MPI_DD.out)

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                // PRODUCTION Workflow starting from the first DI step //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
workflow ProductionPipelineFromFirstDI {
    //Check that all params are valid
    //save them in a json file for modification if needed and for posterity
    VET_PARAMS()

    //Generate mpirun and pssh hosts files
    SET_NODES(VET_PARAMS.out)

    //write out a txt file with all 002 MS files (A single one located in the master node data path)
    allMsetsPerStageNumber("002", params.data.sub_bands_per_node)

    // run mpi DI
    SAGECAL_MPI_DI(SET_NODES.out[0])

    //run bandpass
    SAGECAL_BANDPASS(SAGECAL_MPI_DI.out)

    // average to 003
    GEN_003_VIS(SAGECAL_BANDPASS.out)

    //write out a txt file with all 003 MS files (A single one located in the master node data path)
    allMsetsPerStageNumber("003", params.data.sub_bands_per_node)

    //Run sagecal mpi DD
    SAGECAL_MPI_DD(GEN_003_VIS.out)

    //Plot sagecal DD gains
    ANALYSE_GAINS(SAGECAL_MPI_DD.out, "mpi_dd")

    //run pspipe
    POWER_SPECTRUM(SAGECAL_MPI_DD.out)

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                // PRODUCTION Workflow starting from the second DI step //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
workflow ProductionPipelineFromBandpass {
    //Check that all params are valid
    //save them in a json file for modification if needed and for posterity
    VET_PARAMS()

    //Generate mpirun and pssh hosts files
    SET_NODES(VET_PARAMS.out)

    //write out a txt file with all 002 MS files (A single one located in the master node data path)
    allMsetsPerStageNumber("002", params.data.sub_bands_per_node)


    //run bandpass
    SAGECAL_BANDPASS(SET_NODES.out[0])

    // average to 003
    GEN_003_VIS(SAGECAL_BANDPASS.out)

    //write out a txt file with all 003 MS files (A single one located in the master node data path)
    allMsetsPerStageNumber("003", params.data.sub_bands_per_node)

    //Run sagecal mpi DD
    SAGECAL_MPI_DD(GEN_003_VIS.out)

    //Plot sagecal DD gains
    ANALYSE_GAINS(SAGECAL_MPI_DD.out, "mpi_dd")

    //run pspipe
    POWER_SPECTRUM(SAGECAL_MPI_DD.out)

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        // PRODUCTION Workflow starting from the averaging step after DI calibration //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
workflow ProductionPipelineFromPostDIAveraging {
    //Check that all params are valid
    //save them in a json file for modification if needed and for posterity
    VET_PARAMS()

    //Generate mpirun and pssh hosts files
    SET_NODES(VET_PARAMS.out)

    //write out a txt file with all 002 MS files (A single one located in the master node data path)
    allMsetsPerStageNumber("002", params.data.sub_bands_per_node)

    // average to 003
    GEN_003_VIS(SET_NODES.out[0])

    //write out a txt file with all 003 MS files (A single one located in the master node data path)
    allMsetsPerStageNumber("003", params.data.sub_bands_per_node)

    //Run sagecal mpi DD
    SAGECAL_MPI_DD(GEN_003_VIS.out)

    //Plot sagecal DD gains
    ANALYSE_GAINS(SAGECAL_MPI_DD.out, "mpi_dd")

    //run pspipe
    POWER_SPECTRUM(SAGECAL_MPI_DD.out)

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        // PRODUCTION Workflow starting from DD calibration //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
workflow ProductionPipelineFromDD {
    //Check that all params are valid
    //save them in a json file for modification if needed and for posterity
    VET_PARAMS()

    //Generate mpirun and pssh hosts files
    SET_NODES(VET_PARAMS.out)

    //write out a txt file with all 002 MS files (A single one located in the master node data path)
    // allMsetsPerStageNumber("002", params.data.sub_bands_per_node)

    //write out a txt file with all 003 MS files (A single one located in the master node data path)
    allMsetsPerStageNumber("003", params.data.sub_bands_per_node)

    //Run sagecal mpi DD
    SAGECAL_MPI_DD(SET_NODES.out[0])

    //Plot sagecal DD gains
    ANALYSE_GAINS(SAGECAL_MPI_DD.out, "mpi_dd")

    //run pspipe
    POWER_SPECTRUM(SAGECAL_MPI_DD.out)

}

workflow VET_PARAMS {
    main:
        log.info """[nextleap]: initialising params"""
        params_valid = checkParams() //CheckingParams()

        log.info """[nextleap]: verifying nodes list"""
        nodes_list  = parseNodes()

        log.info """[nextleap]: verifying data distribution"""
        checkInitialDataDistribution(nodes_list)
        writeMSListPerNodeStageNumber()

        //we can now generate the sagecal command for each mode
        generateSagecalCommand()
        

        //All params are available now
        SavingParams(params_valid)

        //only show params info after validating all params //TODO:
        // if (params.verbose){
        //     paramsInfo()
        // }

    emit:
        SavingParams.out.params_json_written
}


workflow SET_NODES {

    take:
        params_json_written

    main:
        if (params.tasks.contains("mpi_di") || params.tasks.contains("mpi_dd")) {
            mpirun_ch = MpiRunNodesList(params_json_written, nodes_list, params.cluster.slots, params.cluster.mpirun_hosts_txt_file)
            PsshNodesList(mpirun_ch.mpi_standby, nodes_list, params.cluster.pssh_hosts_txt_file)
        }
        else { //if (params.tasks.contains("bandpass"))
            PsshNodesList(params_json_written, nodes_list, params.cluster.pssh_hosts_txt_file)
        }

    emit:
        PsshNodesList.out.all_nodes_standby
        PsshNodesList.out.pssh_hosts_txt
}


workflow FLAG_GEN_002_VIS  {
    take:
        all_nodes_standby
        pssh_hosts_txt

    main:
        FlagAndAverageVisTo002(all_nodes_standby, params.cluster.pssh_hosts_txt_file, params.gen_002_vis.nf_module, params.data.path, "${params.data.path}/ms_files_001.txt")
    
    emit:
        FlagAndAverageVisTo002.out

}

workflow SAGECAL_MPI_DI {
    take:
        gen_002_ready

     main:

        cp_ch = GetModels(gen_002_ready, params.cluster.pssh_hosts_txt_file, params.shapelets.modes, params.data.path, params.sim)
        add_cols_ch = AddColumnToMeasurementSet(cp_ch.collect(), params.cluster.pssh_hosts_txt_file, params.data.path, "${params.data.path}/ms_files_002.txt", params.mpi_di.output_column)
        di_preprocess_ch = ApplyPreDIFlag(add_cols_ch.collect(), params.cluster.pssh_hosts_txt_file, params.mpi_di.preprocessing_file, params.data.path, "${params.data.path}/ms_files_002.txt", params.sim)
        sage_mpi_di_ch = SagecalMPI(di_preprocess_ch.collect(), params.mpi_di.sagecal_command, params.data.path, params.cluster.pssh_hosts_txt_file, params.mpi_di.solsdir, params.mpi_di.ms_pattern)

        eff_ch = MakeEffectiveClustersNumberFile(sage_mpi_di_ch.sagecal_complete, params.mpi_di.clusters_file)

        conv_ch = ConvertSagecalSolutions(eff_ch, params.data.obsid, params.mpi_di.ms_pattern, params.cluster.nodes, params.cluster.pssh_hosts_txt_file, params.mpi_di.solsdir, params.data.path, params.mpi_di.stage_number)

        ConvertSagecalGlobalSolutions(conv_ch.ready, eff_ch, params.data.obsid, params.mpi_di.solsdir, params.data.path)
        PlotSagecalDISolutions(conv_ch.npy, conv_ch.npz, eff_ch, params.data.obsid, params.mpi_di.solsdir, params.gains.fmin, params.gains.fmax)
    emit:
        PlotSagecalDISolutions.out
}


workflow SAGECAL_BANDPASS {
    take:
        sage_mpi_di_done

    main:
        add_cols_ch = AddColumnToMeasurementSet(sage_mpi_di_done, params.cluster.pssh_hosts_txt_file, params.data.path, "${params.data.path}/ms_files_002.txt", params.bandpass.output_column)

        bandpass_ch = SagecalStandalone(add_cols_ch.collect(), params.bandpass.sagecal_command, params.cluster.pssh_hosts_txt_file, params.data.path, params.bandpass.nf_module, "${params.data.path}/ms_files_002.txt", params.bandpass.solsdir, params.bandpass.time_limit)

        wsc_ch = ImageWithWSClean(bandpass_ch.sagecal_complete, params.wsclean.scale, params.wsclean.size, params.wsclean.weight, params.wsclean.minuv_lambda, params.wsclean.maxuv_lambda, params.wsclean.polarisation, params.wsclean.threads, params.bandpass.output_column, "${params.wsclean.dir}_DI", "${params.data.path}/all_ms_files_002.txt")

        eff_ch = MakeEffectiveClustersNumberFile(bandpass_ch.sagecal_complete, params.bandpass.clusters_file) //wsc_ch.wsclean_complete

        conv_ch = ConvertSagecalSolutions(eff_ch, params.data.obsid, params.bandpass.ms_pattern, params.cluster.nodes, params.cluster.pssh_hosts_txt_file, params.bandpass.solsdir, params.data.path, params.bandpass.stage_number)

        PlotSagecalDISolutions(conv_ch.npy, conv_ch.npz, eff_ch, params.data.obsid, params.bandpass.solsdir, params.gains.fmin, params.gains.fmax)


    emit:
        PlotSagecalDISolutions.out
}


workflow GEN_003_VIS  {
    take:
        sagecal_bandpass_done

    main:
        AverageVisTo003(sagecal_bandpass_done, params.cluster.pssh_hosts_txt_file, params.gen_003_vis.nf_module, params.data.path, params.gen_003_vis.msin.datacolumn, "${params.data.path}/ms_files_002.txt")

    emit:
        AverageVisTo003.out
}


workflow SAGECAL_MPI_DD {
    take:
        gen_003_complete

     main:
        //  the di run already copied the models so no need to copy again
        cp_ch = GetModels(gen_003_complete, params.cluster.pssh_hosts_txt_file, params.shapelets.modes, params.data.path, params.sim)
        add_cols_ch = AddColumnToMeasurementSet(cp_ch.collect(), params.cluster.pssh_hosts_txt_file, params.data.path, "${params.data.path}/ms_files_003.txt", params.mpi_dd.output_column)
        dd_preprocess_ch = ApplyPreDDFlag(add_cols_ch.collect(), params.cluster.pssh_hosts_txt_file, params.mpi_dd.preprocessing_file, params.data.path, "${params.data.path}/ms_files_003.txt", params.sim)
        SagecalMPI(dd_preprocess_ch.collect(), params.mpi_dd.sagecal_command, params.data.path, params.cluster.pssh_hosts_txt_file, "${params.data.path}/${params.mpi_dd.solsdir}", params.mpi_dd.ms_pattern)
        ImageWithWSClean(SagecalMPI.out, params.wsclean.scale, params.wsclean.size, params.wsclean.weight, params.wsclean.minuv_lambda, params.wsclean.maxuv_lambda, params.wsclean.polarisation, params.wsclean.threads, params.mpi_dd.output_column, "${params.wsclean.dir}_DD", "${params.data.path}/all_ms_files_003.txt")
    emit:
        SagecalMPI.out
}


workflow ANALYSE_GAINS {
    
    take:
        sagecal_mpi_dd_complete
        gains_type

    main:

        eff_ch = MakeEffectiveClustersNumberFile(sagecal_mpi_dd_complete, params."${gains_type}".clusters_file)

        conv_ch = ConvertSagecalSolutions(eff_ch, params.data.obsid, params."${gains_type}".ms_pattern, params.cluster.nodes, params.cluster.pssh_hosts_txt_file, params."${gains_type}".solsdir, params.data.path, params."${gains_type}".stage_number)


        if (gains_type=="mpi_di") {
            ConvertSagecalGlobalSolutions(conv_ch.ready, eff_ch, params.data.obsid, params.mpi_di.solsdir, params.data.path)
            PlotSagecalDISolutions(conv_ch.npy, conv_ch.npz, eff_ch, params.data.obsid, params.mpi_di.solsdir, params.gains.fmin, params.gains.fmax)
        }
        else if (gains_type=="mpi_dd") {
                ConvertSagecalGlobalSolutions(conv_ch.ready, eff_ch, params.data.obsid, "${params.data.path}/${params.mpi_dd.solsdir}", params.data.path)
                PlotSagecalDDSolutions(conv_ch.npy, conv_ch.npz, eff_ch, params.data.obsid, "${params.data.path}/${params.mpi_dd.solsdir}", params.gains.fmin, params.gains.fmax, params.gains.dd_clusters_to_plot, params.gains.dd_cluster_names, params.gains.stations_to_plot)
            }

        else if (gains_type=="bandpass") {
                PlotSagecalDISolutions(conv_ch.npy, conv_ch.npz, eff_ch, params.data.obsid, params.bandpass.solsdir, params.gains.fmin, params.gains.fmax)
            }
}


workflow POWER_SPECTRUM {
    take:
        sagecal_complete

    main:
        // init_ch = InitPSDB(sagecal_complete, params.data.path)
        ps_dir = "${params.data.path}/${params.pspipe.dir}"
        rev_ch = AddRevision(sagecal_complete, params.data.obsid, params.mpi_dd.output_column, params.data.path, ps_dir, nodes_list[-1], params.pspipe.max_concurrent, params.pspipe.revision, params.pspipe.merge_ms, params.pspipe.aoflag_after_merge_ms) // TODO: stop using only the final node 
        ps_ch = RunPSPIPE(ps_dir, rev_ch.toml_file, params.data.obsid, "${params.data.path}/all_ms_files_003.txt", params.pspipe.merge_ms, params.pspipe.delay_flagger, params.pspipe.vis_flagger, params.pspipe.gpr, params.pspipe.ml_gpr)
        // PlotPowerSpectrum(ps_ch.ready, params.pspipe.dir, rev_ch.toml_file, params.data.obsid)
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Processes
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// process CheckingParams {

//     output:
//     val true, emit: all_params_valid

//     exec:
//     checkParams()
// }


process SavingParams {
    debug true
    publishDir params.logs_dir, mode: 'move'

    input:
    val all_params_valid

    output:
        file 'params.json'
        val true, emit: params_json_written

    script:
    def edited_params =params
    edited_params.remove("init")
    """
        echo '${JsonOutput.prettyPrint(JsonOutput.toJson(edited_params))}' > params.json
    """
}


/*
A task to write out the input MPIRUN text file specifying the nodes to run on 
Also specifies the GPUs pto be used per node (`slots`)

Parameters
----------
ready : string
    A `go-ahead` signal/connector to a previous task
nodes : list
    A list of nodes to be used
slots : int
    Number of GPUs per to be used per node
mpirun_hosts_txt_file : string
    The path of the output txt file

Returns
-------
mpirun_hosts_txt_file : path
    The path to the output file

*/
process MpiRunNodesList {
    // publishDir params.logs_dir

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


/*
A task to write out the input PSSH (Parallel SSH tool) text file specifying the nodes to run on.
See `https://manpages.ubuntu.com/manpages/trusty/man1/parallel-ssh.1.html`

Calls the `write_nodes_for_pssh` function

Parameters
----------
ready : string
    A `go-ahead` signal/connector to a previous task
nodes : list
    A list of nodes to be used
pssh_hosts_txt_file : string
    The path of the output txt file

Returns
-------
pssh_hosts_txt_file : path
    The path to the output file

*/
process PsshNodesList {
    debug true
    // publishDir params.logs_dir

    input:
    val ready
    val nodes
    val pssh_hosts_txt_file

    output:
    val "${pssh_hosts_txt_file}", emit:  pssh_hosts_txt
    val true, emit: all_nodes_standby

    exec:

    write_nodes_for_pssh(nodes, pssh_hosts_txt_file)
}


/*
A task to distribute shapelet files required for the current NCP sky model

Parameters
----------
ready : string
    A `go-ahead` signal/connector to a previous task
pssh_hosts_txt_file : path
    pssh input nodes list file
shapelets_modes : string
    A path to the shapelets files
    default : "${projectDir}/models/shapelets_modes/*modes"
datapath : string
    A path to the directory containing the measurement sets
    It is assumed to be the same for all the nodes used
sim: boolean
    if `true`, then we do not need the shapelets
Returns
-------
true : bool
    The task was successful

*/
process GetModels {
    debug true

    input:
    val ready
    path pssh_hosts_txt_file
    val shapelets_modes
    val datapath
    val sim

    output:
    val true

    script:
    if (sim)

    """
    echo "[nextleap]: We assume the simulation does not need the NCP shapelets model. Shout if otherwise!"
    """

    else

    """
    pssh -v -i -h ${pssh_hosts_txt_file} -t 0 cp ${shapelets_modes} ${datapath} > ${params.logs_dir}/models_copying.log 2>&1
    """
}


/*
A task to apply pre-DI flagging on the data 
Spawns a standalone nextflow job on each node using PSSH

Parameters
----------
ready : string
    A `go-ahead` signal/connector to a previous task
pssh_hosts_txt_file : path
    pssh input nodes list file
preprocessing_file : path
    A standalone nextflow file that does the actual flagging
    default : "${projectDir}/modules/di_high_reg_preprocessing.nf"
datapath : string
    A path to the directory containing the measurement sets
    It is assumed to be the same for all the nodes used
ms_files: path
    A path to a file listing all the measurement sets in a node
    Should exist for each node separately
sim: boolean
    if `true`, then we do not need the shapelets

Returns
-------
true : bool
    The task was successful

*/
process ApplyPreDIFlag {
    debug true

    input:
    val ready
    val pssh_hosts_txt_file
    val preprocessing_file
    val datapath
    val ms_files
    val sim

    output:
    val true

    script:
    if (sim)

    """
    echo "[nextleap]: We assume the simulation does not need any visibilities flagging before DI calibration."
    """

    else

    """
    pssh -v -i -h ${pssh_hosts_txt_file} -t 0 -x "cd ${datapath}; bash" ${params.nextflow_executable} run ${preprocessing_file} --ms_files ${ms_files} > ${params.logs_dir}/preprocessing_di.log 2>&1
    """
}


/*
A task to add a column to a measurement set 
calls `${projectDir}/modules/add_ms_col.nf` using PSSH

Parameters
----------
ready : string
    A `go-ahead` signal/connector to a previous task
pssh_hosts_txt_file : path
    pssh input nodes list file
datapath : string
    A path to the directory containing the measurement sets
    It is assumed to be the same for all the nodes used
ms_files: path
    A path to a file listing all the measurement sets in a node
    Should exist for each node separately
column: string
    Name of the column to be added
    
Returns
-------
true : bool
    The task was successful

*/
process AddColumnToMeasurementSet {
    debug true

    input:
    val ready
    val pssh_hosts_txt_file
    val datapath
    val ms_files
    val column
  

    output:
    val true

    script:

    """
    pssh -v -i -h ${pssh_hosts_txt_file} -t 0 -x "cd ${datapath}; bash" ${params.nextflow_executable} run ${projectDir}/modules/add_ms_col.nf --ms_files ${ms_files} --col ${column} > ${params.logs_dir}/add_${column}_column.log 2>&1
    """
}


/*
A task to apply pre-DD flagging on the data 
Spawns a standalone nextflow job on each node using PSSH

Parameters
----------
ready : string
    A `go-ahead` signal/connector to a previous task
pssh_hosts_txt_file : path
    pssh input nodes list file
preprocessing_file : path
    A standalone nextflow file that does the actual flagging
    default : "${projectDir}/modules/dd_high_reg_preprocessing.nf"
datapath : string
    A path to the directory containing the measurement sets
    It is assumed to be the same for all the nodes used
ms_files: path
    A path to a file listing all the measurement sets in a node
    Should exist for each node separately
sim: boolean
    if `true`, then we do not need the shapelets

Returns
-------
true : bool
    The task was successful

*/
process ApplyPreDDFlag {
    debug true

    input:
    val ready
    val pssh_hosts_txt_file
    val preprocessing_file
    val datapath
    val ms_files
    val sim

    output:
    val true

    script:
    if (sim)

    """
    echo "[nextleap]: We assume the simulation does not need any visibilities flagging before DD calibration."
    """

    else

    """
    pssh -v -i -h ${pssh_hosts_txt_file} -t 0 -x "cd ${datapath}; bash" ${params.nextflow_executable} run ${preprocessing_file} --ms_files ${ms_files} > ${params.logs_dir}/preprocessing_dd.log 2>&1
    """
}


/*
A task to flag and average data before DI calibration
Spawns a standalone nextflow job on each node using PSSH

Parameters
----------
ready : string
    A `go-ahead` signal/connector to a previous task
pssh_hosts_txt_file : path
    pssh input nodes list file
dp3_average_to_002_file : path
    A standalone nextflow file that does the actual averaging and flagging
    default : "${projectDir}/modules/gen_002_vis.nf"
datapath : string
    A path to the directory containing the measurement sets
    It is assumed to be the same for all the nodes used
ms_files: path
    A path to a file listing all the measurement sets in a node
    Should exist for each node separately

Returns
-------
true : bool
    The task was successful

*/
process FlagAndAverageVisTo002 {
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
    pssh -v -i -h ${pssh_hosts_txt_file} -t 0 -x "cd ${datapath}; bash" ${params.nextflow_executable} run ${dp3_average_to_002_file} --path ${datapath} > ${params.logs_dir}/flag_average_to_002.log 2>&1
    """
}


/*
A task to average data before DD calibration
Spawns a standalone nextflow job on each node using PSSH

Parameters
----------
ready : string
    A `go-ahead` signal/connector to a previous task
pssh_hosts_txt_file : path
    pssh input nodes list file
dp3_average_to_002_file : path
    A standalone nextflow file that does the actual averaging and flagging
    default : "${projectDir}/modules/gen_003_vis.nf"
datapath : string
    A path to the directory containing the measurement sets
    It is assumed to be the same for all the nodes used
ms_files: path
    A path to a file listing all the measurement sets in a node
    Should exist for each node separately

Returns
-------
true : bool
    The task was successful

*/
process AverageVisTo003 {
    debug true

    input:
    val ready
    val pssh_hosts_txt_file
    val dp3_average_to_003_file
    val datapath
    val data_column
    val ms_files

    output:
    val true

    script:
    //--ms_files ${ms_files}
    """
    pssh -v -i -h ${pssh_hosts_txt_file} -t 0 -x "cd ${datapath}; bash" ${params.nextflow_executable} run ${dp3_average_to_003_file} --path ${datapath} --data_column ${data_column} > ${params.logs_dir}/average_to_003.log 2>&1
    """
}


/*
A task to run sagecal in standalone mode (No regularisation)
Spawns a standalone nextflow job on each node using PSSH


Parameters
----------
ready : string
    A `go-ahead` signal/connector to a previous task
command : string
    The sagecal command to run
    Made using `generateSagecalCommand()`
pssh_hosts_txt_file : path
    pssh input nodes list file
datapath : string
    A path to the directory containing the measurement sets
    It is assumed to be the same for all the nodes used
standalone_sage_nf_file : path
    A standalone nextflow file that runs the actual sagecal job
    default : "${projectDir}/standalone.nf"
ms_files: path
    A path to a file listing all the measurement sets in a node
    Should exist for each node separately
solsdir: path
    A path to where the solutions will be stored

Returns
-------
true : bool
    The task was successful

*/
process SagecalStandalone {
    debug true
    // errorStrategy 'ignore'

    input:
    val ready
    val command
    path pssh_hosts_txt_file
    val datapath
    val standalone_sage_nf_file
    val ms_files
    val solsdir
    val time_limit

    output:
    val true, emit: sagecal_complete

    script:
    """
    pssh -v -i -h ${pssh_hosts_txt_file} -t 0 -x "cd ${datapath}; bash" ${params.nextflow_executable} run ${standalone_sage_nf_file} --ms_files ${ms_files} --solsdir ${solsdir} --time_limit ${time_limit} --command "'${command}'"  > ${params.logs_dir}/sagecal_standalone.log 2>&1
    """
}

/*
A task to run sagecal in MPI mode (With regularisation, Either DI or DD)
Spawns a standalone nextflow job on each node using PSSH


Parameters
----------
ready : string
    A `go-ahead` signal/connector to a previous task
sagecal_command : string
    The sagecal command to run
    Made using `generateSagecalCommand()`
datapath : string
    A path to the directory containing the measurement sets
    It is assumed to be the same for all the nodes used
pssh_hosts_txt_file : path
    pssh input nodes list file
standalone_sage_nf_file : path
    A standalone nextflow file that runs the actual sagecal job
    default : "${projectDir}/standalone.nf"
solsdir: path
    A path to where the solutions will be stored
ms_pattern: string
    A glob pattern for the measurement sets
    Made by `lofarDefaultMsnamePattern(stage_number)` for a given processing stage
    used to access sagecal `MS.solutions` files and move them to `solsdir`

Returns
-------
true : bool
    The task was successful

*/
process SagecalMPI {
    debug true
    cpus params.mpi_dd.number_of_threads
    // time '1hour 30minutes'

    input:
    val ready
    val sagecal_command
    val datapath
    val pssh_hosts_txt_file
    val solsdir
    val ms_pattern

    output:
    val true , emit: sagecal_complete
    //path(consensus_solutions_file), path(logfile),

    script:
    //when done we remove the many shapelets .modes files
    """
    cd ${datapath}
    time ${sagecal_command}
    pssh -v -i -h ${pssh_hosts_txt_file} -t 0 "mkdir -p ${solsdir}; mv ${datapath}/${ms_pattern}.solutions ${solsdir}"
    #rm *.modes
    """
    //the rm modes should give the absolute path
}


/*
A task to run wsclean
#TODO: add more options
#TODO: make all_sky_plotting and stats optional

Parameters
----------
ready : string
    A `go-ahead` signal/connector to a previous task
scale : string
    -scale
size : string
    size
datacolumn : 
    -data-column
name : path
    -name
ms_fyl: path
    Measurement set

Returns
-------
true : bool
    The task was successful
*.fits: paths
    fits images
*.png : paths
    all-sky plot

*/
process ImageWithWSClean { 
    publishDir "${params.data.path}/${name}_images", pattern: "*.fits", mode: "move", overwrite: true
    publishDir "${params.data.path}/${name}_images", pattern: "*.png", mode: "move", overwrite: true

    input:
    val ready
    val scale
    val size
    val weight
    val minuv_lambda
    val maxuv_lambda
    val polarisation
    val threads
    val datacolumn
    val name
    val msfyl

    output:
    val true , emit: wsclean_complete
    path "*.fits"
    path "*.png"

    script:
    List mslist = file(msfyl).readLines()
    String mses = mslist.collect {"${it}"}.join(" ")
    """
    wsclean -interval 0 361 -data-column ${datacolumn} -gridder wgridder -wgridder-accuracy 1e-5 -reorder -make-psf -scale ${scale} -size ${size} ${size} -weight ${weight} -minuv-l ${minuv_lambda} -maxuv-l ${maxuv_lambda} -pol ${polarisation} -name ${name} -j ${threads} ${mses}  > ${params.logs_dir}/wsclean.log 2>&1
    python3 ${projectDir}/templates/read_fits_image.py -i ${name}-image.fits
    """
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Closures
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def checkParams() {
    //check obsid and nodes
    assert params.data.obsid : log.error("Error:'Observation ID is required")
    assert params.data.label : log.error("Error:'Observation ID label is required")
    assert params.data.path : log.error("Error:'--data.path' needed")
    assert params.cluster.nodes : log.error("Error: Please provide a comma separated string of nodes. Try --cluster.nodes '100,101,102,103,104'")

    // if (params.tasks.contains("ps")){
    //     dir = file("${params.data.path}/${params.pspipe.dir}", type: 'dir', checkIfExists: false) //I don't know how to make it fail if it already exists.
    //     dir.exists() ? log.info("[nextleap]: ${params.data.path}/${params.pspipe.dir} ps dir already exists. will be overwritten!") : log.info("[nextleap]: pspipe dir ${params.pspipe.dir} ")
    // }
    return true
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
    else {
        log.error("Error: The `cluster.nodes` parameter is not valid. Got `--params.cluster.nodes=${params.cluster.nodes}`")
        exit 0
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
    assert dir.isDirectory() : log.error("Error: datapath param is not a valid path on ${node}. Got '--datapath=${params.data.path}'")

    //check that we have a least one file with the given glob pattern
    ms_files_path = "${files_path}/${files_glob_pattern}"
    list_of_msets = file(ms_files_path, type: 'dir', glob: true, checkIfExists: true)
    assert list_of_msets.size()>0 : log.error("Error: No measurement sets matching pattern '${ms_files_path}' found.")

    return list_of_msets
}

/*
Given the tasks requested using 'params.task' return a map of each task with its respective 'stage_number'
This is for the conventional MS naming
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
    log.info("[nextleap]: initial stage number: ${initial_stage_number}")
    return initial_stage_number 
}


/*
Some proesses may need a list of all measurement sets across all nodes. e.g pspipe
Make such a file given a map containing the subbands per node and the stage number of the data needed
This list is written once in the master node
*/
def allMsetsPerStageNumber(stage_number, sub_bands_per_node){
    mses = []
    for (node in sub_bands_per_node) {
        for (sb in node.value) {
            sb = String.format("%03d", sb)
            fyl = lofarDefaultMsnamePattern(stage_number).replace("???", "${sb}")
            full_fyl_path = "/net/${node.key}${params.data.path}/${fyl}"
            mses.add(full_fyl_path)
        }
    }

    txt = "/net/${params.cluster.masternode}/${params.data.path}/all_ms_files_${stage_number}.txt"

    writeListToTxt(txt, mses)

    return txt
}

/*
 This is similar to the `allMsetsPerStageNumber` but it makes a list in each node containing the ms files in that node.
 It is necessary fo rtasks that are to be run in parallel individually in each node e.g preprocessing each msfile, bandpass calibration, data averagineg etc...
*/
def writeMSListPerNodeStageNumber() {
    requested_tasks = getRequiredtasksStageNumbers()
    stage_num_per_task = requested_tasks.values() as List

    for (n in stage_num_per_task.unique()) {
        _writeMSListPerNodeStageNumber(n, params.data.sub_bands_per_node)
    }
}

//dummy method for the above one
def _writeMSListPerNodeStageNumber(stage_number, sub_bands_per_node){
    stage_number = String.format("%03d", stage_number)
    for ( node in sub_bands_per_node ) {
        mses = []
        txt = "/net/${node.key}/${params.data.path}/ms_files_${stage_number}.txt"
        for (sb in node.value) {
            sb = String.format("%03d", sb)
            prefix = lofarDefaultMsnamePattern(stage_number).replace("???", "${sb}")
            ms = "/net/${node.key}/${params.data.path}/${prefix}"
            mses.add(ms)
        }
        writeListToTxt(txt, mses)
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
    log.info("[nextleap]: initial ms pattern: ${initial_ms_pattern}")

    params.data.sub_bands_per_node = [:]

    def isSubBand = { it.toString().split("SB")[1].split("_uv_${initial_stage_number}_${params.data.label}.MS")[0].toInteger()}

    // def isSubBand = { it.toString().split("_")[-4].split("SB")[1].toInteger()}

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


/*
Generate a sagecal command for the different LOFAR EoR calibration stages
*/
def generateSagecalCommand() {
    if (params.tasks.contains("mpi_di")){
        params.mpi_di.ms_pattern = lofarDefaultMsnamePattern(params.mpi_di.stage_number)
        params.mpi_di.sagecal_command =  SagecalMPIDICommand()
    }
    if (params.tasks.contains("bandpass")){
        params.bandpass.ms_pattern = lofarDefaultMsnamePattern(params.bandpass.stage_number)
        params.bandpass.sagecal_command =  SagecalBandpassCommand()
    }
    if (params.tasks.contains("mpi_dd")){
        params.mpi_dd.ms_pattern = lofarDefaultMsnamePattern(params.mpi_dd.stage_number)
        params.mpi_dd.sagecal_command =  SagecalMPIDDCommand()
    }
}


def SagecalBandpassCommand() {
    // if (params.bandpass.sagecal_command) {
    //     log.info("[nextleap]: using user provided sagecal bandpass command")
    //     return params.bandpass.sagecal_command
    // }
    // else {
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
-k ${params.bandpass.gains_correction_cluster_id} \
-a ${params.bandpass.action} \
-D ${params.bandpass.enable_diagnostics} \
-E ${params.bandpass.use_gpu}
""".stripIndent()
    
        bandpass_command = "sagecal_gpu " + sage_bandpass_command.strip()
        bandpass_command= bandpass_command.stripIndent()

        return bandpass_command
    // }
}

def SagecalMPIDICommand(){
    // if (params.mpi_di.sagecal_command) {
    //     log.info("[nextleap]: using user provided sagecal MPI DI command")
    //     return params.mpi_di.sagecal_command
    // }
    // else{
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
-k ${params.mpi_di.gains_correction_cluster_id} \
-A ${params.mpi_di.admm_iterations} \
-P ${params.mpi_di.consensus_polynomial_terms} \
-G ${params.mpi_di.admm_rho_file} \
-p ${params.mpi_di.solutions_file} \
-f \"${params.mpi_di.ms_pattern}\" \
-E ${params.mpi_di.use_gpu} \
-V > ${params.mpi_di.logfile}  2>&1 
""".stripIndent()

        if (params.mpi_di.constant_rho_value){
        sage_mpi_di_command = sage_mpi_di_command.replace("-G ${params.mpi_di.admm_rho_file}", "-r ${params.mpi_di.constant_rho_value}")
    }

        mpi_command = make_mpi_command()

        mpi_di_command = mpi_command + " " +  sage_mpi_di_command.strip()
        mpi_di_command= mpi_di_command.stripIndent()

        return mpi_di_command
    // }
}


def SagecalMPIDDCommand(){
    // if (params.mpi_dd.sagecal_command) {
    //     log.info("[nextleap]: using user provided sagecal MPI DD command")
    //     return params.mpi_dd.sagecal_command
    // }
    // else {
    logfile="${params.logs_dir}/sagecal_mpi_dd.log"

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
-V > ${logfile}  2>&1 
""".stripIndent()

    if (params.mpi_dd.constant_rho_value){
        sage_mpi_dd_command = sage_mpi_dd_command.replace("-G ${params.mpi_dd.admm_rho_file}", "-r ${params.mpi_dd.constant_rho_value}")
    }

    if (params.mpi_dd.ntimeslots_to_calibrate){
        sage_mpi_dd_command = sage_mpi_dd_command.replace("-V > ${logfile}", "-T ${params.mpi_dd.ntimeslots_to_calibrate} -V > ${logfile}")
    }

    if (params.mpi_dd.enable_spatial_reg){
        if (params.mpi_dd.constant_alpha_value){
            sage_mpi_dd_command = sage_mpi_dd_command.replace("-V > ${logfile}", "-X ${params.mpi_dd.spatial_reg.lambda},${params.mpi_dd.spatial_reg.mu},${params.mpi_dd.spatial_reg.n0},${params.mpi_dd.spatial_reg.fista_maxiter},${params.mpi_dd.spatial_reg.candence} -u ${params.mpi_dd.spatial_reg.constant_alpha_value} -V > ${logfile}")
        }

        else {
            sage_mpi_dd_command = sage_mpi_dd_command.replace("-V > ${logfile}", "-X ${params.mpi_dd.spatial_reg.lambda},${params.mpi_dd.spatial_reg.mu},${params.mpi_dd.spatial_reg.n0},${params.mpi_dd.spatial_reg.fista_maxiter},${params.mpi_dd.spatial_reg.candence} -V > ${logfile}")
        }

    }

        mpi_command = make_mpi_command()

        mpi_dd_command = mpi_command.strip() + " " +  sage_mpi_dd_command.strip()
        mpi_dd_command = mpi_dd_command.stripIndent()

        return mpi_dd_command
    // }
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


def helpMessage() {
log.info """
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|

DESCRIPTION! 
    The workflow configuration in the current version runs the end to end LOFAR production pipeline.
    It includes ALL analysis tasks (plus some diagnostics steps) in the following sequence:
        1: gen002vis
        2: mpi_di
        3: bandpass
        4: wsclean*
        5: gen003vis
        6: mpi_dd
        7: gains
        8: wsclean*
        9: ps
    Different tasks combination can also be used.
    However, the possibilities supported might not be exhaustive.
    Let us know if your disired pcombination is not covered.
    For sugestions/questions talk to the author of NextLEAP
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|


TYPICAL USAGE COMMAND:
    This pipeline runs in 2 commands:

    |-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
    |nextflow run main.nf -profile standard,gen002vis,mpi_di,mpi_dd,bandpass,gains,gen003vis,wsclean,ps --data.obsid L254871 --data.path /path/to/MS/files/ --data.label R --cluster.nodes 126,127,128,129 -entry init      |
    |-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
    |nextflow run main.nf -profile standard,gen002vis,mpi_di,mpi_dd,bandpass,gains,gen003vis,wsclean,ps -params-file logs/params.json                                                                                       |
    |-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|

    Command 1:  Initiates the minimum required inputs.

            A json configuration file named `params.json` is output with more variables.
            These variables are set to the default settings in the LOFAR EoR production pipeline.
            Take some time to confirm/change the config settings in the file.

    Command 2:  Runs the pipeline
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|

MANDATORY ARGUMENTS:
    --obsid                         Observation ID. e.g L254874
    --data.path                     The path where data is stored. Should be valid for all the nodeswhere data is distributed.
    --cluster.nodes                 Node numbers where data is distributed as a comma separated string. e.g. "125,126,127,128,129" or a single number if using a single node. The last node in the list is the masternode.
    --cluster.slots                         Number of slots per node for mpi run. For now, all given nodes are assigned the same slots number so maybe don't use nodes 119 and 124 on DAWN.
    -profile                        The configuration profile lists the tasks to be run and where to run them.                                               
    --label                         An extra data string label, like "R" or '2022' previously used for LOFAR NCP data to refer to a specific processing season.


OPTIONAL ARGUMENTS:
    --help                          Display this usage statement, and exit.
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|

"""
}