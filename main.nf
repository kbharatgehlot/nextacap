#!/usr/bin/env nextflow

// Show help message and exit
if (params.help) {
    helpMessage()
    exit 0
}

//import workflows
include {
    VET_PARAMS;
    SET_NODES;
    Production;
    ProductionPipelineFromFirstDI;
    ProductionPipelineFromBandpass;
    ProductionPipelineFromPostDIAveraging;
    ProductionPipelineFromDD;
} from './nextleap_workflows.nf'

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
workflow init {
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
            Production()
        }
        else if (requested_tasks_set == pipeline_from_first_di_set) {
            defaultLOFArEoRPipelineFrom002Message()
            params.tasks = PRODUCTION_PIPELINE_FROM_FIRST_DI.split(",")
            ProductionPipelineFromFirstDI()
        }
        else if (requested_tasks_set == pipeline_from_bandpass_set) {
            defaultLOFArEoRPipelineFromBandpassMessage()
            params.tasks = PRODUCTION_PIPELINE_FROM_BANDPASS.split(",")
            ProductionPipelineFromBandpass()
        }
        else if (requested_tasks_set == pipeline_from_post_di_averaging_set) {
            defaultLOFArEoRPipelineFromPostDIAveraging()
            params.tasks = PRODUCTION_PIPELINE_FROM_POST_DI_AVERAGING.split(",")
            ProductionPipelineFromPostDIAveraging()
        }
        else if (requested_tasks_set == pipeline_from_dd_set) {
            defaultLOFArEoRPipelineFromDDMessage()
            params.tasks = PRODUCTION_PIPELINE_FROM_DD.split(",")
            ProductionPipelineFromDD()
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
Info: Running default LOFAR EoR KSP production pipeline:
Info: ${PRODUCTION_PIPELINE}
"""
}

def defaultLOFArEoRPipelineFrom002Message() {
    log.info """
Info: Running default LOFAR EoR KSP production pipeline from the the first Direction Independent calibration step (fast time-varying effects with spectral regularisation):
Info: ${PRODUCTION_PIPELINE_FROM_FIRST_DI}
"""
}

def defaultLOFArEoRPipelineFromBandpassMessage() {
    log.info """
Info: Running default LOFAR EoR KSP production pipeline from the second Direction Independent calibration step (bandpass without spectral regularisation):
Info: ${PRODUCTION_PIPELINE_FROM_FIRST_DI}
"""
}

def defaultLOFArEoRPipelineFromPostDIAveraging() {
    log.info """
Info: Running default LOFAR EoR KSP production pipeline from an averaging step towards Direction Dependent Calibration:
Info: ${PRODUCTION_PIPELINE_FROM_POST_DI_AVERAGING}
"""
}

def defaultLOFArEoRPipelineFromDDMessage() {
    log.info """
Info: Running default LOFAR EoR KSP production pipeline from the Direction Dependent calibration step:
Info: ${PRODUCTION_PIPELINE_FROM_DD}
"""
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
