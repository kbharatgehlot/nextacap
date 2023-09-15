#!/usr/bin/env nextflow

/*
useful wsclean related functionalities
see: https://wsclean.readthedocs.io/en/latest/index.html
*/

params.scale=null
params.size=null
params.weight=null
params.minuv_lambda=null
params.maxuv_lambda=null
params.polarisation=null
params.threads=null
params.output_column=null
params.name=null
params.mslist=null
params.extra_options=null

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
    publishDir "${params.data.path}/${fits_filename}_images", pattern: "*.fits", mode: "move", overwrite: true

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
    val fits_filename
    val msfyl
    val extra_options

    output:
    val true , emit: wsclean_complete
    path "*.fits"
    

    script:
    List mslist = file(msfyl).readLines()
    String measurement_sets = mslist.collect {"${it}"}.join(" ")
    logs = "${params.logs_dir}/wsclean_image.log"

    if (extra_options)

    """
    wsclean -data-column ${datacolumn} -gridder wgridder -wgridder-accuracy 1e-5 -reorder -make-psf -scale ${scale} -size ${size} ${size} -weight ${weight} -minuv-l ${minuv_lambda} -maxuv-l ${maxuv_lambda} -pol ${polarisation} -name ${fits_filename} -j ${threads} ${extra_options} ${measurement_sets}  > ${logs} 2>&1
    """

    else:

    """
    wsclean -data-column ${datacolumn} -gridder wgridder -wgridder-accuracy 1e-5 -reorder -make-psf -scale ${scale} -size ${size} ${size} -weight ${weight} -minuv-l ${minuv_lambda} -maxuv-l ${maxuv_lambda} -pol ${polarisation} -name ${fits_filename} -j ${threads} ${measurement_sets}  > ${logs} 2>&1
    """
}

//https://wsclean.readthedocs.io/en/latest/prediction.html#prediction
process PredictWithWSClean {
    
    input:
    val fits_filename
    path ms
    val nchans_out

    output:
    val true , emit: wsclean_complete
    path "*.fits"

    script:
    logs = "${params.logs_dir}/wsclean_predict.log"
    """
    wsclean -predict -channels-out ${nchans_out} -name ${fits_filename} ${ms} > {logs} 2>&1
    """
}


process PlotFitsImage {
    publishDir "${params.data.path}/${name}_images", pattern: "*.fits", mode: "move", overwrite: true

    input:
    val fits_filename

    output:
    path "*.png"

    script:
    logs = "${params.logs_dir}/plot_fits_image.log"
    """
    python3 ${projectDir}/templates/read_fits_image.py -i ${fits_filename} > ${logs} 2>&1
    """
}