#!/usr/bin/env nextflow

process makeEffNr {
    publishDir params.data.path
    errorStrategy 'ignore'

    input:
    val ready
    path clusterfile

    output:
    path "eff_nr_${clusterfile.Name}.npy"

    script:
    """
    python3 ${projectDir}/templates/create_eff_nr.py --eff_nr eff_nr_${clusterfile.Name} -c ${clusterfile}
    """
}

process convertSageSolutions {
    debug true
    publishDir params.data.path
    errorStrategy 'ignore'

    input:
    path eff_nr_file
    val obsid
    val ms_pattern
    val nodes
    val pssh_hosts_txt_file
    val solsdir
    val datapath
    val stage_number

    output:
        val("${obsid}.npz"), emit: npz
        val("${obsid}.npy"), emit: npy
        val(true), emit: ready
        
    script:
    allnodes = nodes.split(",").collect{"node${it}"}.join(" ")
    logs = "${params.logs_dir}/convert_sols.log"
    """
    #pssh -v -i -h ${pssh_hosts_txt_file} -t 0 "mkdir -p ${solsdir}; mv ${datapath}/${ms_pattern}.solutions ${solsdir}"

    python3 ${projectDir}/templates/convert_sage_solutions.py -o ${obsid} -m ${datapath}/*.MS -p ${solsdir} -d \${PWD} -n ${allnodes} -c 0 2000 --eff_nr ${eff_nr_file} --pid ${stage_number} > ${logs}
    mv *.npz *.npy ${solsdir}
    """
}


process convertSageZSol {
    debug true
    publishDir params.data.path
    errorStrategy 'ignore'

    input:
    val ready
    path eff_nr_file
    val obsid
    val solsdir
    val datapath

    output:
        path "${obsid}_zsol", optional: true

    script:
    logs = "${params.logs_dir}/convert_Zsols_${obsid}.log"
    """
    mv ${datapath}/*Zsol ${solsdir}

    python3 ${projectDir}/templates/convert_sage_zsol.py --output_name ${obsid}_zsol --eff_nr ${eff_nr_file} ${solsdir}/sagecal_Zsol ${solsdir}/${obsid}.npz > ${logs} 2>&1

    """
}


process plotDDSageSols {
    debug true
    publishDir params.data.path
    errorStrategy 'ignore'

    input:
    val sols_npy
    val sols_npz
    path eff_nr_file
    val obsid
    val solsdir
    val fmin
    val fmax

    output:
        // path "${obsid}_zsol"
        val true

    script:
    logs = "${params.logs_dir}/plot_dd_sols_${obsid}.log"
    """
    python3 ${projectDir}/templates/plot_dd_cal_solutions.py --fmin ${fmin} --fmax ${fmax} --out_dir ${solsdir} ${solsdir}/${sols_npy} --eff_nr ${eff_nr_file} > ${logs} 2>&1
    """
}

process plotDISageSols {
    debug true
    publishDir params.data.path
    errorStrategy 'ignore'

    input:
    val sols_npy
    val sols_npz
    path eff_nr_file
    val obsid
    val solsdir
    val fmin
    val fmax

    output:
        // path "${obsid}_zsol"
        val true

    script:
    logs = "${params.logs_dir}/plot_di_sols_${obsid}.log"
    """
    python3 ${projectDir}/templates/plot_di_cal_solutions.py --fmin ${fmin} --fmax ${fmax} --out_dir ${solsdir} --cluster 0 ${solsdir}/${sols_npy} > ${logs} 2>&1
    """
}

// Not sure what this does
// ./convert_sage_solutions_back.py -o ${MS}  -p ${datapath}solutions_DI${label}_high_reg_${step}/ -d  /net/${masternode}${datapath}solutions_DI${label}_high_reg_${step}/${MS}_zsol.npy -l ${label}  --pid 2 -n 101..115 -m  /net/${masternode}${datapath}solutions_DI${label}_high_reg_${step}/${MS}.npz
