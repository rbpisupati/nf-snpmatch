#!/usr/bin/env nextflow
/*
Simply run nextflow pipeline as:
// Better to use imputed dataset for snpmatch cross

nextflow run simulate.nf --f1 num_snps_for_f1 --accList "accList.txt" --err_rate 0.01 --db hdf5_file --db_acc hdf5_acc_file --outdir output_folder
*/
/*
 * SET UP CONFIGURATION VARIABLES
*/
params.accList = false
params.f1 = false  // provide number of SNPs to be used here for F1
params.het_fraction = 0.1 // only for simulating F1 

params.outdir = 'simulate'
params.err_rate = 0.01
params.skip_db_hets = false

// databases
params.db = "/groups/nordborg/projects/the1001genomes/scratch/rahul/101.VCF_1001G_1135/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.hdf5"

db_file = Channel
    .fromPath ( params.db )
    .map{ [ file("${it}"), file("${it.parent}/${it.baseName}.acc.hdf5") ] }
    .ifEmpty { exit 1, "please provide hdf5 files as a database" }

/*
  Do all pairwise combinations for an F1
*/

if (params.f1 ){
    //input acc list
    ecotype_ids = Channel
        .fromPath ( params.accList )
        .splitCsv(header: false)
        .ifEmpty { exit 1, "Cannot find any input files matching: ${params.input}\nNB: Path needs to be enclosed in quotes!\n" }
        .map{ vcf -> vcf }
        .collect()

    input_accs = ecotype_ids
        .map { vcf -> [vcf, vcf].combinations().findAll { x, y -> x < y } }

    input_files_dbs = db_file.spread(input_accs) //.combine(db_file).println()

} else {

    num_positions = Channel
        .from(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,2000000)

    input_accs = Channel
        .fromPath ( params.accList )
        .splitCsv()
        .ifEmpty { exit 1, "Cannot find any input files matching" }

    input_files_dbs = db_file.combine(input_accs).spread(num_positions)

}


process simulateSNPs {
    tag { "${acc_id}_x${num_snps_or_f1_father}" }
    publishDir "$params.outdir", mode: 'copy', saveAs: {filename ->
            if ( filename.indexOf("bed") > 0 ) "snps_${num_snps_or_f1_father}/$filename"
            else "snpmatch/$filename"  }
    errorStrategy { task.exitStatus in [143,137] ? 'retry' : 'ignore' }

    input:
    set file(f_db), file(f_db_acc), val(acc_id), val(num_snps_or_f1_father) from input_files_dbs

    output:
    file "${acc_id}*.bed" into simulated_bed
    file "*match_${acc_id}*" into snpmatch_output

    script:
    skip_db_hets = params.skip_db_hets != false ? "--skip_db_hets" : ''
    if (params.f1 ){
        """
        snpmatch simulate --f1 --het_frac $params.het_fraction -v -d $f_db -e $f_db_acc -a "${acc_id}x${num_snps_or_f1_father}" -n $params.f1 -o ${acc_id}x${num_snps_or_f1_father}.bed -p $params.err_rate
        mkdir -p csmatch_${acc_id}x${num_snps_or_f1_father}
        snpmatch cross -v $skip_db_hets -d $f_db -e $f_db_acc -i ${acc_id}x${num_snps_or_f1_father}.bed -o csmatch_${acc_id}x${num_snps_or_f1_father}/${acc_id}x${num_snps_or_f1_father}.snpmatch
        """
    } else {
        """
        snpmatch simulate -v -d $f_db -e $f_db_acc -a $acc_id -n $num_snps_or_f1_father -o ${acc_id}_${num_snps_or_f1_father}.bed -p $params.err_rate
        mkdir -p snpmatch_${acc_id}_${num_snps_or_f1_father}
        snpmatch inbred --refine -v $skip_db_hets -d $f_db -e $f_db_acc -i ${acc_id}_${num_snps_or_f1_father}.bed -o  snpmatch_${acc_id}_${num_snps_or_f1_father}/${acc_id}_${num_snps_or_f1_father}.snpmatch
        """
    }
}

input_csv = snpmatch_output.collect()

process make_csv_simulate {
    publishDir "$params.outdir", mode: 'copy'

    input:
    file "*" from input_csv

    output:
    file "intermediate_modified.csv" into output_csv

    script:
    if (params.f1 ){
        """
        python  $workflow.projectDir/bin/02_makeCSVTable_csmatch.py --dirs -i ./ -o intermediate_modified.csv -f $params.outdir
        """
    } else {
        """
        python $workflow.projectDir/bin/01_makeCSVTable_inbred.py --dirs -i ./ -o intermediate_modified.csv -f $params.err_rate
        """
    }
}
