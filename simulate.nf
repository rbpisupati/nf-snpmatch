#!/usr/bin/env nextflow
/*
Simply run nextflow pipeline as:
// Better to use imputed dataset for snpmatch cross

nextflow run simulate.nf --input_acclist "*csv" --err_rate 0.01 --db hdf5_file --db_acc hdf5_acc_file --outdir output_folder
*/
/*
 * SET UP CONFIGURATION VARIABLES
*/
params.input_acclist = false

params.outdir = 'simulate'
params.err_rate = 0.01

// databases
params.db = "/groups/nordborg/projects/the1001genomes/scratch/rahul/101.VCF_1001G_1135/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.hdf5"
params.db_acc= "/groups/nordborg/projects/the1001genomes/scratch/rahul/101.VCF_1001G_1135/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.acc.hdf5"

num_positions = Channel
    .from(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,2000000)

//input acc list
input_accs = Channel
    .fromPath ( params.input_acclist )
    .splitCsv()
    .ifEmpty { exit 1, "Cannot find any input files matching" }

db_file = Channel
    .fromPath ( params.db )
    .ifEmpty { exit 1, "please provide hdf5 file as a database" }

db_acc_file = Channel
    .fromPath ( params.db_acc )
    .ifEmpty { exit 1, "please provide hdf5 acc file as a database" }


input_files_dbs = input_accs.combine(db_file).combine(db_acc_file).spread(num_positions)

process identify_libraries {
    tag { "${acc_id}_${num_snps}" }
    publishDir "$params.outdir/snps_${num_snps}", mode: 'copy', saveAs: {filename -> 
            if ( filename.indexOf("bed") > 0 ) "$filename"
            else null  }
    errorStrategy { task.exitStatus in [143,137] ? 'retry' : 'ignore' }

    input:
    set val(acc_id), file(f_db), file(f_db_acc), val(num_snps) from input_files_dbs

    output:
    file "${acc_id}_${num_snps}.bed" into simulated_bed
    file "${acc_id}_${num_snps}.snpmatch*" into snpmatch_output

    script:
    """
    snpmatch simulate -v -d $f_db -e $f_db_acc -a $acc_id -n $num_snps -o ${acc_id}_${num_snps}.bed -p $params.err_rate
    snpmatch inbred --refine  -v -d $f_db -e $f_db_acc -i ${acc_id}_${num_snps}.bed -o ${acc_id}_${num_snps}.snpmatch
    """
}

input_csv = snpmatch_output.collect()

process make_csv_simulate {
    publishDir "$params.outdir", mode: 'copy'

    input:
    file "*" from input_csv

    output:
    file "intermediate_modified.csv" into output_csv

    """
    python $workflow.projectDir/scripts/01_makeCSVTable_inbred.py -o intermediate_modified.csv -f $params.err_rate
    """
}
