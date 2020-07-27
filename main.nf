#!/usr/bin/env nextflow
/*
Simply run nextflow pipeline as:
// Better to use imputed dataset for snpmatch cross

nextflow run submit.snpmatch.nf --input "*vcf" --db hdf5_file --outdir output_folder
*/
/*
 * SET UP CONFIGURATION VARIABLES
*/
params.input = false
// change input_column when providing bed file as input to read nth column
params.input_column = false
params.project = "the1001genomes"
params.outdir = 'snpmatch_1135g'
params.func = 'inbred'
params.genome = "athaliana_tair10"
// databases
params.db = "/groups/nordborg/projects/the1001genomes/scratch/rahul/101.VCF_1001G_1135/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.hdf5"


//input files
input_files = Channel
    .fromPath ( params.input )
    .map { [ "$it.baseName", "${it.getParent()}", file("$it") ] }
    .ifEmpty { exit 1, "Cannot find any input files matching: ${params.input}\nNB: Path needs to be enclosed in quotes!\n" }

db_file = Channel
    .fromPath ( params.db )
    .map{ [ file("${it}"), file("${it.parent}/${it.baseName}.acc.hdf5") ] }
    .ifEmpty { exit 1, "please provide hdf5 files as a database" }


process parse_inputfiles {
  tag { "${prefix}" }
  // publishDir "${input_folder}", mode: 'copy', overwrite: false
  storeDir "${input_folder}"

  input:
  set val(prefix), val(input_folder), file(input_file) from input_files

  output:
  set val(prefix), file("${input_file}.snpmatch.npz") into input_gzips
  file "${input_file}.snpmatch.stats.json" into input_stats

  script:
  if (params.input_column){
    """
    awk '{print \$1 "\t" \$2 "\t" \$${params.input_column}}' $input_file > ${input_file}.col.bed
    snpmatch parser -v -i ${input_file}.col.bed -o ${input_file}.snpmatch
    """
  } else{
    """
    snpmatch parser -v -i $input_file -o ${input_file}.snpmatch
    """
  }
}

input_files_dbs = input_gzips.combine(db_file)

if (params.func == 'inbred'){
  process identify_libraries {
    tag { "${prefix}" }
    publishDir "$params.outdir", mode: 'copy'
    errorStrategy { task.exitStatus in [143,137] ? 'retry' : 'ignore' }

    input:
    set val(prefix), file(input_npz), file(f_db), file(f_db_acc) from input_files_dbs

    output:
    file "snpmatch_${prefix}" into snpmatch_output

    script:
    """
    mkdir -p snpmatch_${prefix}
    snpmatch inbred -v -d $f_db -e $f_db_acc -i $input_npz -o snpmatch_${prefix}/${prefix} --refine
    """
  }

  input_csv = snpmatch_output.collect()

  process make_csv_inbred {
    publishDir "$params.outdir", mode: 'copy'

    input:
    file "*" from input_csv

    output:
    file "intermediate_modified.csv" into output_csv

    """
    mkdir all_results
    ln -s -r snpmatch_*/* all_results
    python $workflow.projectDir/scripts/01_makeCSVTable_inbred.py -i all_results -o intermediate_modified.csv -f $params.outdir
    """
  }

}

if (params.func == 'cross'){
  process cross_libraries {
    tag { "${prefix}" }
    publishDir "$params.outdir", mode: 'copy'
    errorStrategy { task.exitStatus in [143,137] ? 'retry' : 'ignore' }
    // cross generally puts out many errors based on the number of SNPs in a window and chromosome

    input:
    set val(prefix), file(input_npz), file(f_db), file(f_db_acc) from input_files_dbs

    output:
    file "csmatch_${prefix}" into snpmatch_output

    script:
    """
    mkdir -p csmatch_${prefix}
    snpmatch cross -v -d $f_db -e $f_db_acc -i $input_npz -o csmatch_${prefix}/${prefix}.csmatch --genome $params.genome
    """
  }

  input_csv = snpmatch_output.collect()

  process make_csv_cross {
    publishDir "$params.outdir", mode: 'copy'

    input:
    file "*" from input_csv

    output:
    file "intermediate_modified.csv" into output_csv

    """
    mkdir all_results
    ln -s -r csmatch_*/* all_results
    python  $workflow.projectDir/scripts/02_makeCSVTable_csmatch.py -i all_results -o intermediate_modified.csv -f $params.outdir
    """
  }
}
