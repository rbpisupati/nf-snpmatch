#!/usr/bin/env nextflow
/*
Simply run nextflow pipeline as:

nextflow run submit.snpmatch.nf --input "*vcf" --db hdf5_file --db_acc hdf5_acc_file --outdir output_folder
*/
/*
 * SET UP CONFIGURATION VARIABLES
*/
params.input = false
params.parents=false
params.windows=300000
params.outdir = 'snpmatch_1135g'

params.project = "the1001genomes"
// databases
params.db = "/lustre/scratch/projects/the1001genomes/rahul/101.VCF_1001G_1135/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.hdf5"
params.db_acc= "/lustre/scratch/projects/the1001genomes/rahul/101.VCF_1001G_1135/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.acc.hdf5"

//input files
input_files = Channel
    .fromPath ( params.input )
    .map { [ "$it.baseName", "${it.getParent()}", file("$it") ] }
    .ifEmpty { exit 1, "Cannot find any input files matching: ${params.input}\nNB: Path needs to be enclosed in quotes!\n" }

db_file = Channel
    .fromPath ( params.db )
    .ifEmpty { exit 1, "please provide hdf5 file as a database" }

db_acc_file = Channel
    .fromPath ( params.db_acc )
    .ifEmpty { exit 1, "please provide hdf5 acc file as a database" }


process parse_inputfiles {
  tag { "${prefix}" }
  storeDir "${input_folder}"

  input:
  set val(prefix), val(input_folder), file(input_file) from input_files

  output:
  set val(prefix), file("${input_file}.snpmatch.npz") into input_gzips
  file "${input_file}.snpmatch.stats.json" into input_stats

  script:
  """
  snpmatch parser -v -i $input_file
  """
}

input_files_dbs = input_gzips.combine(db_file).combine(db_acc_file)

process genotype_cross {
  tag { "${prefix}" }
  publishDir "$params.outdir", mode: 'copy'
  errorStrategy { task.exitStatus in [143,137] ? 'retry' : 'ignore' }

  input:
  set val(prefix), file(input_npz), file(f_db), file(f_db_acc) from input_files_dbs

  output:
  file "${prefix}.snpmatch*" into snpmatch_output

  script:
  """
  snpmatch genotype_cross -v -e $f_db_acc -i $input_npz -p "$params.parents" -b "$params.windows" -o ${prefix}.genotyper.txt
  """
}
