#!/usr/bin/env nextflow
/*
Simply run nextflow pipeline as:

nextflow run genocross.nf --input "*vcf" --parents "5856x9452" --db hdf5_file --db_acc hdf5_acc_file --outdir output_folder
*/
/*
 * SET UP CONFIGURATION VARIABLES
*/
params.input = false
params.parents=false
params.windows=300000
params.outdir = 'genotype_cross'

params.project = "the1001genomes"
// databases
params.db = "/lustre/scratch/projects/the1001genomes/rahul/107.VCF_1001G_imputed/the1001genomes_all_chromosomes_binary.hdf5"
params.db_acc= "/lustre/scratch/projects/the1001genomes/rahul/107.VCF_1001G_imputed/the1001genomes_all_chromosomes_binary.acc.hdf5"

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
  file "${prefix}.genotyper.txt" into snpmatch_output

  script:
  """
  snpmatch genotype_cross -v -e $f_db_acc -i $input_npz -p "$params.parents" -b "$params.windows" -o ${prefix}.genotyper.txt --hmm
  """
}

snpmatch_output
    .collect()
    .into{ input_csv; input_hdf5 }

process genotyper_csv {
  publishDir "$params.outdir", mode: 'copy'

  input:
  file "*" from input_csv

  output:
  file "genotyper.csv" into output_table

  """
  python $workflow.projectDir/scripts/03_makeCSVTable_CrossGenotyper.py -b $params.windows -o genotyper.csv -i ./
  """
}

process genotyper_hdf5 {
  publishDir "$params.outdir", mode: 'copy'

  input:
  file "*" from input_hdf5

  output:
  file "genotyper.hdf5" into output_h5py

  """
  bshap generate_h5_1001g -i ./*genotyper.txt -d 6 -o genotyper.hdf5 -v
  """
}
