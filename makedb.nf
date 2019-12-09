#!/usr/bin/env nextflow
/*
Simply run nextflow pipeline as:

nextflow run makedb.nf --input "*vcf"  --outdir output_folder
*/
/*
 * SET UP CONFIGURATION VARIABLES
*/
params.input = false
params.outdir = 'db' 

//input files
input_files = Channel
    .fromPath ( params.input )
    .map { [ "$it.baseName", file("$it") ] }
    .ifEmpty { exit 1, "Cannot find any input files matching: ${params.input}\nNB: Path needs to be enclosed in quotes!\n" }


process makeDB {
  tag { "${prefix}" }
  publishDir "${params.outdir}", mode:"copy"

  input:
  set val(prefix), file(input_file) from input_files

  output:
  set val(prefix), file("${prefix}*.hdf5") into outdb

  script:
  """
  snpmatch makedb -v -i $input_file -o $prefix
  """
}

