nextflow.enable.dsl=2

params.skip_db_hets = false
params.genome = "athaliana_tair10"


workflow SNPMATCH_INBRED {
    take:
    ch_input_vcf
    input_hdf5
    output_folder

    main:

    // ch_input_vcf = Channel
    //     .fromPath ( input_vcf )
    //     .map { [ "$it.baseName", "${it.getParent()}", file("$it") ] }
    //     .ifEmpty { exit 1, "Cannot find any input files matching: ${input_vcf}\nNB: Path needs to be enclosed in quotes!\n" }

    ch_input_hdf5 = Channel
        .fromPath ( input_hdf5 )
        .map{ [ file("${it}"), file("${it.parent}/${it.baseName}.acc.hdf5") ] }
        .ifEmpty { exit 1, "please provide hdf5 files as a database" }
    

    identify_libs = identify_libraries( ch_input_vcf, ch_input_hdf5.collect(), output_folder )
    out_csv = make_csv_inbred( identify_libs.snpmatch_out_for_csv.collect(), output_folder )


    emit:
    output_csv = out_csv
}


process identify_libraries {
    tag { "${sample_id}" }
    conda "/users/rahul.pisupati/.conda/envs/snpmatch-5.0/"


    publishDir "$outdir", mode: 'copy'
    errorStrategy { task.exitStatus in [143,137] ? 'retry' : 'ignore' }

    input:
    tuple val(sample_id), path(vcf), path(vcf_idx)
    tuple path(f_db), path(f_db_acc)
    val(outdir)

    output:
    tuple val(sample_id), path("snpmatch_${sample_id}"), emit: snpmatch_out
    path("snpmatch_${sample_id}"), emit: snpmatch_out_for_csv

    script:
    skip_db_hets = params.skip_db_hets != false ? "--skip_db_hets" : ''
    """
    mkdir -p snpmatch_${sample_id}
    snpmatch inbred -v  $skip_db_hets -d $f_db -e $f_db_acc -i $vcf -o snpmatch_${sample_id}/${sample_id} --refine
    """
}


process make_csv_inbred {
    conda "/users/rahul.pisupati/.conda/envs/snpmatch-5.0/"

    publishDir "$output_id", mode: 'copy'

    input:
    path all_results
    val(output_id)

    output:
    path "intermediate_modified.csv"

    """
    01_makeCSVTable_inbred.py --dirs -i ./ -o intermediate_modified.csv -f $output_id
    """
}

