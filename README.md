# Nextflow pipeline for running SNPmatch

## Installation

This pipeline itself needs no installation - NextFlow will automatically fetch it from GitHub rbpisupati/nf-snpmatch

```bash
git clone https://github.com/rbpisupati/nf-snpmatch.git
```

You need to have SNPmatch installed in your system. Please visit (https://github.com/Gregor-Mendel-Institute/SNPmatch.git) for detailed instructions on installing SNPmatch.

## Running the pipeline

### Genotyping inbred lines
```bash
nextflow run main.nf --func "inbred" --input "*.vcf" --outdir "snpmatch" --db "hdf5" --db_acc "hdf5_acc"
## or run cross genotyper
nextflow run main.nf --func "cross" --input "*.vcf" --outdir "snpmatch" --db "hdf5" --db_acc "hdf5_acc"
```

### Genotyping crosses in windows given parents
```bash
nextflow run genocross.nf --input "*.vcf" --parents "6191x6046" --window 300000 --outdir "snpmatch" --db "hdf5_path" --db_acc "hdf5_acc_path"
```
Input, parents, outdir are required arguments.


## Configuration

The pipeline is written mainly to run SNPmatch on GMI HPC mendel which is PBS system. Please change config file accordingly to run it on your system.

## Credits

- Rahul Pisupati (rahul.pisupati[at]gmi.oeaw.ac.at)

## Citation
Cite the paper below if you use SNPmatch tool.
Pisupati, R. *et al.*. Verification of *Arabidopsis* stock collections using SNPmatch, a tool for genotyping high-plexed samples.  *Nature Scientific Data*  **4**, 170184 (2017).
[doi:10.1038/sdata.2017.184](https://www.nature.com/articles/sdata2017184)
