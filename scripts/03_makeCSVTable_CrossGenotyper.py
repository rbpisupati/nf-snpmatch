# Generate a CSV file for the cross genotyper
# Input:
# 1. Output files (genotyper.txt)
#!/usr/bin/python
#----------------
## given list of genes, two accessions, identify genes that are differentially methylated
# Written by Rahul Pisupati 04.07.2018
#
#----------------
import pandas as pd
import numpy as np
import argparse
import logging
from snpmatch.core import csmatch
from bshap.core import genome, the1001g
tair10 = genome.ArabidopsisGenome()
from glob import glob
import os.path

logging.basicConfig(format='%(levelname)s:%(asctime)s: %(name)s:  %(message)s', level=logging.DEBUG)
log = logging.getLogger(__name__)

inOptions = argparse.ArgumentParser(description='get a genotype matrix')
inOptions.add_argument("-i", dest="files_path", default = ".", help="path for input files")
inOptions.add_argument("-b", dest="binLen", default = 300000, type=int, help="length for windows")
inOptions.add_argument("-o", dest="out_prefix", default = "genotyper", help="output file prefix")

args = inOptions.parse_args()

log.info("reading input files")
input_files = glob(args.files_path + "/" + "*.genotyper.txt")
input_ids = [ os.path.basename(efile).split("_")[0] for efile in input_files]
input_ids = np.array(input_ids)
redundant_inputs = np.max(np.unique(input_ids, return_counts=True)[1])

log.info("done! %s files" % len( input_files ))

log.info("iterating over genome")
iter_winds = tair10.iter_windows(args.binLen)
iter_winds_str = []
for ei in iter_winds:
    iter_winds_str.append( "%s\t%s\t%s" % (ei[0], ei[1], ei[2]) )
log.info("done")

all_genotypes = np.zeros((len(input_files), len(iter_winds_str)), dtype="string" )
all_genotypes_ids = np.zeros(0, dtype="string")
for efile_ix in range(len(input_files)):
    epd = pd.read_table(input_files[efile_ix], header = None)
    e_id = os.path.basename(input_files[efile_ix]).split("_")
    if redundant_inputs > 1:
        e_id = e_id[0] + e_id[1]
    else:
        e_id = e_id[0]
    e_out = open( args.out_prefix + "." + e_id + ".txt", 'w' )
    for ewind_ix in range(len(iter_winds_str)):
        t_gen = epd.iloc[:,3][np.where( epd.iloc[:,0] == ewind_ix + 1 )[0]]
        if t_gen.shape[0] > 0 and t_gen.dropna().shape[0] > 0:
            t_gen = str(int(float(t_gen) * 2))
            e_out.write("%s\t%s\n" % ( iter_winds_str[ewind_ix], t_gen ))
            all_genotypes[efile_ix, ewind_ix] = t_gen
        else:
            e_out.write("%s\tnan\n" % ( iter_winds_str[ewind_ix] ))
            all_genotypes[efile_ix, ewind_ix] = '-'
    all_genotypes_ids = np.append(all_genotypes_ids, e_id)
    e_out.close()

the1001g.generate_h5_1001g(args['file_paths'], args['output_file'])
#all_genotypes = pd.DataFrame(all_genotypes, )
#all_genotypes.insert(0, column = "pheno", value = np.zeros( all_genotypes.shape[0] ) )
#all_genotypes = pd.concat([pd.DataFrame(np.zeros( all_genotypes.columns.values.shape[0], dtype="string"), index= all_genotypes.columns.values).T, all_genotypes], ignore_index=True)
