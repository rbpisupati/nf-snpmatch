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
inOptions.add_argument("-o", dest="output_file", default = "genotyper.csv", type=str, help="output file prefix")

args = inOptions.parse_args()

log.info("reading input files")
input_files = glob(args.files_path + "/" + "*.genotyper.txt")
input_ids = pd.Series([ os.path.basename(efile) for efile in input_files]).str.split("_", expand=True)

if len(np.unique(input_ids.iloc[:,0])) != len(input_files):
    input_ids = np.array(input_ids.iloc[:,0] + input_ids.iloc[:,1], dtype="string")
else:
    input_ids = np.array(input_ids.iloc[:,0], dtype="string")
log.info("done! %s files" % len( input_files ))


base_epd = pd.read_table(input_files[0], header = None)

all_genotypes = np.zeros((base_epd.shape[0], len(input_files)), dtype="float" )
all_genotypes[:,0] = base_epd.iloc[:,5]
for efile_ix in range(len(input_files) - 1 ):
    epd = pd.read_table(input_files[efile_ix + 1], header = None)
    if base_epd.shape[0] != epd.shape[0]:
        log.error("Lines are not same in all the files, please check")
    all_genotypes[:,efile_ix + 1] = epd.iloc[:,5]


est_cm = base_epd.iloc[:,0].astype(str) + "," + base_epd.iloc[:,1].astype(str) + "," + base_epd.iloc[:,2].astype(str)
est_cm = est_cm.apply(tair10.estimated_cM_distance)
est_cm = est_cm - est_cm[0]

log.info("writing csvr file")
out_csvr = open( args.output_file, 'w')
out_csvr.write( "%s,%s\n" % ('id,,', str(pd.Series(input_ids).str.cat(sep=","))))
out_csvr.write( "%s,%s\n" % ('pheno,', ',0' * len(input_ids)))

for em_ix in range(base_epd.shape[0]):
    echr_ix = tair10.get_chr_ind(base_epd.iloc[em_ix,0])
    em_str = str( echr_ix + 1 ) + ":" + str( base_epd.iloc[em_ix,1] ) + "-" + str(base_epd.iloc[em_ix,2]) + "," + str(echr_ix + 1) + "," + str(est_cm[em_ix])
    em_geno_str = pd.Series([ '%.0f' % ef  for ef in all_genotypes[0,:]]).astype(str).str.cat(sep=",")
    out_csvr.write( "%s,%s\n" % ( em_str, em_geno_str ) )

out_csvr.close()
log.info("done!")
