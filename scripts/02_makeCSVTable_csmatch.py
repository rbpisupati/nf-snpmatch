# Generate a CSV file based on the ScoreAcc files generated before
# Input:
# 1. Output file
# 2. Folder ID that needs to be added
import numpy as np
import pandas as pd
import argparse
import logging
from glob import glob
from snpmatch.core import snpmatch
from snpmatch.core import snp_genotype
import os.path
import json
import re

logging.basicConfig(format='%(levelname)s:%(asctime)s: %(name)s:  %(message)s', level=logging.DEBUG)
log = logging.getLogger(__name__)

inOptions = argparse.ArgumentParser(description='get a genotype matrix')
inOptions.add_argument("-i", dest="files_path", default = "./", help="path to the snpmatch output files")
inOptions.add_argument("--file_spearation", dest="file_sep", default = ".", type = str, help="file name to file id (is it . separated or _)")
inOptions.add_argument("--dirs", action="store_true", dest="dirs", default=False, help="Input files are in folders")
inOptions.add_argument("-d", dest="accDB", help="hdf5 file for database accessions")
inOptions.add_argument("-e", dest="expParents", help="expected parents, ex: 6090x9416")
inOptions.add_argument("-f", dest="folder_id", help="Plate ID to be added in the CSV file")
inOptions.add_argument("-o", dest="output_file", default = "final_csmatch_genotyper", type=str, help="output file")

args = inOptions.parse_args()

log.info("reading input files")
if args.dirs:
    input_files = glob(args.files_path + "/csmatch_*/" + "*.windowscore.txt")
else:
    input_files = glob(args.files_path + "/" + "*.windowscore.txt")

input_files = pd.Series(input_files).astype(str)
assert len(input_files) > 0, "there are  no files with snpmatch scores in the folder."
input_ids = input_files.apply( os.path.basename ).str.split(args.file_sep, expand=True)
if np.unique(input_ids.iloc[:,0]).shape[0] == input_ids.shape[0]:
    input_ids = np.array(input_ids.iloc[:,0])
elif np.unique(input_ids.iloc[:,0] + args.file_sep + input_ids.iloc[:,1]).shape[0] == input_ids.shape[0]:
    input_ids = np.array(input_ids.iloc[:,0] + args.file_sep + input_ids.iloc[:,1])
else:
    input_ids = input_files.apply( os.path.basename ).str.replace( ".windowscore.txt", "" ).values

input_files = input_files.values

main_output_cols = ['FOL', "FILENAME", "ExpectedParents", 'TopHit', 'NextHit', 'Score', 'FracScore', 'SNPsCalled', 'LikelihoodRatio','NextHitLLR', 'TopHitsNumber', "ObservedParent1", "ObservedParent2", "CountP1", "CountP2", 'percent_heterozygosity', 'Overlap', "IdenticalWindows", "NumberofWindows", "HomozygousWindow", "HomozygousWindowCount"]
main_output = pd.DataFrame(columns = main_output_cols, index = input_ids)
main_output['FOL'] = args.folder_id
main_output['ExpectedParents'] = args.expParents

snp_db = snp_genotype.load_genotype_files(args.accDB)
sample_ancestry = pd.DataFrame(index = snp_db.accessions, columns =  input_ids)

assert input_ids.shape == np.array(input_files).shape, "Something is wrong"

for ef_ix in range(len(input_ids)):
    window_data = pd.read_csv(input_files[ef_ix],sep = "\t", dtype = {'acc': str})
    ScoreAcc = pd.read_csv(re.sub("windowscore.txt$", "scores.txt", input_files[ef_ix]) , header = None, sep = "\t")
    ScoreAcc = ScoreAcc.sort_values([5, 3], ascending=[True, False])
    main_output.loc[input_ids[ef_ix], "FILENAME"] = re.sub(".windowscore.txt$", "", input_files[ef_ix])
    main_output.loc[input_ids[ef_ix], "TopHit"] = ScoreAcc.iloc[0,0]
    main_output.loc[input_ids[ef_ix], "NextHit"] = ScoreAcc.iloc[1,0]
    main_output.loc[input_ids[ef_ix], "Score"] = ScoreAcc.iloc[0,3]
    main_output.loc[input_ids[ef_ix], "FracScore"] = ScoreAcc.iloc[1,3] / ScoreAcc.iloc[0,3]
    main_output.loc[input_ids[ef_ix], "SNPsCalled"] = ScoreAcc.iloc[0, 6]
    main_output.loc[input_ids[ef_ix], "LikelihoodRatio"] = ScoreAcc.iloc[0,4]
    main_output.loc[input_ids[ef_ix], "NextHitLLR"] = ScoreAcc.iloc[1,5]
    main_output.loc[input_ids[ef_ix], "TopHitsNumber"] = np.where(ScoreAcc.iloc[:,5] < snpmatch.lr_thres)[0].shape[0]
    num_winds = np.unique(window_data['window_index']).shape[0]
    main_output.loc[input_ids[ef_ix], "NumberofWindows"] = num_winds
    identical_wind = np.where(window_data.groupby('window_index').max()['identical'] == 1)[0].shape[0]
    sample_ancestry[input_ids[ef_ix]] = window_data.groupby('acc').mean()['identical'].reindex(sample_ancestry.index)
    main_output.loc[input_ids[ef_ix], "IdenticalWindows"] = snpmatch.get_fraction( identical_wind, num_winds )
    if not os.path.isfile(re.sub("windowscore.txt$", "matches.json", input_files[ef_ix])):
        continue
    with open(re.sub("windowscore.txt$", "matches.json", input_files[ef_ix])) as json_out:
        jsonstat = json.load(json_out)
    main_output.loc[input_ids[ef_ix], "ObservedParent1"] = jsonstat['parents']['mother'][0]
    main_output.loc[input_ids[ef_ix], "CountP1"] = jsonstat['parents']['mother'][1]
    main_output.loc[input_ids[ef_ix], "ObservedParent2"] = jsonstat['parents']['father'][0]
    main_output.loc[input_ids[ef_ix], "CountP2"] = jsonstat['parents']['father'][1]
    main_output.loc[input_ids[ef_ix], "percent_heterozygosity"] = jsonstat['percent_heterozygosity']
    main_output.loc[input_ids[ef_ix], "Overlap"] = jsonstat['overlap'][0]
    if len(jsonstat['matches']) > 0:
        main_output.loc[input_ids[ef_ix], "HomozygousWindow"] = pd.DataFrame(jsonstat['matches']).iloc[0:20,0].astype(str).str.cat(sep=',')
        main_output.loc[input_ids[ef_ix], "HomozygousWindowCount"] = pd.DataFrame(jsonstat['matches']).iloc[0:20,1].astype(str).str.cat(sep=',')
    if (ef_ix + 1) % 100 == 0:
        log.info("progress: %s files, %s%%" % (ef_ix, (ef_ix * 100 / len(input_ids)) ))


log.info("writing output into %s" % args.output_file )
main_output.to_csv(args.output_file + ".csv")
sample_ancestry.to_csv(args.output_file + "iden_ancestry.csv")
log.info("finished!")
