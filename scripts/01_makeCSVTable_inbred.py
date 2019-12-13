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
import os.path
import json
import re


logging.basicConfig(format='%(levelname)s:%(asctime)s: %(name)s:  %(message)s', level=logging.DEBUG)
log = logging.getLogger(__name__)

inOptions = argparse.ArgumentParser(description='get a genotype matrix')
inOptions.add_argument("-i", dest="files_path", default = "./", help="path to the snpmatch output files")
inOptions.add_argument("--file_spearation", dest="file_sep", default = ".", type = str, help="file name to file id (is it . separated or _)")
inOptions.add_argument("-f", dest="folder_id", help="Plate ID to be added in the CSV file")
inOptions.add_argument("-o", dest="output_file", default = "final_snpmatch_genotyper.csv", type=str, help="output file")

args = inOptions.parse_args()

log.info("reading input files")
input_files = glob(args.files_path + "/" + "*.scores.txt")
input_files = pd.Series(input_files).astype(str)
input_files = input_files[~input_files.str.contains('refined.scores.txt',  regex=True)]
input_files = np.array(input_files)
assert len(input_files) > 0, "there are  no files with snpmatch scores in the folder."
input_ids = pd.Series([ os.path.basename(efile) for efile in input_files]).str.split(args.file_sep, expand=True)
if np.unique(input_ids.iloc[:,0]).shape[0] == input_ids.shape[0]:
    input_ids = np.array(input_ids.iloc[:,0])
elif np.unique(input_ids.iloc[:,0] + args.file_sep + input_ids.iloc[:,1]).shape[0] == input_ids.shape[0]:
    input_ids = np.array(input_ids.iloc[:,0] + args.file_sep + input_ids.iloc[:,1])


main_output_cols = ['FOL', "FILENAME", 'TopHitAccession', 'NextHit', 'ThirdHit','Score', 'FracScore', 'SNPsinfoAcc', 'SNPsCalled', 'MeanDepth', 'LikelihoodRatio','NextHitLLR', 'percent_heterozygosity', 'Overlap', 'TopHitsNumber','TopHits', "RefinedTopHit", "RefinedScore", "RefinedSNPs", "RefinedTopHitNumber"]
main_output = pd.DataFrame(columns = main_output_cols, index = input_ids)
main_output['FOL'] = args.folder_id
assert input_ids.shape == np.array(input_files).shape, "cannot determine ids for %s files, try changing --file_separation" % (input_ids.shape[0], np.array(input_files).shape[0] )

for ef_ix in range(len(input_ids)):
    ScoreAcc = pd.read_csv(input_files[ef_ix], header = None, sep = "\t")
    ScoreAcc = ScoreAcc.sort_values([5, 3], ascending=[True, False])
    if not os.path.isfile(re.sub("scores.txt$", "matches.json", input_files[ef_ix])):
        log.info("skipping: missing a json file for %s" % input_files[ef_ix])
        continue
    with open(re.sub("scores.txt$", "matches.json", input_files[ef_ix])) as json_out:
            jsonstat = json.load(json_out)
    main_output.loc[input_ids[ef_ix], "FILENAME"] = re.sub(".scores.txt$", "", input_files[ef_ix])
    main_output.loc[input_ids[ef_ix], "TopHitAccession"] = ScoreAcc.iloc[0,0]
    main_output.loc[input_ids[ef_ix], "NextHit"] = ScoreAcc.iloc[1,0]
    main_output.loc[input_ids[ef_ix], "ThirdHit"] = ScoreAcc.iloc[2,0]
    main_output.loc[input_ids[ef_ix], "Score"] = ScoreAcc.iloc[0,3]
    main_output.loc[input_ids[ef_ix], "FracScore"] = ScoreAcc.iloc[1,3] / ScoreAcc.iloc[0,3]
    main_output.loc[input_ids[ef_ix], "SNPsinfoAcc"] = ScoreAcc.iloc[0,2]
    main_output.loc[input_ids[ef_ix], "SNPsCalled"] = jsonstat['overlap'][1]
    main_output.loc[input_ids[ef_ix], "MeanDepth"] = ScoreAcc.iloc[0,7]
    main_output.loc[input_ids[ef_ix], "LikelihoodRatio"] = ScoreAcc.iloc[0,4]
    main_output.loc[input_ids[ef_ix], "NextHitLLR"] = ScoreAcc.iloc[1,5]
    main_output.loc[input_ids[ef_ix], "percent_heterozygosity"] = jsonstat['percent_heterozygosity']
    main_output.loc[input_ids[ef_ix], "Overlap"] = jsonstat['overlap'][0]
    main_output.loc[input_ids[ef_ix], "TopHitsNumber"] = len(jsonstat['matches'])
    if len(jsonstat['matches']) == 1:
        main_output.loc[input_ids[ef_ix], "TopHits"] = ""
    elif len(jsonstat['matches']) <= 10:
        main_output.loc[input_ids[ef_ix], "TopHits"] = pd.DataFrame(jsonstat['matches']).iloc[:,0].str.cat(sep=',')
    else:
        main_output.loc[input_ids[ef_ix], "TopHits"] = "TooManytoPrint"
    if os.path.isfile(re.sub("scores.txt$", "refined.scores.txt", input_files[ef_ix])):
        refined_scores = pd.read_csv(re.sub("scores.txt$", "refined.scores.txt", input_files[ef_ix]), header = None, sep = "\t")
        refined_scores = refined_scores.sort_values([5, 3], ascending=[True, False])
        main_output.loc[input_ids[ef_ix], "RefinedTopHit"] = str(refined_scores.iloc[0,0])
        main_output.loc[input_ids[ef_ix], "RefinedScore"] = refined_scores.iloc[0,3]
        main_output.loc[input_ids[ef_ix], "RefinedSNPs"] = int(refined_scores.iloc[0,2])
        main_output.loc[input_ids[ef_ix], "RefinedTopHitNumber"] = int(np.where(refined_scores.iloc[:,5] <= snpmatch.lr_thres)[0].shape[0])
    if (ef_ix + 1) % 100 == 0:
        log.info("progress: %s files, %s%%" % (ef_ix, (ef_ix * 100 / len(input_ids)) ))


log.info("writing output into %s" % args.output_file )
main_output.to_csv(args.output_file, sep = "\t")
log.info("finished!")
#
#
# write.csv(DF, file = outFile)
