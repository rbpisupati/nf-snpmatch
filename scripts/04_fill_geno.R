library(qtl)
library("optparse")
options(warn=-1)



option_list = list(
  make_option(c("-i", "--inFile"), type="character", default=NULL, help="input csv file of genotypes per windows", metavar="character"),
  make_option(c("-t", "--error_thres"), type="numeric", default=NULL, help="error rate for windows", metavar="numeric"),
  make_option(c("-o", "--outFile"), type="character", help="output file for filled genotypes", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

input.geno.file <- opt$inFile
output.geno.file <- opt$outFile

geno_t = read.cross("csvr", dirname(input.geno.file), basename(input.geno.file), genotypes=c(0,1,2), na.strings = c('NA', "nan"))

geno_fill = fill.geno(geno_t, error.prob = opt$error_thres, method = "argmax")

print(geno_fill)

write.cross(geno_fill, format = "csvr", output.geno.file)
