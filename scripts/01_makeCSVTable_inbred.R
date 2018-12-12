# Generate a CSV file based on the ScoreAcc files generated before
# Input:
# 1. Output file
# 2. Folder ID that needs to be added

library("optparse")
library("jsonlite")
option_list = list(
  make_option(c("-o", "--outFile"), type="character", default="intermediate_modified.csv", help="Output file name [default= %default]", metavar="character"),
  make_option(c("-f", "--folID"), type="character", default=NULL, help="Plate ID to be added in the CSV file", metavar="character")
);
opt = parse_args(OptionParser(option_list=option_list));

outFile = opt$outFile
folid = opt$folID

# Correction for changing the rows to coloumns
rows <- c("A","B","C","D","E","F","G","H")
cols <- c("1","2","3","4","5","6","7","8","9","10","11","12")
if(file.exists(outFile)){
  file.remove(outFile)
}

#-------
allScoreFiles <- list.files("./", pattern = "[.]scores.txt$")
Names <- character()
sampleRows <- character()
sampleCol <- character()
TopHitAcc <- character()
NextHitAcc <- character()
ThirdHit <- character()
TopHitScore <- numeric()
FracScore <- numeric()
TopHitAccSNPs <- numeric()
TopHitMatchedSNPs <- numeric()
NextHitMatchedSNPs <- numeric()
MeanDepth <- numeric()
SNPscalled <- numeric()
AFfreq <- numeric()
LikeLihoodTopHit <- numeric()
LLRNextHit <- numeric()
TopHitsNumber <- numeric()
PerHeterozygosity <- numeric()
TopHits <- character()
ChoiceAcc <- numeric()
Overlap <- numeric()
for (file in allScoreFiles){
  ScoreAcc <- read.table(file, header = F)
  name <- sub(".scores.txt","",file)
  jsonstat <- try(fromJSON(paste(name, ".matches.json", sep = "")))
  if(inherits(jsonstat, 'try-error')){
    next
  }
  name <- sub(".snpmatch", "", name)
  snps <- as.integer(jsonstat$overlap[2])
  overlap <- jsonstat$overlap[1]
  casenu <- jsonstat$interpretation$case
  interpreter <- jsonstat$interpretation$text
  perhet = jsonstat$percent_heterozygosity

  #Converting the name of file into rows and columns
  nums <- sub(folid, "", name)
  tempnum <- as.numeric(strsplit(nums, "_")[[1]][1])-700
  #tempidnum <- as.numeric(sub("50","",strsplit(name, "_")[[1]][2]))
  tempid <- chartr("12345678","ABCDEFGH",sub("50","",strsplit(nums, "_")[[1]][2]))
  sampleRows <- c(sampleRows, tempid)
  sampleCol <- c(sampleCol, tempnum)

## Changed the sorting order based on the score to the likelihood ratio
#  ranks <- order(-ScoreAcc$V4, na.last = TRUE, partial = ScoreAcc$V6)
  ranks <- order(ScoreAcc$V6, na.last = TRUE, partial = ScoreAcc$V4)

  topscore <- ScoreAcc$V4[ranks[1]]
  topacc <- as.character(ScoreAcc$V1[ranks[1]])
  snps <- ScoreAcc$V7[ranks[1]]
  depth <- ScoreAcc$V8[ranks[1]]
  maxlike <- ScoreAcc$V5[ranks[1]]

  newLike <- ScoreAcc$V5[ranks]/ScoreAcc$V5[ranks[1]]
  nextlike <- newLike[2]

  nextacc <- as.character(ScoreAcc$V1[ranks[2]])
  nextscore <- ScoreAcc$V4[ranks[2]]
  frac <- nextscore/topscore

  topaccsnps <- ScoreAcc$V3[ranks[1]]
  topmatchsnps <- ScoreAcc$V2[ranks[1]]
  nexthitmatchsnps <- ScoreAcc$V2[ranks[2]]
  Names <- c(Names, name)

  thirdHit <- as.character(ScoreAcc$V1[ranks[3]])
#   if(file.exists(paste(name, ".refScore.txt", sep = ""))){
#      if(length(readLines(paste(name, ".refScore.txt", sep = ""))) > 0){
#        refScore <- read.table(paste(name,".refScore.txt", sep = ""), header = F)
#        refRanks <- order(refScore$V5)
#        topacc <- refScore$V1[refRanks[1]]
#        nextacc <- refScore$V1[refRanks[2]]
#        if(length(which(refScore$V1[refRanks] == accass))){
#          choicenum <- which(refScore$V1[refRanks] == accass)
#        }
#     }
#   }
  topnum <- length(which(ScoreAcc$V6 < 3.841))
  TopHitsNumber <- c(TopHitsNumber, topnum)
  if(topnum > 20){
    TopHits <- c(TopHits, "TooManytoPrint")
  } else if(topnum > 3) {
    alltops <- paste(ScoreAcc$V1[ranks[4:topnum]], collapse = ":")
    TopHits <- c(TopHits, alltops)
  } else {
    TopHits <- c(TopHits, "NA")
  }

#  choiceacc <- which(ScoreAcc$V1[ranks] == name)
#  if(length(choiceacc)){
#    ChoiceAcc <- c(ChoiceAcc, choiceacc[1])
#  } else {
#    ChoiceAcc <- c(ChoiceAcc, -1)
#  }

  Overlap <- c(Overlap, overlap)
  FracScore <- c(FracScore, frac)
  TopHitAcc <- c(TopHitAcc, topacc)
  NextHitAcc <- c(NextHitAcc, nextacc)
  ThirdHit <- c(ThirdHit, thirdHit)
  PerHeterozygosity = c(PerHeterozygosity, perhet)
  TopHitScore <- c(TopHitScore, topscore)
  TopHitAccSNPs <- c(TopHitAccSNPs, topaccsnps)
  TopHitMatchedSNPs <- c(TopHitMatchedSNPs, topmatchsnps)
  NextHitMatchedSNPs <- c(NextHitMatchedSNPs, nexthitmatchsnps)
  SNPscalled <- c(SNPscalled, snps)
  MeanDepth <- c(MeanDepth, depth)
  LikeLihoodTopHit <- c(LikeLihoodTopHit, maxlike)
  LLRNextHit <- c(LLRNextHit, nextlike)
}

fol = rep(folid, length(Names))

DF <- cbind(FOL = fol, FILENAME = Names, TopHitAccession = TopHitAcc, NextHit = NextHitAcc, ThirdHit = ThirdHit, Score = as.numeric(TopHitScore), FracScore = as.numeric(FracScore), MatchedSNPs = as.numeric(TopHitMatchedSNPs), NextHitMatchedSNPs = as.numeric(NextHitMatchedSNPs), SNPsinfoAcc  = as.numeric(TopHitAccSNPs), SNPsCalled = as.numeric(SNPscalled), MeanDepth = MeanDepth, LikelihoodRatio = LikeLihoodTopHit, NextHitLLR = LLRNextHit, percent_heterozygosity = PerHeterozygosity, Overlap = Overlap, TopHitsNumber = TopHitsNumber, TopHits = TopHits)

summary(warnings())

write.csv(DF, file = outFile)
