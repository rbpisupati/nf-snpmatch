# Generate a CSV file based on the ScoreAcc files generated before
# Input:
# 1. Output file
# 2. Folder ID that needs to be added

library("optparse")
library("jsonlite")
option_list = list(
  make_option(c("-o", "--outFile"), type="character", default="intermediate_modified.csv", help="Output file name [default= %default]", metavar="character")
);
opt = parse_args(OptionParser(option_list=option_list));

outFile = opt$outFile

#-------
allScoreFiles <- list.files("./", pattern = "[.]scores.txt$")
AssignedAcc <- character()
TopHitAcc <- character()
NextHitAcc <- character()
ThirdHit <- character()
TopHitScore <- numeric()
fol <- character()
FracScore <- numeric()
LikeLihoodTopHit <- numeric()
LLRNextHit <- numeric()
TopHitsNumber <- numeric()
TopHits <- character()
ChoiceAcc <- numeric()
for (file in allScoreFiles){
  ScoreAcc <- read.table(file, header = F)
  name <- sub(".snpmatch.scores.txt","",file)
  if(inherits(ScoreAcc, 'try-error')){
    next
  }

  folid <- unlist(strsplit(name, "_"))[2]
  assacc <- unlist(strsplit(name, "_"))[1]
  AssignedAcc <- c(AssignedAcc, assacc)
## Changed the sorting order based on the score to the likelihood ratio
#  ranks <- order(-ScoreAcc$V4, na.last = TRUE, partial = ScoreAcc$V6)
  ranks <- order(ScoreAcc$V6, na.last = TRUE, partial = ScoreAcc$V4)

  topscore <- ScoreAcc$V4[ranks[1]]
  topacc <- as.character(ScoreAcc$V1[ranks[1]])
  maxlike <- ScoreAcc$V5[ranks[1]]

#  snps <- as.numeric(numSNPs$V1[which(numSNPs$V2 == name)])
  nextacc <- as.character(ScoreAcc$V1[ranks[2]])
  nextscore <- ScoreAcc$V4[ranks[2]]
  nextlike <- ScoreAcc$V6[ranks[2]]
  frac <- nextscore/topscore

  thirdHit <- as.character(ScoreAcc$V1[ranks[3]])
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

  choiceacc <- which(ScoreAcc$V1[ranks] == assacc)
  ChoiceAcc <- c(ChoiceAcc, choiceacc)

  fol <- c(fol, folid)
  FracScore <- c(FracScore, frac)
  TopHitAcc <- c(TopHitAcc, topacc)
  NextHitAcc <- c(NextHitAcc, nextacc)
  ThirdHit <- c(ThirdHit, thirdHit)
  TopHitScore <- c(TopHitScore, topscore)
  LikeLihoodTopHit <- c(LikeLihoodTopHit, maxlike)
  LLRNextHit <- c(LLRNextHit, nextlike)
}

DF <- cbind(FOL = fol, AssignedAcc = AssignedAcc, TopHitAccession = TopHitAcc, NextHit = NextHitAcc, ThirdHit = ThirdHit, Score = as.numeric(TopHitScore), FracScore = as.numeric(FracScore), LikelihoodRatio = LikeLihoodTopHit, NextHitLLR = LLRNextHit, TopHitsNumber = TopHitsNumber, ChoiceAcc = ChoiceAcc, TopHits = TopHits)

#----

write.csv(DF, file = outFile, quote = F)
