# Generate a CSV file based on the ScoreAcc files generated before
# Input:
# 1. Output file
# 2. Folder ID that needs to be added
# 3. expParents


library("optparse")
library("jsonlite")

option_list = list(
  make_option(c("-o", "--outFile"), type="character", default="intermediate_modified.csv", help="Output file name [default= %default]", metavar="character"),
  make_option(c("-f", "--folID"), type="character", default=NULL, help="Plate ID to be added in the CSV file", metavar="character"),
  make_option(c("-e", "--expParents"), type="character", default=".x.", help="expected parents, ex: 6090x9416", metavar="character")
);
opt = parse_args(OptionParser(option_list=option_list));

outFile = opt$outFile
id = opt$folID
expParents = opt$expParents

if(file.exists(outFile)){
  file.remove(outFile)
}

wind.filelist <- list.files("./", pattern = "[.]windowscore.txt$")
Names <- character()
plate <- character()
GivParents <- character()
FouParents1 <- character()
FouParents2 <- character()
TopHit <- character()
TopScore <- numeric()
NextHit <- character()
FracScore <- numeric()
LikeLihoodTopHit <- numeric()
LLRNextHit <- numeric()
TopHitsNumber <- numeric()
ScoreP1 <- numeric()
ScoreP2 <- numeric()
NumWindows <- numeric()
CountP1 <- numeric()
CountP2 <- numeric()
SNPscalled <- numeric()
PerHeterozygosity <- numeric()
Overlap <- numeric()
AmbiWindows <- character()
AmbiWindowCount <- character()
TopHitsNumber <- numeric()
for (echscore in wind.filelist){
  filename <- sub("[.]windowscore[.]txt", "", echscore)
  name <- sub("[.]csmatch", "", sub("[.]filter", "", sub("[.]snpmatch", "", filename)))
  likeLiwind <- try(read.table(echscore, header = F, as.is = TRUE))
  ScoreAcc <- try(read.table(paste(filename, ".scores.txt", sep =""), header = FALSE, as.is = TRUE))
  jsonstat_per <- try(fromJSON(paste(filename, ".scores.txt.matches.json", sep = "")))
  if(inherits(likeLiwind, 'try-error') | inherits(ScoreAcc, 'try-error') |  inherits(jsonstat_per, 'try-error')){
    next
  }
  
  plate <- c(plate, id)
  ranks <- order(-ScoreAcc$V4, na.last = TRUE, partial = ScoreAcc$V6)
  topscore <- ScoreAcc$V4[ranks[1]]
  topacc <- as.character(ScoreAcc$V1[ranks[1]])
  nextacc <- as.character(ScoreAcc$V1[ranks[2]])
  nextscore <- ScoreAcc$V4[ranks[2]]
  frac <- nextscore/topscore
  snps <- ScoreAcc$V7[ranks[1]]
  topnum <- length(which(ScoreAcc$V6 < 3.841))
  maxlike <- ScoreAcc$V5[ranks[1]]
  nextlike <- ScoreAcc$V6[ranks[2]]
  PerHeterozygosity = c(PerHeterozygosity, jsonstat_per$percent_heterozygosity)
  Overlap = c(Overlap, jsonstat_per$overlap[1])

  if (file.exists( paste(filename, ".matches.json", sep = "") )){
    jsonstat <- try(fromJSON(paste(filename, ".matches.json", sep = "")))
    FouParents1 <- c(FouParents1, as.character( jsonstat$parents$mother[1] ))
    FouParents2 <- c(FouParents2, as.character( jsonstat$parents$father[1] ))
    CountP1 <- c(CountP1, as.numeric( jsonstat$parents$mother[2] ))
    CountP2 <- c(CountP2, as.numeric( jsonstat$parents$father[2] ))
  } else {
    FouParents1 <- c(FouParents1, "")
    FouParents2 <- c(FouParents2, "")
    CountP1 <- c(CountP1, 0)
    CountP2 <- c(CountP2, 0)
  }

  ambWind <- which(likeLiwind$V7 < 20)
  nclean <- sort(-table(likeLiwind$V1[ambWind]))
  nclean <- as.data.frame(nclean)
  nclean <- nclean[which(nclean$nclean <= -1),]
  AmbiWindows <- c(AmbiWindows, paste(rownames(nclean), collapse = ":"))
  AmbiWindowCount <- c(AmbiWindowCount, paste(as.numeric(-nclean), collapse = ":"))
  NumWindows <- c(NumWindows, max(likeLiwind$V8))

  LikeLihoodTopHit <- c(LikeLihoodTopHit, maxlike)
  LLRNextHit <- c(LLRNextHit, nextlike)
  SNPscalled <- c(SNPscalled, snps)
  TopHit <- c(TopHit, topacc)
  NextHit <- c(NextHit, nextacc)
  TopScore <- c(TopScore, topscore)
  FracScore <- c(FracScore, frac)
  topnum <- length(which(ScoreAcc$V6 < 3.841))
  TopHitsNumber <- c(TopHitsNumber, topnum)
  Names <- c(Names, name)
  GivParents <- c(GivParents, expParents)
}

DF <- cbind(PLATE = plate, FILENAME = Names, ExpectedParents = GivParents, TopHit = TopHit, NextHit = NextHit, TopScore = TopScore, FracScore = FracScore, TopHitsNumber = TopHitsNumber, LikelihoodRatio = LikeLihoodTopHit, NextHitLLR = LLRNextHit, SNPsCalled = SNPscalled, ObservedParent1 = FouParents1, ObservedParent2 = FouParents2, CountP1 = CountP1, CountP2 = CountP2, ScoreP1 = ScoreP1, ScoreP2 = ScoreP2, Overlap = Overlap, Heterozygosity = PerHeterozygosity, AmbigousWindows = AmbiWindows, AmbigousWindowCount = AmbiWindowCount)

warnings()

write.csv(DF, file = outFile)
