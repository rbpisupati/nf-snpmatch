# Generate a CSV file based on the ScoreAcc files generated before
# Input:
# 1. Output file
# 2. Folder ID that needs to be added
# 3. expParents


library("optparse")
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
AmbiWindows <- character()
AmbiWindowCount <- character()
TopHitsNumber <- numeric()
for (echscore in wind.filelist){
  name <- sub("[.]windowscore[.]txt", "", echscore)
  likeLiwind <- try(read.table(echscore, header = F, as.is = TRUE))
  ScoreAcc <- try(read.table(paste(name, ".scores.txt", sep =""), header = FALSE, as.is = TRUE))
  if(inherits(likeLiwind, 'try-error') | inherits(ScoreAcc, 'try-error')){
    next
  }
  filename <- unlist(strsplit(name, "[.]"))[1]
#  givparents <- sub("F2_p1_#", "", filename)
#  expP1 <- unlist(strsplit(givparents, "x"))[1]
#  expP2 <- unlist(strsplit(givparents, "x"))[2]

  homowind <- which(likeLiwind$V7 == 1)
  clean <- sort(table(likeLiwind$V1[homowind]))
  clean <- as.data.frame(clean)
  topHits <- rownames(tail(clean, 2))

  if (!is.na(topHits[2])){
    ScoreP2 <- c(ScoreP2, mean(likeLiwind$V4[homowind[which(likeLiwind$V1[homowind] == topHits[1])]], na.rm = TRUE))
    ScoreP1 <- c(ScoreP1, mean(likeLiwind$V4[homowind[which(likeLiwind$V1[homowind] == topHits[2])]], na.rm = TRUE))
    FouParents1 <- c(FouParents1, as.character(topHits[2]))
    FouParents2 <- c(FouParents2, as.character(topHits[1]))
    CountP2 <- c(CountP2, as.numeric(tail(clean, 2)[1,]))
    CountP1 <- c(CountP1, as.numeric(tail(clean, 2)[2,]))
  } else {
    ScoreP1 <- c(ScoreP1, mean(likeLiwind$V4[homowind[which(likeLiwind$V1[homowind] == topHits[1])]], na.rm = TRUE))
    ScoreP2 <- c(ScoreP2, NA)
    FouParents1 <- c(FouParents1, as.character(topHits[1]))
    FouParents2 <- c(FouParents2, NA)
    CountP1 <- c(CountP1, as.numeric(tail(clean, 2)[1,]))
    CountP2 <- c(CountP2, NA)
  }

  ambWind <- which(likeLiwind$V7 < 20)
  nclean <- sort(-table(likeLiwind$V1[ambWind]))
  nclean <- as.data.frame(nclean)
  nclean <- nclean[which(nclean$nclean <= -1),]

  AmbiWindows <- c(AmbiWindows, paste(rownames(nclean), collapse = ":"))
  AmbiWindowCount <- c(AmbiWindowCount, paste(as.numeric(-nclean), collapse = ":"))

  NumWindows <- c(NumWindows, max(likeLiwind$V8))
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

DF <- cbind(PLATE = plate, FILENAME = Names, ExpectedParents = GivParents, TopHit = TopHit, NextHit = NextHit, TopScore = TopScore, FracScore = FracScore, TopHitsNumber = TopHitsNumber, LikelihoodRatio = LikeLihoodTopHit, NextHitLLR = LLRNextHit, SNPsCalled = SNPscalled, ObservedParent1 = FouParents1, ObservedParent2 = FouParents2, CountP1 = CountP1, CountP2 = CountP2, ScoreP1 = ScoreP1, ScoreP2 = ScoreP2, AmbigousWindows = AmbiWindows, AmbigousWindowCount = AmbiWindowCount)

write.csv(DF, file = outFile)
