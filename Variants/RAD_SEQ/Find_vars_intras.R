library(zoo)
fm <- read.table("fleximat.txt", sep="\t", header=TRUE)
sm <- read.table("strictmat.txt", sep="\t", header=TRUE)

prange <- 1:750

p3 <- c(30,31,32,33,34,35,36,37,38,39,41,42,43,44,45)
p3mat <- p3 + 2
p1 <- c(1,4,5,6,7,8,15,16,20,21,22,24,27,28,29)
p1mat <- p1 + 2
p2 <- c(2,3,10,11,13,14,17,18,19,23,25,26)
p2mat <- p2 + 2

#correcting the scaffold IDS & loci
scs <- read.table("sclist.txt", sep=" ", header=FALSE)
loci <- read.table("loci.txt", sep=" ", header=FALSE)
fm[,1] <- scs
sm[,1] <- scs
fm[,2] <- loci
sm[,2] <- loci

bins <- read.table("../bins.txt", sep="\t", header=FALSE)
ravg <- rollmean(bins[,2], 500, fill=0.4, na.pad = TRUE)

flex_ids <- read.table("flex_locs2.txt", sep="\t", header=FALSE)
med_ids <- read.table("med_locs2.txt", sep="\t", header=FALSE)
idx <- read.table("../ag.idx", sep=" ", header=FALSE)

p1negs <- vector()
p2negs <- vector()
p3negs <- vector()

fmp1 <- fm[,p1mat]
fmp2 <- fm[,p2mat]
fmp3 <- fm[,p3mat]

#this takes ~10-15min
for(i in 1:nrow(fm))
{
	if(i %% 10000 == 0)
	{
		print(i)
	}
	p1negs[i] <- length(fmp1[i,][fmp1[i,] == 0])
	p2negs[i] <- length(fmp2[i,][fmp2[i,] == 0])
	p3negs[i] <- length(fmp3[i,][fmp3[i,] == 0])
}

trimFM <- fm[p1negs < 1 & p2negs < 1 & p3negs < 1,]
trimSM <- sm[p1negs < 1 & p2negs < 1 & p3negs < 1,]

#64.25% of matrix kept 

xvals <- trimSM[,2]

#Aggregate normalisation per population
trimFC <- (trimFM[,3:47] - trimSM[,3:47])

trimFC[,p1] <- trimFC[,p1]/rowSums(trimFM[,p1mat])
trimFC[,p2] <- trimFC[,p2]/rowSums(trimFM[,p2mat])
trimFC[,p3] <- trimFC[,p3]/rowSums(trimFM[,p3mat])

#per individual normalisation
trimFC2 <- (trimFM[,3:47] - trimSM[,3:47])/trimFM[,3:47]

trimFC2[trimFC == -Inf] <- 0
trimFC2[trimFC == Inf] <- 0
trimFC2[is.na(trimFC)] <- 0

trimFC[trimFC == -Inf] <- 0
trimFC[trimFC == Inf] <- 0
trimFC[is.na(trimFC)] <- 0

stp1 <- apply(trimFC[,p1], 1, sd)
stp2 <- apply(trimFC[,p2], 1, sd)
stp3 <- apply(trimFC[,p3], 1, sd)

stp1_2 <- apply(trimFC2[,p1], 1, sd)
stp2_2 <- apply(trimFC2[,p2], 1, sd)
stp3_2 <- apply(trimFC2[,p3], 1, sd)

rmst1 <- rollmean(stp1, 20)
rmst2 <- rollmean(stp2, 20)
rmst3 <- rollmean(stp3, 20)

rmst1_2 <- rollmean(stp1_2, 20)
rmst2_2 <- rollmean(stp2_2, 20)
rmst3_2 <- rollmean(stp3_2, 20)

Mscfs <- vector()
Mlocs <- vector()
Mresults <- vector()

Fscfs <- vector()
Flocs <- vector()
Fresults <- vector()

refIndex1 <- 1:length(rmst1)

idx2 <- idx[idx[,2] > 100000,]

lvls <- idx2[,1]
cnt <- 1

medIntras <- vector()
flexIntras <- vector()

#LOOP BEGINS HERE
for( i in lvls)
{
	subs <- bins[bins[,1] == i,]
	ravg <- rollmean(subs[,2], 250, fill=99, na.pad = TRUE)
	ravg <- ravg * 2

	bigIndex <- refIndex1[trimFM[,1] == i]

	intraPop <- (rmst1[bigIndex] + rmst2[bigIndex] + rmst3[bigIndex])/3

	meds <- med_ids[med_ids[,1] == i,2]
	flexs <- flex_ids[flex_ids[,1] == i,2]
	mInds <- floor(meds/50);
	fInds <- floor(flexs/50);

	totRow <- length(ravg)
	totRow2 <- length(intraPop)

	mInds <- sapply(mInds, FUN=function(x){max(x,1)})
	fInds <- sapply(fInds, FUN=function(x){max(x,1)})
	mInds <- sapply(mInds, FUN=function(x){min(x,totRow)})
	fInds <- sapply(fInds, FUN=function(x){min(x,totRow)})

	mIndsIntra <- ceiling(meds/2000)
	fIndsIntra <- ceiling(flexs/2000)

	mIndsIntra <- sapply(mIndsIntra, FUN=function(x){max(x,1)})
	fIndsIntra <- sapply(fIndsIntra, FUN=function(x){max(x,1)})
	mIndsIntra <- sapply(mIndsIntra, FUN=function(x){min(x,totRow2)})
	fIndsIntra <- sapply(fIndsIntra, FUN=function(x){min(x,totRow2)})

	Mscfs <- c(Mscfs, rep(i, length(mInds)))
	Mlocs <- c(Mlocs, meds)
	Mresults <- c(Mresults, ravg[mInds])
	medIntras <- c(medIntras, intraPop[mIndsIntra])

	Fscfs <- c(Fscfs, rep(i, length(fInds)))
	Flocs <- c(Flocs, flexs)
	Fresults <- c(Fresults, ravg[fInds])
	flexIntras <- c(flexIntras, intraPop[fIndsIntra])

	print(i)
}

Fout <- data.frame(Fscfs, Flocs, Fresults, flexIntras)
Mout <- data.frame(Mscfs, Mlocs, Mresults, medIntras)


