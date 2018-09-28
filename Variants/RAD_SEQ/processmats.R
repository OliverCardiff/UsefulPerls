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
#trimFC[,1] <- scs[p1negs < 1 & p2negs < 1 & p3negs < 1,1]
#trimFC[,2] <- loci[p1negs < 1 & p2negs < 1 & p3negs < 1,1]

disto <- rowSums(fm[,3:47])
indivs <- colSums(fm[,3:47])
indivs2 <- colSums(sm[,3:47])

png("CutoffQC.png", width=900, height=700)
plot(1:length(indivs), indivs, type="h", lwd=5, main="Read alignment Counts Per Individual", xlab="Individual Worms", ylab="Read Alignments")
abline(h=500000, lty=2, col="red", lwd=4)
legend("topright", lty=c(1,2), lwd=c(5,4), col=c("red", "black"), legend=c("QC Cutoff", "Levels"))
dev.off()

disto1 <- rowSums(fm[,p1mat])
disto2 <- rowSums(fm[,p2mat])
disto3 <- rowSums(fm[,p3mat])

fdisto1 <- rowSums(sm[,p1mat])
fdisto2 <- rowSums(sm[,p2mat])
fdisto3 <- rowSums(sm[,p3mat])

png("Reference.png", width=1020, height=2200)
par(mfrow=c(3,1), cex=1.5)
plot(250:500, disto1[250:500], type="l")
points(250:500, fdisto1[250:500], type="l", lty=2, lwd=2, col="black", main="Population - M", ylab="Read Depth", xlab="2Kbp Bins")
plot(250:500, disto2[250:500], type="l")
points(250:500, fdisto2[250:500], type="l", lty=2, lwd=2, col="black", main="Population - F", ylab="Read Depth", xlab="2Kbp Bins")
plot(250:500, disto3[250:500], type="l")
points(250:500, fdisto3[250:500], type="l", lty=2, lwd=2, col="black", main="Population - S", ylab="Read Depth", xlab="2Kbp Bins")
dev.off()

comp1 <- fdisto1 - disto1
comp2 <- fdisto2 - disto2
comp3 <- fdisto3 - disto3

comp1 <- comp1/mean(comp1)
comp2 <- comp2/mean(comp2)
comp3 <- comp3/mean(comp3)

comp12 <- (comp1-comp2) * (1/disto) * 300
comp23 <- (comp2-comp3) * (1/disto) * 300
comp13 <- (comp1-comp3) * (1/disto) * 300

diverge <- (abs(comp12) + abs(comp23) + abs(comp13)) / 3
diverge[diverge == Inf] <- 0
diverge[is.na(diverge)] <- 0

diverge_roll <- rollmean(diverge, 20)

divXvals <- fm[,2]

png("Reference2a.png", width=1020, height=2200)
par(mfrow=c(3,1), cex=1.5)
plot(250:500, comp1[250:500], type="l", main="F, M & S Populations read alignment differences", col="darkorange2", xlab="2Kbp Positions", ylab="Coverage Differences")
points(250:500, comp2[250:500], type="l", col="darkolivegreen3")
points(250:500, comp3[250:500], type="l", col="cornflowerblue")
legend("topright", lty=c(1,1,1), lwd=c(1,1,1), col=c("darkorange2", "darkolivegreen3", "cornflowerblue"), legend=c("Pop. M", "Pop. F", "Pop. S"))
plot(250:500, comp12[250:500], type="l", main="Pairwise Population read Alignment difference comparisons", xlab="2Kbp Positions", ylab="Read Depth % difference of differences")
points(250:500, comp23[250:500], type="l", col="red")
points(250:500, comp13[250:500], type="l", col="blue")
legend("topright", lty=c(1,1,1), lwd=c(1,1,1), col=c("black", "red", "blue"), legend=c("M vs F", "F vs S", "M vs S"))
plot(250:500, diverge[250:500], type="l", main="Normalised Population separation by read depth differnece", xlab="2Kbp Positions", ylab="Pairwise alignment depth divergence %, Normalised by absolute levels")
dev.off()

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

plot(xvals[prange], stp1[prange], type="l", col=rgb(0.1,0.1,0.1,0.3), lty=2)
points(xvals[prange], stp2[prange], type="l", col=rgb(1,0.1,0.1,0.3), lty=2)
points(xvals[prange], stp3[prange], type="l", col=rgb(0.1,0.1,1,0.3), lty=2)
points(xvals[prange], rmst1[prange], type="l", col="black", lwd=2)
points(xvals[prange], rmst2[prange], type="l", col="red", lwd=2)
points(xvals[prange], rmst3[prange], type="l", col="blue", lwd=2)

rmst1 <- rollmean(stp1, 10)
rmst2 <- rollmean(stp2, 10)
rmst3 <- rollmean(stp3, 10)

rmst1_2 <- rollmean(stp1_2, 10)
rmst2_2 <- rollmean(stp2_2, 10)
rmst3_2 <- rollmean(stp3_2, 10)

bins <- read.table("../bins.txt", sep="\t", header=FALSE)
ravg <- rollmean(bins[,2], 500, fill=0.4, na.pad = TRUE)
idx <- read.table("../ag.idx", sep=" ", header=FALSE)
covs <- read.table("../wigmat.txt", sep="\t", header=FALSE)
orfs <- read.table("../orfs.txt", sep="\t", header=FALSE)
forfs <- orfs[orfs[,2] > 350,]


idx[,4] <- (idx[,3] * 150) / idx[,2]
idxC <- idx[idx[,4] < 200 & idx[,2] > 500000,]
lvls <- idxC[,1]

refIndex1 <- 1:length(rmst1)
refIndex2 <- 1:nrow(fm)

#LOOP BEGINS HERE
for( i in lvls)
{
subs <- bins[bins[,1] == i,]

#only 850k+ scaffolds
if(nrow(subs) > 17000)
{

bigIndex <- refIndex1[trimFM[,1] == i]
divIndex <- refIndex2[fm[,1] == i]

cnt <- 1

scov <- covs[covs[,1] == i,]
dlevs <- seq(10, 10* nrow(scov), 10)

sorf <- forfs[forfs[,1] == i,]

xlevs <- seq(50, 50 * nrow(subs), 50)
ravg <- rollmean(subs[,2], 500, fill=99, na.pad = TRUE)
ravg <- ravg * 2
subs[,2] <- subs[,2] * 2

png(paste("scaffgraph/plot for ", i, ".png", sep=""), width=2000, height=1350)

#Intra-population variation
par(mfrow = c(4,1), cex=2, fig=c(0,1,0.7,1), mar=c(0,4,2,2))
plot(xvals[bigIndex], stp1[bigIndex], type="l", col=rgb(0.1,0.1,0.1,0.3), ylim=c(0,0.08), lty=2, xaxt="n", main=paste("Polymorphism Rates Across ", i, " Read Pileup, Plus RADtag stacks", sep=""), ylab="Intra-Pop. Divergence", xlim=c(0, max(xlevs)))
points(xvals[bigIndex], stp2[bigIndex], type="l", col=rgb(1,0.1,0.1,0.3), lty=2)
points(xvals[bigIndex], stp3[bigIndex], type="l", col=rgb(0.1,0.1,1,0.3), lty=2)
points(xvals[bigIndex], rmst1[bigIndex], type="l", col="black", lwd=2)
points(xvals[bigIndex], rmst2[bigIndex], type="l", col="red", lwd=2)
points(xvals[bigIndex], rmst3[bigIndex], type="l", col="blue", lwd=2)
abline(h=0.04, lty=2)
legend("topright", lwd=c(2,2,2), col=c("black", "red", "blue"), legend=c("Population M", "Population F", "Population S"))

##Divergence between populations
par(cex=2, fig=c(0,1,0.5,0.7), mar=c(0,4,0,2), new=TRUE)
plot(divXvals[divIndex], diverge_roll[divIndex], lwd=2, type="l", xaxt="n", ylab="Inter-Pop. DNA diverge.", xlim=c(0, max(xlevs)))

par(cex=2, fig=c(0,1,0.2,0.5), mar=c(0,4,0,2), new=TRUE)
plot(xlevs, subs[,2], xlim=c(0, max(xlevs)), type="l", ylim=c(0.2,7), col=rgb(0.1,0.1,0.1,0.2), xaxt='n', xlab="", ylab="% Polymorphisms")
points(xlevs, ravg, type="l", col="orange", lwd=2)
symbols(sorf[,3], sorf[,4], circles=sorf[,2], inches=0.1, add=TRUE, fg="cyan", bg=rgb(0.1,0.5,1,0.3))
legend("topright", lwd=c(1,3, 2), lty=c(1,1,NA), pch=c(NA,NA,21), col=c(rgb(0.1,0.1,0.1), "orange", rgb(0.1,0.5,1,0.3)), legend=c("50bp Mean", "25Kbp Rolling Mean", "ORFS"))

par(cex=2, fig=c(0,1,0,0.2), mar=c(3,4,0,2), new=TRUE)
plot(dlevs, scov[,2], type="h", ylim=c(0,160), xlim=c(0, max(dlevs)), col="green4", xlab="DNA Bases in Scaffold", main="", ylab="Read pileup")
abline(h=85, lty=2)
abline(h=42.5, lty=2)

dev.off()
}
}

