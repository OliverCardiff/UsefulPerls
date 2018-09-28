library(zoo)
library(reshape2)
setwd("C:/cygwin64/home/Oliver Rimington/Documents/Variants/")
bins <- read.table("bins.txt", sep="\t", header=FALSE)
covs <- read.table("wigmat.txt", sep="\t", header=FALSE)
patches <- read.table("mats2.txt", sep="\t", header=FALSE)

orfs <- read.table("orfs.txt", sep="\t", header=FALSE)
forfs <- orfs[orfs[,2] > 350,]

idx <- read.table("ag.idx", sep="\t", header=FALSE)

idx[,4] <- (idx[,3] * 150) / idx[,2]

idxC <- idx[idx[,4] < 200 & idx[,2] > 500000,]

lvls <- idxC[,1]

#Load fullMat
fmat <- read.table("fullmatrix.txt", sep="\t", header=TRUE)
filmat <- fmat[!is.na(fmat[,48]),]

#LOOP BEGINS HERE
for( i in lvls)
{
subs <- bins[bins[,1] == i,]

#only 850k+ scaffolds
if(nrow(subs) > 17000)
{

cnt <- 1

sfil <- filmat[filmat[,1] == i,]

###########################################################

p1 <- c(1,4,5,6,7,8,15,16,20,21,22,24,27,28,29)
p1 <- p1 + 2
p2 <- c(2,3,9,10,11,12,13,14,17,18,19,23,25,26)
p2 <- p2 + 2
p3 <- 30:45
p3 <- p3 + 2

fp1 <- sfil[,p1]
fp2 <- sfil[,p2]
fp3 <- sfil[,p3]

cs1 <- rowSums(fp1)
cs2 <- rowSums(fp2)
cs3 <- rowSums(fp3)

fullOc <- cs1 + cs2 + cs3

sub12 <- abs(cs1 - cs2)
sub13 <- abs(cs1 - cs3)
sub23 <- abs(cs2 - cs3)

Invar <- sub12 + sub13 + sub23

Invar <- Invar[order(sfil[,2])]
fullOc <- fullOc[order(sfil[,2])]
sfil <- sfil[order(sfil[,2]),]

Invar <- rollmean(Invar, 12, fill=0.4, na.pad = TRUE)
fullOc <- rollmean(fullOc, 12, fill=0.4, na.pad = TRUE)

###########################################################

fp1 <- cbind(locus=sfil[,2], fp1)
fp2 <- cbind(locus=sfil[,2], fp2)
fp3 <- cbind(locus=sfil[,2], fp3)

mfp1 <- melt(fp1, id.var="locus")
mfp2 <- melt(fp2, id.var="locus")
mfp3 <- melt(fp3, id.var="locus")

yv1 <- vector()
yv2 <- vector()
yv3 <- vector()

for(j in 1:15)
{
	yv1 <- c(yv1, rep(j,nrow(fp1)))
}
for(j in 16:29)
{
	yv2 <- c(yv2, rep(j,nrow(fp2)))
}
for(j in 30:45)
{
	yv3 <- c(yv3, rep(j,nrow(fp3)))
}

###########################################################
spatch <- patches[patches[,1] == i,]

sp1 <- spatch[spatch[,2] < 16,]
sp2 <- spatch[spatch[,2] > 15 & spatch[,2] < 31,]
sp3 <- spatch[spatch[,2] > 30,]

scov <- covs[covs[,1] == i,]
dlevs <- seq(10, 10* nrow(scov), 10)

sorf <- forfs[forfs[,1] == i,]

xlevs <- seq(50, 50 * nrow(subs), 50)
ravg <- rollmean(subs[,2], 500, fill=0.4, na.pad = TRUE)
ravg <- ravg * 2
subs[,2] <- subs[,2] * 2

png(paste("scaffgraph7/plot for ", i, ".png", sep=""), width=2000, height=1350)

par(mfrow = c(4,1), cex=2, fig=c(0,1,0.7,1), mar=c(0,4,2,2))
symbols(mfp1[,1], yv1, squares=mfp1[,3], inches=0.08, ylim=c(0,45), xlim=c(0, max(dlevs)), fg=rgb(0.54,0.13,0.13,0.5), bg=rgb(0.54,0.13,0.13,0.5),xaxt='n', main=paste("Polymorphism Rates Across ", i, " Read Pileup, Plus RADtag stacks", sep=""), ylab="Individuals")
symbols(mfp2[,1], yv2, squares=mfp2[,3], inches=0.08, fg=rgb(0.43,0.545,0.23,0.5), bg=rgb(0.43,0.545,0.23,0.5), add=TRUE)
symbols(mfp3[,1], yv3, squares=mfp3[,3], inches=0.08, fg=rgb(0.8, 0.35, 0.18, 0.5), bg=rgb(0.8, 0.35, 0.18, 0.5), add=TRUE)
abline(h=15.5, lwd=2, lty=2)
abline(h=29.5, lwd=2, lty=2)
#legend("topright", pch=c(22,22,22),lwd=c(3,3,3), col=c("azure3", "darkolivegreen4", "coral3"), legend=c("Population 1", "Population 2", "Population 3"))

##Stacks figures
par(cex=2, fig=c(0,1,0.5,0.7), mar=c(0,4,0,2), new=TRUE)
plot(sfil[,2], fullOc, xlim=c(0, max(xlevs)), type="l", lwd=2, col="red", xaxt='n', xlab="", ylab="Stack Count")
points(sfil[,2], Invar, type="l", col="blue", lwd=2)
legend("topright", lwd=c(2,2), lty=c(1,1), col=c("red", "blue"), legend=c("Stack Occupancy", "Inter-pop Variation"))

par(cex=2, fig=c(0,1,0.2,0.5), mar=c(0,4,0,2), new=TRUE)
plot(xlevs, subs[,2], xlim=c(0, max(xlevs)), type="l", ylim=c(0.2,7), col=rgb(0.1,0.1,0.1,0.2), xaxt='n', xlab="", ylab="% Polymorphisms")
points(xlevs, ravg, type="l", col="orange", lwd=2)
symbols(sorf[,3], sorf[,4], circles=sorf[,2], inches=0.1, add=TRUE, fg="cyan", bg=rgb(0.1,0.5,1,0.3))
legend("topright", lwd=c(1,3, 2), lty=c(1,1,NA), pch=c(NA,NA,21), col=c(rgb(0.1,0.1,0.1), "orange", rgb(0.1,0.5,1,0.3)), legend=c("50bp Mean", "25Kbp Rolling Mean", "ORFS"))

par(cex=1.5, fig=c(0,1,0,0.2), mar=c(3,4,0,2), new=TRUE)
plot(dlevs, scov[,2], type="h", ylim=c(0,160), xlim=c(0, max(dlevs)), col="green4", xlab="DNA Bases in Scaffold", main="", ylab="Read pileup")
abline(h=85, lty=2)
abline(h=42.5, lty=2)

dev.off()
}
}
