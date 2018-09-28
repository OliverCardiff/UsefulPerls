library(zoo)
setwd("C:/cygwin64/home/Oliver Rimington/Documents/Variants/")
bins <- read.table("bins.txt", sep="\t", header=FALSE)
covs <- read.table("wigmat_fix.txt", sep="\t", header=FALSE)
patches <- read.table("mats2.txt", sep="\t", header=FALSE)

orfs <- read.table("orfs.txt", sep="\t", header=FALSE)
forfs <- orfs[orfs[,2] > 350,]

idx <- read.table("ag.idx", sep=" ", header=FALSE)

idx[,4] <- (idx[,3] * 150) / idx[,2]

idxC <- idx[idx[,4] < 200 & idx[,2] > 500000,]

lvls <- idxC[,1]

subs <- bins[bins[,1] == lvls[1],]

xlevs <- seq(50, 50 * nrow(subs), 50)
ravg <- rollmean(subs[,2], 1000, fill=0.4, na.pad = TRUE)
ravg <- ravg * 2

plot(xlevs, subs[,2], type="l", ylim=c(0,8), col=rgb(0.1,0.1,0.1,0.2))
points(xlevs, ravg, type="l", col="orange", lwd=2)

for( i in lvls)
{
subs <- bins[bins[,1] == i,]

#only 850k+ scaffolds
if(nrow(subs) > 17000)
{

cnt <- 1

scov <- covs[covs[,1] == i,]
dlevs <- seq(10, 10* nrow(scov), 10)

sorf <- forfs[forfs[,1] == i,]

xlevs <- seq(50, 50 * nrow(subs), 50)
ravg <- rollmean(subs[,2], 500, fill=0.4, na.pad = TRUE)
ravg <- ravg * 2
subs[,2] <- subs[,2] * 2

png(paste("scaffgraph_fin/plot for ", i, ".png", sep=""), width=900, height=600)

par(mfrow = c(2,1), cex=1.4, fig=c(0,1,0.3,1), mar=c(0,4,2,2))

plot(xlevs, subs[,2], xlim=c(0, max(xlevs)), type="l", ylim=c(0.2,6), col=rgb(0.1,0.1,0.1,0.2), main=paste("Polymorphism Rates Across ", i, sep=""), xaxt='n', xlab="", ylab="% Polymorphisms")
points(xlevs, ravg, type="l", col="orange", lwd=2)
symbols(sorf[,3], sorf[,4], circles=sorf[,2], inches=0.1, add=TRUE, fg="cyan", bg=rgb(0.1,0.5,1,0.3))
legend("topright", lwd=c(1,3, 2), lty=c(1,1,NA), pch=c(NA,NA,21), col=c(rgb(0.1,0.1,0.1), "orange", rgb(0.1,0.5,1,0.3)), legend=c("50bp Mean", "25Kbp Rolling Mean", "ORFS"))

par(cex=1.4, fig=c(0,1,0,0.3), mar=c(3,4,0,2), new=TRUE)
plot(dlevs, scov[,2], type="h", ylim=c(0,140), xlim=c(0, max(dlevs)), col=rgb(0.05,0.7,0.2,0.1), xlab="DNA Bases in Scaffold", main="", ylab="Read pileup")
abline(h=85, lty=2)
abline(h=42.5, lty=2)

dev.off()
}
}
