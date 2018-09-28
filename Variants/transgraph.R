vcfs <- read.table("vcf_bins.txt", sep="\t", header=TRUE)
vcSUm <- read.table("vcf_counts.txt", sep="\t", header=TRUE)
bins <- read.table("bins.txt", sep="\t", header=FALSE)

bins2 <- bins
bins2[,2] <- rollmean(bins[,2], 1000, fill=0.4, na.pad = TRUE)
bins2[,2] <- bins2[,2] * 2

png("simplePie.png", height=700, width=800)
par(cex=1.2)
pie(colSums(vcSUm[,3:6]), labels=c("CT/GA Transition 58%", "GT/CA Transversion 21%", "AT Transversion 13%", "CG Transversion 6%"), main="Transitions and Transversion Proportions in A.gracilis Genome")
dev.off()
vcfsRoll <- vcfs

lvls <- levels(vcfs[,1])

vcfsRoll[,3] <- rollmean(vcfs[,3], 1000, fill=0.4, na.pad = TRUE)
vcfsRoll[,4] <- rollmean(vcfs[,4], 1000, fill=0.4, na.pad = TRUE)
vcfsRoll[,5] <- rollmean(vcfs[,5], 1000, fill=0.4, na.pad = TRUE)
vcfsRoll[,2] <- rollmean(vcfs[,2], 1000, fill=0.4, na.pad = TRUE)

rsum <- rowSums(vcfsRoll[,2:5])
rsum[rsum == 0] <- 0.0001


png("PieDensity.png", width=1000, height=800)
par(cex=1.4)
plot(density((vcfsRoll[,2]/rsum)), type="l", col="firebrick", lwd=3, ylim=c(0,20), xlab="25kbp Mean Variant Proportion", main=" Total A.gracilis Genome SNP type Proportion")
points(density((vcfsRoll[,3]/rsum)), type="l", col="yellow2", lwd=3)
points(density((vcfsRoll[,4]/rsum)), type="l", col="green3", lwd=3)
points(density((vcfsRoll[,5]/rsum)), type="l", col="pink4", lwd=3)
legend("topright", lwd=c(3,3,3,3), col=c("firebrick", "yellow2", "green3", "pink4"), legend=c("CT/GA Transition 58%", "GT/CA Transversion 21%", "AT Transversion 13%", "CG Transversion 6%"))
dev.off()

vcfsRoll[,2] <- vcfsRoll[,2]/rsum
vcfsRoll[,3] <- vcfsRoll[,3]/rsum
vcfsRoll[,4] <- vcfsRoll[,4]/rsum
vcfsRoll[,5] <- vcfsRoll[,5]/rsum

bins

xlevs <- seq(50, 50 * nrow(subs), 50)
for( i in lvls)
{
	subs <- bins[bins[,1] == i,]

	#only 850k+ scaffolds
	if(nrow(subs) > 17000)
	{
		svcf <- vcfsRoll[vcfsRoll[,1] == i,]
		cnt <- 1

		xlevs <- seq(50, 50 * nrow(subs), 50)
		ravg <- rollmean(subs[,2], 500, fill=0.4, na.pad = TRUE)
		ravg <- ravg * 2
		subs[,2] <- subs[,2] * 2

		png(paste("scaffgraph4/plot for ", i, ".png", sep=""), width=1500, height=900)

		par(mfrow = c(2,1), cex=1.6, fig=c(0,1,0.4,1), mar=c(0,4,2,2))
		plot(xlevs, subs[,2], xlim=c(0, max(xlevs)), type="l", ylim=c(0.2,8), col=rgb(0.1,0.1,0.1,0.2), main=paste("Polymorphism Rates Across ", i, sep=""), xaxt='n', xlab="", ylab="% Polymorphisms")
		points(xlevs, ravg, type="l", col="orange", lwd=2)
		legend("topright", lwd=c(1,3,4,4,4,4), col=c(rgb(0.1,0.1,0.1), "orange", "firebrick", "yellow2", "green3", "pink4"), legend=c("50bp Mean", "25Kbp Rolling Mean", "CT/GA Transition 58%", "GT/CA Transversion 21%", "AT Transversion 13%", "CG Transversion 6%"))

		xlevs2 <- seq(50, 50 * length(svcf[,2]), 50)
		par(cex=1.6, fig=c(0,1,0,0.4), mar=c(3,4,0,2), new=TRUE)
		plot(xlevs2, svcf[,2], ylim=c(0,1), xlim=c(0, max(xlevs)), col="firebrick", lwd=4, xlab="DNA Bases in Scaffold", type="l", main="", ylab="SNP Ratio")
		points(xlevs2, svcf[,3], type="l", col="yellow2", lwd=4)
		points(xlevs2, svcf[,4], type="l", col="green3", lwd=4)
		points(xlevs2, svcf[,5], type="l", col="pink4", lwd=4)
		grid(col="black")
		dev.off()
		}
	}
}

vvv <- vector()
ct <- vector()
gt <- vector()
at <- vector()
cg <- vector()

for( i in lvls)
{
	subs <- bins[bins[,1] == i,]
	ravg <- rollmean(subs[,2], 1000, fill=0.4, na.pad = TRUE)
	svcf <- vcfsRoll[vcfsRoll[,1] == i, 2:5]
	
	toplev <- min(length(vcfsRoll[,2]), length(subs[,2]))
	
	vvv <- c(vvv, ravg[1:toplev])
	ct <- c(ct, svcf[,1][1:toplev])
	gt <- c(gt, svcf[,2][1:toplev])
	at <- c(at, svcf[,3][1:toplev])
	cg <- c(cg, svcf[,4][1:toplev])
}

x <- data.frame(vars = vvv, ct, gt, at, cg)
x[is.na(x)] <- 0

x2 <- x[sample(1:nrow(x), 100000,replace=FALSE),]

png("CTloess.png", height=700, width=700)
par(cex=1.2)
scatter.smooth(x2$vars*2, x2$ct, lpars = list(col = "firebrick", lwd = 3, lty = 3), xlim=c(0,5), ylim=c(0,1), main="LOESS, Polymorphism Rate vs Proportion of CT/GA Transitions", xlab="% Polymorphims", ylab="Variant Ratio")
abline(h=mean(x2$ct), col="green", lty=2)
legend("topright", lwd=c(3,1), col=c("firebrick", "green"), lty=c(3,2), legend=c("LOESS", "Mean"))
dev.off()

png("Allloess.png", height=1400, width=1400)
par(mfrow=c(2,2), cex=1.8)

scatter.smooth(x2$vars*2, x2$ct, lpars = list(col = "firebrick", lwd = 3, lty = 3), xlim=c(0,5), ylim=c(0,1), main="CT/GA Transitions", xlab="% Polymorphims", ylab="Variant Ratio")
abline(h=mean(x2$ct), col="green", lty=2)
legend("topright", lwd=c(3,1), col=c("firebrick", "green"), lty=c(3,2), legend=c("LOESS", "Mean"))

scatter.smooth(x2$vars*2, x2$gt, lpars = list(col = "yellow2", lwd = 3, lty = 3), xlim=c(0,5), ylim=c(0,1), main="GT/CA Transversions", xlab="% Polymorphims", ylab="Variant Ratio")
abline(h=mean(x2$gt), col="green", lty=2)
legend("topright", lwd=c(3,1), col=c("yellow2", "green"), lty=c(3,2), legend=c("LOESS", "Mean"))

scatter.smooth(x2$vars*2, x2$at, lpars = list(col = "green3", lwd = 3, lty = 3), xlim=c(0,5), ylim=c(0,1), main="AT Transversions", xlab="% Polymorphims", ylab="Variant Ratio")
abline(h=mean(x2$at), col="green", lty=2)
legend("topright", lwd=c(3,1), col=c("green3", "green"), lty=c(3,2), legend=c("LOESS", "Mean"))

scatter.smooth(x2$vars*2, x2$cg, lpars = list(col = "pink4", lwd = 3, lty = 3), xlim=c(0,5), ylim=c(0,1), main="CG Transversions", xlab="% Polymorphims", ylab="Variant Ratio")
abline(h=mean(x2$cg), col="green", lty=2)
legend("topright", lwd=c(3,1), col=c("pink4", "green"), lty=c(3,2), legend=c("LOESS", "Mean"))

dev.off()






