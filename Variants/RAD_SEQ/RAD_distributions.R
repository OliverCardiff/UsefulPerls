library(MASS)
library(RColorBrewer)
library(plot3D)

rf <- colorRampPalette(rev(brewer.pal(11,'PuOr')))
r <- rf(32)

png("Pop_diffs.png", width=1200, height=1000)
par(cex=1.5)
plot(density(rmst3), type="l", col="darkkhaki", lwd=4, main="Whole Genome Intra-population Sequence Variation Estimation", ylab="Density", xlab="Population depth-Normalised pileup SD in stringency difference", xlim=c(0,0.06))
points(density(rmst1), type="l", col="darkolivegreen", lwd=4)
points(density(rmst2), type="l", col="green2", lwd=4)
legend("topright", col=c("darkolivegreen", "green2", "darkkhaki"), lwd=c(4,4,4), legend=c("Pop M", "Pop F", "Pop S"))
dev.off()

similarity <- (abs(rmst1-rmst2) + abs(rmst2-rmst3) + abs(rmst1-rmst3))/3
levelPop <- (rmst1 + rmst2 + rmst3)/3
similarity <- similarity/levelPop

kSimLev <- kde2d(similarity, levelPop)

allvecs1 <- list()
allvecs2 <- list()
allvecs3 <- list()

lvls <- levels(fm[,1])

idx[,4] <- (idx[,3] * 150) / idx[,2]
idxC <- idx[idx[,4] < 200 & idx[,2] > 500000,]
lvls <- idxC[,1]

for(i in 1:46)
{
	allvecs1[[i]] <- vector()
	allvecs2[[i]] <- vector()
	allvecs3[[i]] <- vector()
}

#LOOP BEGINS HERE
for( i in lvls)
{
	subs <- bins[bins[,1] == i,]

	#only 500k+ scaffolds
	if(nrow(subs) > 10000)
	{
		bigIndex <- refIndex1[trimFM[,1] == i]
		divIndex <- refIndex2[fm[,1] == i]
		coords <- trimFM[bigIndex,2]
		rInds <- coords/50
		ravg <- (rollmean(subs[,2], 500, fill=0.4, na.pad = TRUE)) * 2
	
		AvgPoints <- ravg[rInds]
		p1Points <- rmst1_2[bigIndex]
		p2Points <- rmst2_2[bigIndex]
		p3Points <- rmst3_2[bigIndex]
	
		cnt2 <- 1

		for(b in seq(0.5, 5, 0.1))
		{
			b2 <- b - 0.5
	
			vM1 <- p1Points[AvgPoints < b & AvgPoints > b2]
			vM2 <- p2Points[AvgPoints < b & AvgPoints > b2]
			vM3 <- p3Points[AvgPoints < b & AvgPoints > b2]
			
			if(length(vM1) > 0)
			{
				allvecs1[[cnt2]] <- c(allvecs1[[cnt2]], vM1)
			}
			if(length(vM2) > 0)
			{
				allvecs2[[cnt2]] <- c(allvecs2[[cnt2]], vM2)
			}
			if(length(vM3) > 0)
			{
				allvecs3[[cnt2]] <- c(allvecs3[[cnt2]], vM3)
			}
		
			cnt2 <- cnt2 + 1
		}
	}
	print(i)
}

tot <- cnt2 - 1;

for(i in 1:46)
{
	allvecs1[[i]][is.na(allvecs1[[i]])] <- 0.017
	allvecs2[[i]][is.na(allvecs2[[i]])] <- 0.017
	allvecs3[[i]][is.na(allvecs3[[i]])] <- 0.017
	
	#allvecs1[[i]] <- (allvecs1[[i]] / sum(allvecs1[[i]])) * 100
	#allvecs2[[i]] <- (allvecs2[[i]] / sum(allvecs2[[i]])) * 100
	#allvecs3[[i]] <- (allvecs3[[i]] / sum(allvecs3[[i]])) * 100
}

mat <- rainbow(30)

plot(density(allvecs1[[1]]), type="l", col=mat[1], xlim=c(0,0.2), ylim=c(0,14))
for(i in 2:40)
{
	points(density(allvecs1[[i]]), type="l", col=mat[i])
}

allDens1 <- matrix(0,46,30)
allDens2 <- matrix(0,46,30)
allDens3 <- matrix(0,46,30)
for(i in 1:46)
{
	allDens1[i,] <- hist(allvecs1[[i]], breaks=c(seq(0,0.5,0.01),Inf), plot=FALSE)$density[1:30]
	allDens2[i,] <- hist(allvecs2[[i]], breaks=c(seq(0,0.5,0.01),Inf), plot=FALSE)$density[1:30]
	allDens3[i,] <- hist(allvecs3[[i]], breaks=c(seq(0,0.5,0.01),Inf), plot=FALSE)$density[1:30]
}

plot(seq(0,0.29,0.01), allDens1[1,], type="l", col=mat[1], ylim=c(0,13))
for(i in 1:40)
{
	points(seq(0,0.29,0.01), allDens1[i,], type="l", col=mat[i])
}

for(i in 1:360)
{
png(paste("img", i, ".png", sep=""), width=600, height=600)
persp3D(z=t(allDens1), col=ramp.col(c("green4", "gold", "purple")), phi=35, theta=i, shade=0.2, ylab="25Kbp Polymorphism 0-4%", xlab="Intra-pop. Seq Variation 0-5%", zlab="Density", main="A.gracilis Ref. Genome Polymorphisms vs Intra-population Variation (M)")
dev.off()
png(paste("F-plots/img", i, ".png", sep=""), width=600, height=600)
persp3D(z=t(allDens2), col=ramp.col(c("green4", "gold", "purple")), phi=35, theta=i, shade=0.2, ylab="25Kbp Polymorphism 0-4%", xlab="Intra-pop. Seq Variation 0-5%", zlab="Density", main="A.gracilis Ref. Genome Polymorphisms vs Intra-population Variation (F)")
dev.off()
png(paste("S-plots/img", i, ".png", sep=""), width=600, height=600)
persp3D(z=t(allDens3), col=ramp.col(c("green4", "gold", "purple")), phi=35, theta=i, shade=0.2, ylab="25Kbp Polymorphism 0-4%", xlab="Intra-pop. Seq Variation 0-5%", zlab="Density", main="A.gracilis Ref. Genome Polymorphisms vs Intra-population Variation (S)")
dev.off()
}
