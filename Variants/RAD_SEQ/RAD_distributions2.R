
Iallvecs1 <- list()
Iallvecs2 <- list()
Iallvecs3 <- list()

lvls <- levels(fm[,1])

for(i in 1:46)
{
	Iallvecs1[[i]] <- vector()
	Iallvecs2[[i]] <- vector()
	Iallvecs3[[i]] <- vector()
}

refIndex1 <- 1:length(rmst1)
refIndex2 <- 1:nrow(fm)

for( i in lvls)
{
	subs <- bins[bins[,1] == i,]

	#only 500k+ scaffolds
	if(nrow(subs) > 10000)
	{
		bigIndex <- refIndex1[trimFM[,1] == i]
		divIndex <- refIndex2[fm[,1] == i]

		p1Points <- rmst1[bigIndex]
		p2Points <- rmst2[bigIndex]
		p3Points <- rmst3[bigIndex]

		#Converting poly% bins locations to intra-pop indexes
		ravg <- (rollmean(subs[,2], 500, fill=0.4, na.pad = TRUE)) * 2
		polyCoords <- (1:length(ravg)) * 50
		mInds <- floor(polyCoords/2000);
		mInds <- sapply(mInds, FUN=function(x){max(x,1)})
		totRow <- length(p1Points)
		mInds <- sapply(mInds, FUN=function(x){min(x,totRow)})

		p1refs <- p1Points[mInds]
		p2refs <- p2Points[mInds]
		p3refs <- p3Points[mInds]
	
		cnt2 <- 1

		for(b in seq(0.005, 0.05, 0.001))
		{
			b2 <- b - 0.005
	
			vM1 <- ravg[p1refs < b & p1refs > b2]
			vM2 <- ravg[p2refs < b & p2refs > b2]
			vM3 <- ravg[p3refs < b & p3refs > b2]
			
			if(length(vM1) > 0)
			{
				Iallvecs1[[cnt2]] <- c(Iallvecs1[[cnt2]], vM1)
			}
			if(length(vM2) > 0)
			{
				Iallvecs2[[cnt2]] <- c(Iallvecs2[[cnt2]], vM2)
			}
			if(length(vM3) > 0)
			{
				Iallvecs3[[cnt2]] <- c(Iallvecs3[[cnt2]], vM3)
			}
		
			cnt2 <- cnt2 + 1
		}
	}
	print(i)
}

tot <- cnt2 - 1;

plot(density(Iallvecs3[[1]]), type="l", col=mat[1], xlim=c(0,4))
for(i in 2:30)
{
	points(density(Iallvecs3[[i]]), type="l", col=mat[i])
}

for(i in 1:46)
{
	Iallvecs1[[i]][is.na(allvecs1[[i]])] <- 0.0017
	Iallvecs2[[i]][is.na(allvecs2[[i]])] <- 0.0017
	Iallvecs3[[i]][is.na(allvecs3[[i]])] <- 0.0017
	
	#allvecs1[[i]] <- (allvecs1[[i]] / sum(allvecs1[[i]])) * 100
	#allvecs2[[i]] <- (allvecs2[[i]] / sum(allvecs2[[i]])) * 100
	#allvecs3[[i]] <- (allvecs3[[i]] / sum(allvecs3[[i]])) * 100
}

IallDens1 <- matrix(0,30,50)
IallDens2 <- matrix(0,30,50)
IallDens3 <- matrix(0,30,50)
for(i in 1:30)
{
	IallDens1[i,] <- hist(Iallvecs1[[i]], breaks=c(seq(0,5,0.1),Inf), plot=FALSE)$density[1:50]
	IallDens2[i,] <- hist(Iallvecs2[[i]], breaks=c(seq(0,5,0.1),Inf), plot=FALSE)$density[1:50]
	IallDens3[i,] <- hist(Iallvecs3[[i]], breaks=c(seq(0,5,0.1),Inf), plot=FALSE)$density[1:50]
}

