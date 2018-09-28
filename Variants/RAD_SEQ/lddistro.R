library(zoo)
library(MASS)
library(RColorBrewer)

bins <- read.table("../bins.txt", sep="\t", header=FALSE)
ravg <- (rollmean(bins[,2], 500, fill=0.4, na.pad = TRUE)) * 2
idx <- read.table("../ag.idx", sep=" ", header=FALSE)
lds <- read.table("fsts/trim3.ld", sep="\t", header=FALSE)

idx[,4] <- (idx[,3] * 150) / idx[,2]
idxC <- idx[idx[,4] < 200 & idx[,2] > 100000,]
lvls <- idxC[,1]

lev2 <- levels(lds[,1])

sts <- vector()
eds <- vector()
p1 <- vector()
p2 <- vector()
r2 <- vector()


#LOOP BEGINS HERE
for( i in lev2)
{
	subs <- ravg[bins[,1] == i]
	sld <- lds[lds[,1] == i,]
	
	#only 100k+ scaffolds
	if(length(subs) > 2000 & nrow(sld) > 1)
	{
		poly1 <- floor(sld$V2 / 50)
		poly2 <- floor(sld$V3 / 50)
		lim <- min(length(ravg[poly1]),length(ravg[poly2]))
		
		p1 <- c(p1, ravg[poly1][1:lim])
		p2 <- c(p2, ravg[poly2][1:lim])
		
		sts <- c(sts, sld$V2[1:lim])
		eds <- c(eds, sld$V3[1:lim])
		r2 <- c(r2, sld$V5[1:lim])
		
		print(i)
	}
}

allmat <- rbind(sts, eds, p1, p2, r2)
allmat <- t(allmat)
lowmat <- allmat[allmat[,3] < 1.5 & allmat[,4] < 1.5,]
highmat <- allmat[allmat[,3] > 2 & allmat[,4] > 2,]



lowLen <- lowmat[,2] - lowmat[,1]
highLen <- highmat[,2] - highmat[,1]

twoLen <- twomat[,2] - twomat[,1]


for(i in seq(10000,50000,10000))
{
	lowb <- i - 10000
	
	aLen <- allmat[,2] - allmat[,1]
	twomat <- allmat[allmat[,3] > 2 & allmat[,4] > 2 & aLen < i & aLen > lowb,]
	
	td <- density(twomat[,5])
	if(i == 10000)
	{
		plot(td, type="l", lwd=4, col=rgb(0.1,0.1,0.1,0.5), ylim=c(0,7))
	} else if(i==50000)
	{
		points(td, type="l", lwd=4, col=rgb(0.1,0.1,0.1,0.5))
	}
	else{
		points(td, type="l", lwd=2, col=rgb(0.1,0.1,0.1,0.25))
	}
}

