library(zoo)

bins <- read.table("bins.txt", sep="\t", header=FALSE)
ravg <- rollmean(bins[,2], 500, fill=0.4, na.pad = TRUE) * 2

Laccu <- 0
Baccu <- 0

Low_Laccu <- 0
Low_Baccu <- 0

Inblock <- 0

cnt <- 1

Blocks <- list()
Vars <- list()

Low_Blocks <- list()
Low_Vars <- list()

for(k in seq(0.25, 3.5, 0.25))
{
	Blocks[[cnt]] <- vector()
	Vars[[cnt]] <- vector()
	
	Low_Blocks[[cnt]] <- vector()
	Low_Vars[[cnt]] <- vector()

	for(i in 1:length(ravg))
	{
		if(i %% 100000 == 0)
		{
			print(i)
		}
		if(ravg[i] > k)
		{
			if(Inblock == 0 && Low_Laccu > 50000)
			{
				Low_Blocks[[cnt]] <- c(Low_Blocks[[cnt]], Low_Laccu)
				Low_Vars[[cnt]] <- c(Low_Vars[[cnt]], Low_Baccu)
			}
			Laccu <- Laccu + 50
			Baccu <- Baccu + bins[i,2]

			Low_Laccu <- 0
			Low_Baccu <- 0
			Inblock <- 1
		}
		else
		{
			if(Inblock == 1 && Laccu > 50000)
			{
				Blocks[[cnt]] <- c(Blocks[[cnt]], Laccu)
				Vars[[cnt]] <- c(Vars[[cnt]], Baccu)
			}
			Low_Laccu <- Low_Laccu + 50
			Low_Baccu <- Low_Baccu + bins[i,2]

			Laccu <- 0
			Baccu <- 0
			Inblock <- 0
		}
	}
	cnt <- cnt + 1
	print(k)
}

percs <- (Vars[[2]] / Blocks[[2]]) * 100
percs2 <- (Low_Vars[[6]] / Low_Blocks[[6]]) * 100

hibloc <- Blocks[[2]]
lobloc <- Low_Blocks[[6]]
n_hi <- hibloc * -1
n_lo <- lobloc * -1

lim <- 1100000

increment <- 0.04
yseq <- seq(increment, 4, increment)

densi <- vector()
Low_densi <- vector()
cnt <- 1
for(i in yseq)
{
	b <- i - increment
	densi[cnt] <- sum(c(0,hibloc[percs > b & percs < i]))
	Low_densi[cnt] <- sum(c(0,lobloc[percs2 > b & percs2 < i]))
	cnt <- cnt +1
}

Low_densi[is.na(Low_densi)] <- 0
densi[is.na(densi)] <- 0

denom <- max(Low_densi) / 1e06
densi <- densi / denom
Low_densi <- Low_densi / denom

hdx <- c(densi, rep(0, length(densi)))
ldx <- c(Low_densi, rep(0, length(Low_densi)))
yvals <- c(yseq, rev(yseq))

png("Block_Tower.png", width=600, height=1200)
par(cex=1.4)
plot(percs,percs, xlim=c(-lim, lim), ylim=c(0,4), type="n", xlab="Size of Block (bp)                 Total Density % (of genome)", ylab=c("% Polymorphism in Block"), main="Linkage blocks by Polymorphism Rate")
grid(lwd=2, col="grey4")
segments(x0=n_hi, y0=percs, x1=rep(0, length(hibloc)), y1=percs, col=rgb(0.1,0.7,0.9,0.2), lwd=5)
segments(x0=n_lo, y0=percs2, x1=rep(0, length(lobloc)), y1=percs2, col=rgb(0.9,0.7,0.1,0.2), lwd=5)
polygon(x=hdx, y=yvals, col=rgb(0.1,0.7,0.9,0.4))
polygon(x=ldx, y=yvals, col=rgb(0.9,0.7,0.1,0.4))
abline(v=0, lwd=8, col="black")
legend("topright", lwd=c(8,8), col=c(rgb(0.1,0.7,0.9,0.6), rgb(0.9,0.7,0.1,0.6)), legend=c("Consistent > 0.5% Blocks", "Consistent < 1.5% Blocks"))
dev.off()

png("Block_Dist_1.png", width=1200, height=800)
par(cex=1.2)
scatter.smooth(percs,Blocks[[2]], col=rgb(0.3, 0.7, 0.9, 0.5), pch="@", main="20kb+ block regions > 1% polymophism (A.gracilis Whole Genome)", xlab="% Polymorphism In Block", ylab="Block Size (bp)")
legend("topright", lty=1, lwd=1, col="black", legend="LOWESS Trend")
grid()
dev.off()
