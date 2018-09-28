library(zoo)

bins <- read.table("bins.txt", sep="\t", header=FALSE)
ravg <- rollmean(bins[,2], 500, fill=0.4, na.pad = TRUE)

#Load fullMat
fmat <- read.table("fullmatrix.txt", sep="\t", header=TRUE)
filmat <- fmat[!is.na(fmat[,48]),]


###########################################################


fp1 <- filmat[,3:17]
fp2 <- filmat[,18:32]
fp3 <- filmat[,33:47]

cs1 <- rowSums(fp1)
cs2 <- rowSums(fp2)
cs3 <- rowSums(fp3)

fullOc <- cs1 + cs2 + cs3

sub12 <- abs(cs1 - cs2)
sub13 <- abs(cs1 - cs3)
sub23 <- abs(cs2 - cs3)

Invar <- sub12 + sub13 + sub23

###########################################################

rlist <- list()
cnt <- 1
for( i in seq(100, 800, 20))
{
	rlist[[cnt]] <- (rollmean(bins[,2], i, fill=0.4, na.pad = TRUE)) * 2
	cnt <- cnt + 1
}
rm <- data.frame(rlist)
colnames(rm) <- seq(100, 800, 20)

pdf("wispgraph_cols.pdf", width=15, height=10)
par(cex=1.8)
plot(density(rm[,1]), lwd=4, type="l", col=rgb(0.1,0.1,0.5,0.8), xlim=c(0,6.5), main="L.rubellus 'core' genome variant distributions, 5-40Kbp Windows", xlab="Polymorphism %", ylab="Density")
for(i in 2:35)
{
	points(density(rm[,i]), type="l", lwd=2, col=rgb(0.1,0.1,0.1,0.25))
}
points(density(rm[,36]), type="l", lwd=4, col=rgb(0.1,0.5,0.1,0.8))
legend("topright", col=c(rgb(0.1,0.1,0.5,0.5), rgb(0.1,0.5,0.1,0.5)), lwd=c(4,4), legend=c("5Kbp Windows", "40Kbp Windows"))
dev.off()

###########################################################

#vars[,2] <- vars[,2] * 2
ravg <- ravg * 2
rd <- density(ravg)
vfull <- vars[vars[,2] > 10,2]
vd <- density(vars[,2])

pdf("25kbpVar.pdf", width=12, height=7)
plot(rd, type="l", lwd=2, col="orange", main="Prob. Density of Polymorphism % in 25Kbp Windows, Full Genome", xlab="Polymorphism %", ylab="Density")
grid()
polygon(rd, col="orange2", border="orange4")
dev.off()

xvals <- vector()
yvals <- vector()
xvals2 <- vector()
yvals2 <- vector()
cnt <- 1

for(i in seq(0,40, 1))
{
	vfull <- vars[vars[,2] <= i,2]
	xvals[cnt] <- i
	yvals[cnt] <- sum(vfull)
	
	cnt <- cnt + 1	
}

cnt <- 1
for(i in seq(1,40, 1))
{
	vfull <- vars[vars[,2] == i,2]
	xvals2[cnt] <- i
	yvals2[cnt] <- sum(vfull)
	
	cnt <- cnt + 1	
}

pdf("50bpVar.pdf", width=12, height=7)
plot(xvals, yvals, type="l", lwd=2, col="pink2", main="SNP Distribution Across 50bp Variant Density Windows", xlab="Polymorphic Bases Per 50bp Window", ylab="Total Genomic Incidence of SNPs")
points(xvals2, yvals2, type="l", lwd=2, col="orange3")
points(xvals2, yvals2, type="h", lwd=4, col="orange3")
abline(h=yvals[40]/2, lty=2)
abline(v=6.65, lty=2)
grid()
legend("topleft", lwd=c(2,2), col=c("pink2", "orange3"), legend=c("Cumulative Frequency", "Counts Per Density Bin"), bg="white")
dev.off()

allvecs <- matrix(0, 40, 46)
cnt2 <- 1

for(b in seq(0.5, 5, 0.1))
{
	b2 <- b - 0.5
	
	vM1 <- vars[ravg < b & ravg > b2 & vars[,2] > 0,2]

	cnt <- 1

	for(i in seq(1,40, 1))
	{
		vfull1 <- vM1[vM1 == i]
		
		allvecs[cnt, cnt2] <- sum(vfull1)
		
		cnt <- cnt + 1	
	}

	allvecs[,cnt2] <- (allvecs[,cnt2] / sum(allvecs[,cnt2])) * 100

	cnt2 <- cnt2 + 1
}

xvals <- 1:40
yvals <- seq(0.5, 5, 0.1)


ymat <- replicate(40, yvals)
xmat <- replicate(46, xvals)

library(lattice)

x <- data.frame(z = as.vector(allvecs), x = as.vector(xmat), y=as.vector(ymat))

png("plotb.png", width=1000, height = 1000)
par(cex=1.4)
wireframe(allvecs,
  xlab = "Variant Density per 50bp", ylab = "25Kbp Window Range as Polymorphism %",
  zlab= "% Of Profile",
  main = "50bp Variant Density Profiles Across 25Kbp Window Ranges",
	scales = list(arrows=FALSE,cex=1,tick.number = 10),
  light.source = c(10,10,10), drape=T,
  screen = list(z = 250, x = -70, y=-30),
  col.regions = colorRampPalette(c("blue", "orange", "pink", "red"))(100)
)
dev.off()

plot(xvals2, yvals2, type="l")
points(xvalsM1, yvalsM1, type="l", col="blue")
points(xvalsM2, yvalsM2, type="l", col="red")

