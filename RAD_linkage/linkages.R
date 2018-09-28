cwH <- read.table("cwm_high_mat.txt", sep="\t", header=TRUE, row.names=1)
cwL <- read.table("cwm_low_mat.txt", sep="\t", header=TRUE, row.names=1)

cwHDiv <- ((cwH[,2:55] * 100)/cwH[,1]) * 100
cwLDiv <- ((cwL[,2:55] * 100)/cwL[,1]) * 100

rsL <- rowSums(cwLDiv)
rsH <- rowSums(cwHDiv)

cwLDiv <- cwLDiv[rsL < 1000 & rsL > 6,]
cwHDiv <- cwHDiv[rsH < 1000 & rsH > 6,]

plot(density(rowSums(cwLDiv)), xlim=c(0,200))
points(density(rowSums(cwHDiv)), type="l")

cwLC <- cwLDiv
cwHC <- cwHDiv
cwLC[cwLC > 0] <- 1
cwHC[cwHC > 0] <- 1

d <- dist(cwHC)
hcH <- hclust(d)
plot(hc, cex=0.1)

d <- dist(cwLC)
hcL <- hclust(d)
plot(hc, cex=0.1)

accu <- vector()

cwlList <- cwLC[hcL$order,]

cwhList<- cwHC[hcH$order,]

cwhList <- cwHC
cwlList <- cwLC


for(i in 1:(nrow(cwlList) - 1))
{
	accu[i] <- sum(abs(cwlList[i,] - cwlList[i+1,]))
}

sumZero <- rep(0,54)

for(i in 1:length(accu))
{
	sumZero[accu[i]+1] <- sumZero[accu[i]+1] + 1
}

accu2 <- vector()

for(i in 1:(nrow(cwhList) - 1))
{
	accu2[i] <- sum(abs(cwhList[i,] - cwhList[i+1,]))
}

sumZeroH <- rep(0,54)
oldAc <- 999

for(i in 1:length(accu2))
{
	sumZeroH[accu2[i]+1] <- sumZeroH[accu2[i]+1] + 1
}

staZL <- sumZero/sum(sumZero) * 100
staZH <- sumZeroH/sum(sumZeroH) * 100

plot(0:34, staZL[1:35], type="l", col=rgb(0.4,0.22,0.66,0.5), lwd=2)
points(0:34, staZH[1:35], type="l", col=rgb(0.2,0.66,0.21,0.5), lwd=2)

cwLf <- cwL[,colSums(cwL[,2:55]) > 260000]
cwHf <- cwH[,colSums(cwH[,2:55]) > 300000]

