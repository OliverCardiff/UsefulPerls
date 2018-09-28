cwH <- read.table("cwm_high_mat.txt", sep="\t", header=TRUE, row.names=1)
cwL <- read.table("cwm_low_mat.txt", sep="\t", header=TRUE, row.names=1)

cwLf <- cwL[,colSums(cwL[,2:55]) > 260000]
cwHf <- cwH[,colSums(cwH[,2:55]) > 300000]

cwHDiv <- ((cwHf[,2:39] * 100)/cwHf[,1]) * 100
cwLDiv <- ((cwLf[,2:39] * 100)/cwLf[,1]) * 100

rsL <- rowSums(cwLDiv)
rsH <- rowSums(cwHDiv)

cwLDiv <- cwLDiv[rsL < 10000 & rsL > 5,]
cwHDiv <- cwHDiv[rsH < 10000 & rsH > 5,]

plot(density(rowSums(cwLDiv)), xlim=c(0,200))
points(density(rowSums(cwHDiv)), type="l")

cwLC <- cwLDiv
cwHC <- cwHDiv
cwLC[cwLC < 4] <- 0
cwHC[cwHC < 4] <- 0

cwLC[cwLC >= 4] <- 1
cwHC[cwHC >= 4] <- 1

d <- dist(cwHC)
hcH <- hclust(d)
#plot(hc, cex=0.1)

d2 <- dist(cwLC)
hcL <- hclust(d2)
#plot(hc, cex=0.1)

accu <- vector()

cwlList <- cwLC[hcL$order,]

cwhList<- cwHC[hcH$order,]

#cwhList <- cwHC
#cwlList <- cwLC


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

for(i in 1:length(accu2))
{
	sumZeroH[accu2[i]+1] <- sumZeroH[accu2[i]+1] + 1
}

spanL <- GetSpans(accu, 1,4)
spanH <- GetSpans(accu2,1,4)

staZL <- sumZero/sum(sumZero) * 100
staZH <- sumZeroH/sum(sumZeroH) * 100

plot(0:34, staZL[1:35], type="l", col=rgb(0.4,0.22,0.66,0.5), lwd=2)
points(0:34, staZH[1:35], type="l", col=rgb(0.2,0.66,0.21,0.5), lwd=2)

span_mat <- matrix(0,7,10)
for(i in 1:7)
{
  spn <- GetSpans(accu, i,20)
  spn <- spn[rev(order(spn))]
  span_mat[i,] <- spn[1:10]
}
pdf("flex_spans.pdf", width=5, height=4)
par(cex=0.9)
barplot(span_mat, main = "Linkage blocks sizes by model flexibility", 
        beside=TRUE, xlab = "Top 10 Blocks Sizes", ylab="Scaffold Count", 
        legend=c("Flex = 0", "Flex = 1", "Flex = 2", "Flex = 3", 
                 "Flex = 4", "Flex = 5", "Flex = 6"))
dev.off()

GetSpans <- function(x, b, c)
{
	last <- 999
	spanH <- rep(0,length(x))
	spanCnt <- 0

	for(i in 1:length(x))
	{
		if(x[i] < b+1)
		{
			if(last >= b+1)
			{
				spanCnt <- spanCnt + 1
			}

			spanH[spanCnt] = 1 + spanH[spanCnt]
		}
		last <- x[i]
	}
	spanH <- spanH[1:spanCnt]
	
	spanH[spanH > c]
}
