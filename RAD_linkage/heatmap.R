
cwL <- read.table("cwm_all.mat", sep="\t", header=TRUE, row.names=1)

cwLf <- cwL[,colSums(cwL[,1:ncol(cwL)]) > 500000]
csm <- colSums(cwLf[,2:ncol(cwLf)])
csmed <- median(colSums(cwLf[,2:ncol(cwLf)]))

cwLf[cwLf < 4] <- 0

##Read Count Normalisation
for(i in 2:ncol(cwLf))
{
  cwLf[,i] <- (cwLf[,i] / csm[i-1]) * csmed
}

#Coverage Conversion
cwLDiv <- ((cwLf[,2:ncol(cwLf)] * 100)/cwLf[,1]) * 100
rsm <- rowSums(cwLDiv)
cwLDiv <- cwLDiv[rsm < 200 & rsm > 1,]
rsm2 <- rowSums(cwLDiv)

cwLC <- cwLDiv
cwLC[cwLC < 0.5] <- 0
cwLC[cwLC >= 0.5] <- 1

rsm3 <- rowSums(cwLC)

cwLC <- cwLC[rsm3 > 6,]

d2 <- dist(cwLC)
hcL <- hclust(d2)
all_order <- hcL$order
cwlList <- cwLC[all_order,]

d3 <- dist(t(cwLC))
hcL3 <- hclust(d3)
col_order <- hcL3$order
cwlList <- cwlList[,col_order]

cwlList <- cwlList[, c(3:ncol(cwlList),1,2)]

cwlList <- cwlList[,c(1:8,39,40,38,37,35,36, 9:34)]

fx <- vector()
fy <- vector()

for(i in 1:ncol(cwlList))
{
	x <- rep(i, nrow(cwlList))
	y <- cwlList[,i]
	yx <- cbind(y, c(1:nrow(cwlList)))
	y <- yx[yx[,1] == 1, 2] 
	x <- x[yx[,1] == 1]	
	fy <- c(fy, y)
	fx <- c(fx, x)
}

bsm <- rowSums(cwlList[,1:8])
bsm2 <- rowSums(cwlList[,9:ncol(cwlList)])

bd <- bsm - bsm2

Bs <- cwlList[rev(order(bd)),]

fx2 <- vector()
fy2 <- vector()

for(i in 1:ncol(Bs))
{
	x <- rep(i, nrow(Bs))
	y <- Bs[,i]
	yx <- cbind(y, c(1:nrow(Bs)))
	y <- yx[yx[,1] == 1, 2] 
	x <- x[yx[,1] == 1]	
	fy2 <- c(fy2, y)
	fx2 <- c(fx2, x)
}

png("pop_struct.png", width=550, height=1000)
plot(fx, fy, pch="-", col=rgb(0.1,0.1,0.1,0.3), cex=2.8, main="RADseq Presence/Absence Patterns, Clustered by Row", xlab="Individuals", ylab="Scaffolds with Mapped RADseq data")
dev.off()

png("pop_struct_s.png", width=550, height=1000)
plot(fx2, fy2, pch="-", col=rgb(0.1,0.1,0.1,0.3), cex=2.8, main="RADseq Presence/Absence Patterns, Sorted by Lineage-Difference", xlab="Individuals", ylab="Scaffolds with Mapped RADseq data")
dev.off()

png("counts_pop.png", width=1000, height=500)
barplot(colSums(cwlList), names=abs, ylab="Scaffolds with Aligned Stacks", xlab="Individuals by Lineage", main="L. rubellus lineages reference map rates")
abline(v=9.7, lwd=4)
dev.off()

abs <- c("B", "B", "B", "B", "B", "B", "B", "B", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A")

write.table(hcL$order, file=("cwm_all_order_full.txt"), row.names=FALSE, col.names=FALSE)
all_order <- read.table("cwm_all_order_full.txt", header=FALSE)
