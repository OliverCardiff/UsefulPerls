cwL <- read.table("cwm_all.mat", sep="\t", header=TRUE, row.names=1)

cwLf <- cwL[,colSums(cwL[,1:ncol(cwL)]) > 500000]
csm <- colSums(cwLf[,2:ncol(cwLf)])
csmed <- median(colSums(cwLf[,2:ncol(cwLf)]))

##Read Count Normalisation
for(i in 2:ncol(cwLf))
{
  cwLf[,i] <- (cwLf[,i] / csm[i-1]) * csmed
}

#Coverage Conversion
cwLDiv <- ((cwLf[,2:ncol(cwLf)] * 100)/cwLf[,1]) * 100

#Making binary and filtering empties
rsL <- rowSums(cwLDiv)
cwLDiv <- cwLDiv[rsL > 3,]
cwLC <- cwLDiv
cwLC[cwLC < 4] <- 0
cwLC[cwLC >= 4] <- 1
rsL2 <- rowSums(cwLC)
cwLC <- cwLC[rsL2 > 3 & rsL2 < ncol(cwLC) - 1,]
cwT <- t(cwLC)
cw.pca <- prcomp(cwT, center=TRUE)
PcP <- cw.pca$x[,1:2]
A_points <- PcP[PcP[,1] < 0,]
B_points <- PcP[PcP[,1] > 0,]

pdf("pca_sep.pdf", width=7, height=7)
plot(A_points[,1], A_points[,2], main="PCA separation of A & B Lineages", xlab="PC1", ylab="PC2", xlim=c(-30, 80), ylim=c(-50,30), pch="A", col="darkgreen")
points(B_points[,1], B_points[,2], pch="B", col="darkred")
grid()
legend("topright", bty="n", pch=c("A","B"), col=c("darkgreen", "darkred"), legend=c("Lineage A", "Lineage B"))
dev.off()

As <- row.names(A_points)
Bs <- row.names(B_points)

Nams <- c(Bs, As)
ABs <- c(rep("B", length(Bs)), rep("A", length(As)))

both <- cbind(Nams, ABs)

write.table(both, file="AB_ids.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)
