cwL <- read.table("mat_triple/cwm_gen.mat", sep="\t", header=TRUE, row.names=1)
cfL <- read.table("mat_triple/cf_gen.mat", sep="\t", header=TRUE, row.names=1)
dgL <- read.table("mat_triple/dgc_gen.mat", sep="\t", header=TRUE, row.names=1)

cwCom <- findcommon(cwL, 500000, 12)
cfCom <- findcommon(cfL, 200000, 6)
dgCom <- findcommon(dgL, 400000, 10)

allnam <- c(row.names(dgCom), row.names(cfCom), row.names(cwCom))

write.table(allnam, file="scafflist2.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)

findcommon <- function(cwX, thresh, flex)
{
	cwLf <- cwX[,colSums(cwX[,1:ncol(cwX)]) > thresh]
	csm <- colSums(cwLf[,2:ncol(cwLf)])
	csmed <- median(colSums(cwLf))	

	##Read Count Normalisation
	for(i in 2:ncol(cwLf))
	{
	  cwLf[,i] <- (cwLf[,i] / csm[i-1]) * csmed
	}

	#Coverage Conversion
	cwLDiv <- ((cwLf[,2:ncol(cwLf)] * 100)/cwLf[,1]) * 100

	#Making binary and filtering empties
	rsL <- rowSums(cwLDiv)
	cwLDiv <- cwLDiv[rsL > 4,]
	cwLC <- cwLDiv
	cwLC[cwLC < 4] <- 0
	cwLC[cwLC >= 4] <- 1
	rsL2 <- rowSums(cwLC)
	cwLC <- cwLC[rsL2 > ((ncol(cwLDiv) - 1) - flex),]
	cwLC
}

dobars <- function(cw)
{
	barplot(colSums(cw[,2:ncol(cw)]))
}
