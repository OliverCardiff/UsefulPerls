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
cwLC_single <- cwLC[rsL2 < 2 & rsL2 > 0,]
cwLC_sparse <- cwLC[rsL2 < 6 & rsL2 > 1,]
cwLC <- cwLC[rsL2 > 5 & rsL2 < ncol(cwLC) - 1,]

lens <- cwLf[,1]
lens <- lens[rsL > 4]
lens_single <- lens[rsL2 < 2 & rsL2 > 0]
lens_sparse <- lens[rsL2 < 6 & rsL2 > 1]
lens <- lens[rsL2 > 5 & rsL2 < ncol(cwLC) - 1]

#d2 <- dist(cwLC)
#hcL <- hclust(d2)
#all_order <- hcL$order
#write.table(hcL$order, file=("cwm_all_order_full.txt"), row.names=FALSE, col.names=FALSE)
all_order <- read.table("cwm_all_order_full.txt", header=FALSE)

#d_sparse <- dist(cwLC_sparse)
#hcL_sparse <- hclust(d_sparse)
#sparse_order <- hcL_sparse$order
#write.table(hcL_sparse$order, file=("cwm_all_order_sparse.txt"), row.names=FALSE, col.names=FALSE)
sparse_order <- read.table("cwm_all_order_sparse.txt", header=FALSE)

#d_single <- dist(cwLC_single)
#hcL_single <- hclust(d_single)
#single_order <- hcL_single$order
#write.table(hcL_single$order, file=("cwm_all_order_single.txt"), row.names=FALSE, col.names=FALSE)
single_order <- read.table("cwm_all_order_single.txt", header=FALSE)

accu <- vector()
accu_sp <- vector()
accu_si <- vector()

cwlList <- cwLC[all_order$V1,]
cwl_sparseList <- cwLC_sparse[sparse_order$V1,]
cwl_singleList <- cwLC_single[single_order$V1,]

lens <- lens[all_order$V1]
lens_sparse <- lens_sparse[sparse_order$V1]
lens_single <- lens_single[single_order$V1]

for(i in 1:(nrow(cwlList) - 1))
{
  accu[i] <- sum(abs(cwlList[i,] - cwlList[i+1,]))
}

for(i in 1:(nrow(cwl_sparseList) - 1))
{
  accu_sp[i] <- sum(abs(cwl_sparseList[i,] - cwl_sparseList[i+1,]))
}

for(i in 1:(nrow(cwl_singleList) - 1))
{
  accu_si[i] <- sum(abs(cwl_singleList[i,] - cwl_singleList[i+1,]))
}

##Running the functions

fullspan <- GetSpans(cwlList, lens, accu, 1, 20)
sparsespan <- GetSpans(cwl_sparseList, lens_sparse, accu_sp, 0, 30)
singlespan <- GetSpans(cwl_singleList, lens_single, accu_si, 0, 20)

namvec_sparse <- GetScaffs(cwl_sparseList, lens_sparse, accu_sp, 0, 30)
namvec_single <- GetScaffs(cwl_singleList, lens_single, accu_si, 0, 20)
namvec <- GetScaffs(cwlList, lens, accu, 1, 20)

write.table(fullspan, file="span_full.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(sparsespan, file="span_sparse.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(singlespan, file="span_single.txt", sep="\t", quote=FALSE, row.names=FALSE)

write.table(namvec, file="namvec_full.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(namvec_sparse, file="namvec_sparse.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(namvec_single, file="namvec_single.txt", sep="\t", quote=FALSE, row.names=FALSE)

head(namvec)

GetSpans <- function(slist, lens, x, b, c)
{
  last <- 999
  spanH <- rep(0,length(x))
  spanCnt <- 0
  indivs <- rep(0,length(x))
  outLen <- rep(0,length(x))

  for(i in 1:length(x))
  {
    if(x[i] < b+1)
    {
      if(last >= b+1)
      {
        spanCnt <- spanCnt + 1
		indivs[spanCnt] <- sum(slist[i,])
      }
      
      outLen[spanCnt] <- outLen[spanCnt] + lens[i]
      spanH[spanCnt] = 1 + spanH[spanCnt]
    }
    last <- x[i]
  }
  spanH <- spanH[1:spanCnt]
  indivs <- indivs[1:spanCnt]
  outLen <- outLen[1:spanCnt]

  outLen <- outLen[spanH > c]
  indivs <- indivs[spanH > c]
  spanH <- spanH[spanH > c]

  data.frame(SpanFrags=spanH, In_Individuals=indivs, Size_bp=outLen) 
}

GetScaffs <- function(slist, lens, x, b, c)
{
  names <- row.names(slist)
  indivs <- rep(0,length(x))
  last <- 999
  spanCnt <- 0
  outNams <- vector()
  outSpans <- vector()
  outLens <- vector()
  inc <- 1
  indVec <- vector()

  for(i in 1:length(x))
  {
    if(x[i] < b+1)
    {
      if(last >= b+1)
      {
        spanCnt <- spanCnt + 1
		indivs[spanCnt] <- sum(slist[i,])
      }
      outNams[inc] <- names[i]
      outSpans[inc] <- spanCnt
      outLens[inc] <- lens[i]
      indVec[inc] <- i

      inc <- inc + 1	
		      
    }
    last <- x[i]
  }
  indivs <- indivs[1:spanCnt]
  outN2 <- vector()
  outS2 <- vector()
  sumlen <- vector()
  indCnt <- vector()
  spanInds <- vector()

  for(j in 1:spanCnt)
  {
     seg <- outNams[outSpans == j]
     sind <- indVec[outSpans == j]
     if(length(seg) > c)
     {
        indCnt <- c(indCnt, rep(indivs[j], length(seg)))
        sumlen <- c(sumlen, rep(sum(outLens[outSpans == j]), length(seg)))
        outS2 <- c(outS2, outSpans[outSpans == j])
        outN2 <- c(outN2, seg)
		spanInds <- c(spanInds, sind)
     }
  }
  data.frame(Names=outN2, ID=outS2, Size_bp=sumlen, In_Individuals=indCnt, spans=slist[spanInds,])
  #data.frame(outNams, outSpans)
}
