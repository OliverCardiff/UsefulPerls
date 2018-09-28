
ratsCWM <- read.table("ratios.txt.csv", sep="\t", header=FALSE)
ratsCF <- read.table("cf_abs.txt.csv", sep="\t", header=FALSE)
ratsDGC <- read.table("dgc_abs.txt.csv", sep="\t", header=FALSE)

pdf("ratios_all.pdf", width=6.5, height=5)

plot((ratsCWM[,2]*100), (ratsCWM[,3]*100), xlim=c(50, 95), ylim=c(40,63), main="RADseq Data Aln. with mito. lineage A/B Genomes", xlab="% alignment rate to Lineage A Genome", ylab="% alignment rate to Lineage B Genome", pch="C", col="darkred")
points((ratsCF[,2]*100), (ratsCF[,3]*100), pch="F", col="purple")
points((ratsDGC[,2]*100), (ratsDGC[,3]*100), pch="D", col="darkblue")
legend("topright", pch=c("C", "F", "D"), col=c("darkred", "purple", "darkblue"), legend=c("CWM Worms", "CF Worms", "DGC Worms"), bg="white")
grid()
dev.off()
