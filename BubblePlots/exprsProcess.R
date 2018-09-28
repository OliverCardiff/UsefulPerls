names <- colnames(exprs)

png("mqqexprs.png", width=800, height=500)
plot(1:5, meanLvl, type="b", ylim=c(0,0.25), ylab="Per Base Read Coverage Probability",xlab="Gene Expression Ascending 20% Size Rank Groups", lwd=2, col="orangered", main="Intra Genic Me-DIP Read Alignment vs Gene Expression")
points(1:5, q25, type="b", lwd=2, col="aquamarine")
points(1:5, q50, type="b", lwd=2, col="cyan3")
points(1:5, q75, type="b", lwd=2, col="dodgerblue")
legend("topleft", col=c("orangered", "aquamarine", "cyan3", "dodgerblue"), lwd=c(2,2,2,2), legend = c("Means", "Q25", "Median", "Q75"))
dev.off()

png("loess.png", width=800, height=500)
predd <- lowess(bigt[,2], f = 1/3)
plot(bigt[,4], bigt[,2], type="p", pch="*", ylim=c(0.02,0.25), xlab="Gene Expression Ranks", ylab="Per Base Me-DIP Read Alignment Probability", main="LOESS Trend Line of Gene Expression Rank vs Methylation")
lines(predd, col=2, lwd=4)
dev.off()

## genome stuff start here

png("symbs.png", height=1000, width=1000)
symbols(genome2[,4], genome2[,3], circles=genome2[,1], ylim=c(25,180), inches = 1, xlab="% GC Content", ylab="Scaffold Read Coverage", main="Scaffolds OVER 20k base pairs in A.gracilis Genome")
abline(h=50, col="red")
abline(h=100, col="green")
legend("topright", col=c("red", "green"), lty=c(1,1), lwd=c(1,1), legend=c("Half Expected Coverage Level", "Expected Full Coverage Level"))
dev.off()

png("symb2s.png", height=1000, width=1000)
symbols(genome3[,4], genome3[,3], circles=genome3[,1], ylim=c(10,120), inches = 1, xlab="% GC Content", ylab="Scaffold Read Coverage", main="Scaffolds UNDER 20k base pairs in A.gracilis Genome (zoomed in)")
abline(h=50, col="red")
abline(h=100, col="green")
legend("topright", col=c("red", "green"), lty=c(1,1), lwd=c(1,1), legend=c("Half Expected Coverage Level", "Expected Full Coverage Level"))
dev.off()

png("symbsBoth.png", height=1000, width=1000)
symbols(filt4[,4], filt4[,3], circles=filt4[,1], ylim=c(20,160), xlim=c(18,62), inches = 1.5, fg="darkgoldenrod", xlab="% GC Content", ylab="Scaffold Read Coverage", main="Scaffolds: Genuine and Alien in the A.gracilis Assembly")
symbols(ex1[,4], ex1[,3], circles=ex1[,1], inches = 0.5, fg="firebrick", add=TRUE)
legend("topright", col=c("darkgoldenrod", "firebrick"), lty=c(1,1), lwd=c(2,2), pch=3, legend=c("Core Genome: 583Mb", "Duplication and Alien: 241Mb"))
dev.off()

png("genome.png", height=1000, width=1000)
symbols(filt4[,4], filt4[,3], circles=filt4[,1], ylim=c(60,140), xlim=c(35,45), inches = 1.5, fg="darkgoldenrod", xlab="% GC Content", ylab="Scaffold Read Coverage", main="Scaffolds: Core A.gracilis Assembly")
dev.off()

filt1 <- genome[genome[,3] > 75,]
filt2 <- filt1[filt1[,4] < 43,]
filt3 <- filt2[filt2[,4] > 37,]
filt4 <- filt3[filt3[,3] < 120,]

ex1 <- genome[genome[,3] < 75 | genome[,4] > 43 | genome[,4] < 37 | genome[,3] > 120,]
ex2 <- ex1[ex1[,4] > 43,]
ex3 <- ex2[ex2[,4] < 37,]
ex4 <- ex3[ex3[,3] > 120,]

png("genomeNxt.png", height=2000, width=2000)
symbols(genome[,4], genome[,3], circles=genome[,1], ylim=c(20,160), xlim=c(34,56), inches = 1.5, fg=rgb(0.1,0.1,0.4), bg=rgb(0.1,0.1,0.4,0.05), xlab="% GC Content", ylab="Scaffold Read Coverage", main="Scaffolds: A.gracilis Assembly")
symbols(bus_dup[,5], bus_dup[,4], stars=cbind(bus_dup[,2], bus_dup[,2], bus_dup[,2]), inches=0.2, fg=rgb(0.8,0.2,0.2), bg=rgb(0.8,0.2,0.2,0.8), add=TRUE)
symbols(bus_comp[,5], bus_comp[,4], stars=cbind(bus_comp[,2], bus_comp[,2], bus_comp[,2]), inches=0.2, fg=rgb(0.2,0.8,0.2), bg=rgb(0.2,0.8,0.2,0.8), add=TRUE)
segments(netDF[,1], netDF[,2], netDF[,3], netDF[,4], lty=1, lwd=1, col=rgb(0.1,0.1,0.1,0.5))
dev.off()

png("genomeNxt3.png", height=2000, width=2000)
symbols(genome2[,4], genome2[,3], circles=genome2[,1], ylim=c(20,350), xlim=c(34,56), inches = 1.5, fg=rgb(0.1,0.1,0.4), bg=rgb(0.1,0.1,0.4,0.05), xlab="% GC Content", ylab="Scaffold Read Coverage", main="Scaffolds: A.gracilis Assembly")
symbols(bus_dup2[,5], bus_dup2[,4], stars=cbind(bus_dup2[,2], bus_dup2[,2], bus_dup2[,2]), inches=0.2, fg=rgb(0.8,0.2,0.2), bg=rgb(0.8,0.2,0.2,0.8), add=TRUE)
symbols(bus_comp2[,5], bus_comp2[,4], stars=cbind(bus_comp2[,2], bus_comp2[,2], bus_comp2[,2]), inches=0.2, fg=rgb(0.2,0.8,0.2), bg=rgb(0.2,0.8,0.2,0.8), add=TRUE)
segments(netDF2[,1], netDF2[,2], netDF2[,3], netDF2[,4], lty=1, lwd=1, col=rgb(0.1,0.1,0.1,0.5))
dev.off()

png("genomeNxt4.png", height=2000, width=2000)
symbols(bus_both[,5], bus_both[,4], circles=bus_both[,2], ylim=c(20,160), xlim=c(34,56), inches = 1.5, fg=rgb(0.1,0.1,0.4), bg=rgb(0.1,0.1,0.4,0.05), xlab="% GC Content", ylab="Scaffold Read Coverage", main="Scaffolds: A.gracilis Assembly")
symbols(bus_dup2[,5], bus_dup2[,4], stars=cbind(bus_dup2[,2], bus_dup2[,2], bus_dup2[,2]), inches=0.2, fg=rgb(0.8,0.2,0.2), bg=rgb(0.8,0.2,0.2,0.8), add=TRUE)
symbols(bus_comp2[,5], bus_comp2[,4], stars=cbind(bus_comp2[,2], bus_comp2[,2], bus_comp2[,2]), inches=0.2, fg=rgb(0.2,0.8,0.2), bg=rgb(0.2,0.8,0.2,0.8), add=TRUE)
segments(netDF2[,1], netDF2[,2], netDF2[,3], netDF2[,4], lty=1, lwd=1, col=rgb(0.1,0.1,0.1,0.5))
dev.off()