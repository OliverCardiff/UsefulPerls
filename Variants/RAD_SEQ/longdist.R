gm <- read.table("geno.mat", sep="\t", header=FALSE)
sc10 <- gm[gm[,1] == levels(gm$V1)[1],]
sc10$dist <- sc10$V3 - sc10$V2
hls <- sc10[sc10$dist > 100000 & sc10$V5 > 0.25 & sc10$V4 > 30,]

plot(x=c(0, max(sc10$V3)), y=c(0,0), ylim=c(0,1))
segments(x0=hls$V2, x1=hls$V3, y0=hls$V5, y1=hls$V5, col=rgb(0.1,0.1,0.1,0.005))

hla <- sc10[sc10$dist < 10000 & sc10$V5 > 0.05 & sc10$V4 > 30,]

plot(x=hla$V2, y=hla$V5, col=rgb(0.1,0.1,0.1,0.01), pch="*")

plot(hls$dist, hls$V5, col=rgb(0.1,0.1,0.1,0.05))
