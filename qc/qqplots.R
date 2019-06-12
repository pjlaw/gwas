#pvals <- read.table("/run/user/1000/gvfs/sftp:host=dalek/scratch/cancgene/plaw/CRC/VQ/data/VQ_allchr_QC_geno_maf.hwe", header=T)

subsetP = pvals[which(pvals$TEST=="ALL"),9]
subsetP = subsetP[which(subsetP>1e-5)]
observed <- sort(subsetP)
# observed <- sort(pvals$P)

lobs <- -(log10(observed))

expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))

#pdf("qqplot.pdf", width=6, height=6)
plot(c(0,7), c(0,7), col="red", lwd=1, pch='.', type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,7), ylim=c(0,7), las=1, xaxs="i", yaxs="i", bty="l")
points(lexp, lobs, pch=23, cex=.4, bg="black") 
#dev.off()