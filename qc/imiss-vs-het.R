#path to input data
path=""
#root of .bed/.bim/.fam files
expername=""

imiss=read.table(paste0(path,"/",expername,".imiss"),h=T)
imiss$logF_MISS = log10(imiss[,6])
het=read.table(paste0(path,"/",expername,".het"),h=T)
het$meanHet = (het$N.NM. - het$O.HOM.)/het$N.NM.
colors  <- densCols(imiss$logF_MISS,het$meanHet, nbin=128, bandwidth = c(0.5,0.5))
pdf(paste0(path,"/",expername,".imiss-vs-het.pdf"))
plot(imiss$logF_MISS,het$meanHet, col=colors, pch=20, ylim=c(0,0.5), xlim=c(-6,0), xlab="Proportion of missing genotypes", ylab="Heterozygosity rate",axes=F)
axis(2,at=c(0,0.05,0.10,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5),tick=T)
axis(1,at=-6:0, labels=10**(-6:0))
m=mean(het$meanHet)
stdev=sd(het$meanHet)

missingness_threshold=0.05
n_stdevs=3
abline(h=m-(n_stdevs*stdev),col="RED",lty=2)
abline(h=m+(n_stdevs*stdev),col="RED",lty=2)
abline(v=log10(missingness_threshold), col="RED", lty=2)
dev.off()

missing=imiss[which(imiss$logF_MISS>log10(missingness_threshold)),1]
hetero1=het[which(het$meanHet<(m-n_stdevs*stdev)),1]
hetero2=het[which(het$meanHet>(m+n_stdevs*stdev)),1]
exclude=union(missing, hetero1)
exclude=union(exclude, hetero2)
fileConn<-file(paste0(path,"/fail-imisshet-qc.txt"),"w")
for (e in exclude){
  write(paste(e,e), fileConn, append=T)
}
close(fileConn)