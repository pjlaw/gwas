path="/scratch/DGE/MOPOPGEN/plaw/CRC_GWAS/oncoarray/imputation/unphased/redoQC/"
expername="crc_onco_caco"

imiss=read.table(paste0("/run/user/1000/gvfs/sftp:host=davros",path,expername,".imiss"),h=T)
imiss$logF_MISS = log10(imiss[,6])
het=read.table(paste0("/run/user/1000/gvfs/sftp:host=davros",path,expername,".het"),h=T)
het$meanHet = (het$N.NM. - het$O.HOM.)/het$N.NM.
library("geneplotter")
colors  <- densCols(imiss$logF_MISS,het$meanHet, nbin=128, bandwidth = c(0.5,0.5))
# pdf(paste0("/run/user/1000/gvfs/sftp:host=davros",path,expername,".imiss-vs-het.pdf"))
plot(imiss$logF_MISS,het$meanHet, col=colors, pch=20, ylim=c(0,0.5), xlim=c(-6,0), xlab="Proportion of missing genotypes", ylab="Heterozygosity rate",axes=F)
axis(2,at=c(0,0.05,0.10,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5),tick=T)
axis(1,at=-6:0, labels=10**(-6:0))
m=mean(het$meanHet)
stdev=sd(het$meanHet)
abline(h=m-(3*stdev),col="RED",lty=2)
abline(h=m+(3*stdev),col="RED",lty=2)
abline(v=log10(0.03), col="RED", lty=2)
# dev.off()
# 
# missing=imiss[which(imiss$logF_MISS>log10(0.05)),1]
# hetero1=het[which(het$meanHet<(m-3*stdev)),1]
# hetero2=het[which(het$meanHet>(m+3*stdev)),1]
# exclude=union(missing, hetero1)
# exclude=union(exclude, hetero2)
# fileConn<-file(paste0("/run/user/1000/gvfs/sftp:host=davros",path,"fail-imisshet-qc.txt"),"w")
# for (e in exclude){
#   write(paste(e,e), fileConn, append=T)
# }
# close(fileConn)