library("geneplotter")
library("tidyverse")

#path to input data
path=""
#root of .bed/.bim/.fam files
expername=""

imiss = read_table(paste0(path,"/",expername,".imiss"))
imiss = imiss %>% mutate(logF_MISS = log10(F_MISS))
het = read_table(paste0(path,"/",expername,".het"))
het = het %>% mutate(meanHet = (`N(NM)` -  `O(HOM)`)/`N(NM)` )

colors = densCols(imiss$logF_MISS,het$meanHet, nbin=128, bandwidth = c(0.5,0.5))
pdf(paste0(path,"/",expername,".imiss-vs-het.pdf"))
plot(imiss$logF_MISS,het$meanHet, col=colors, pch=20, ylim=c(0,0.5), xlim=c(-6,0), xlab="Proportion of missing genotypes", ylab="Heterozygosity rate",axes=F)
axis(2,at=c(0,0.05,0.10,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5),tick=T)
axis(1,at=-6:0, labels=10**(-6:0))
m = mean(het$meanHet)
stdev = sd(het$meanHet)

missingness_threshold=0.05
n_stdevs=3
abline(h=m-(n_stdevs*stdev),col="RED",lty=2)
abline(h=m+(n_stdevs*stdev),col="RED",lty=2)
abline(v=log10(missingness_threshold), col="RED", lty=2)
dev.off()

missing = imiss %>% filter($logF_MISS>log10(missingness_threshold)) %>% select(FID,IID)
heterozygous = het %>% filter(meanHet<(m-n_stdevs*stdev) | meanHet>(m+n_stdevs*stdev)) %>% select(FID,IID)
exclude = union(missing, heterozygous)

write_delim(exclude, paste0(path,"/fail-imisshet-qc.txt"), col_names=F)
