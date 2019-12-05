library(tidyverse)

#input arguments
file_path="/scratch/DGE/MOPOPGEN/plaw/CRC_GWAS/scotland_meta/output/"
inprefix="studyname"
outprefix="studyname"
flist=Sys.glob(paste0(file_path,inprefix,"chr*.txt.gz"))

# read in all the files individually, using readr
all_dat = flist %>% purrr::map(read_delim, delim=" ") %>% reduce(rbind)

#####################
# Manhattan plot
#####################

filtered_dat=all_dat[which(all_dat$P_value<0.05),]

chrSizes<-c(0,249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566)
for(i in 2:length(chrSizes)){
  chrSizes[i]<-chrSizes[i-1]+chrSizes[i]
}
filtered_dat = filtered_dat %>% mutate(gwPos=chrSizes[chr]+pos, chromcols=c("A", "B")[1+chr%%2])

p=ggplot(data=filtered_dat, aes(x=gwPos, y=-log10(P_value), color=chromcols, alpha=exp(abs(beta)))) +
  geom_point( size=1) +
  theme_classic(base_size = 16) +
  labs(x="Chromosome & position",y=expression('-log'[10]*'(P)')) +
  scale_x_continuous(labels=1:22, breaks=(chrSizes[1:(length(chrSizes)-1)]+chrSizes[2:length(chrSizes)])/2) +
  scale_colour_manual(values = c("#616365", "#F9A100"),guide=FALSE) +
  geom_hline(yintercept = -log10(5e-8), lty=2, color="red")

ggsave(paste0(file_path,outprefix,"_manhattan.png"), p, height=7, width=15)


#####################
# QQ plot
#####################
inDat = all_dat %>% mutate(stat=(frequentist_add_beta_1/frequentist_add_se_1)^2)

#proportion of data to use
proportion=0.9
stat=inDat %>% pull(stat)
stat=stat[which(!is.na(stat))]

ntp = round(proportion * length(stat))
stat = qchisq(stat, 1, lower.tail = FALSE)
stat[which(abs(stat) < 1e-08)] = NA

stat = sort(stat)
ppoi = ppoints(stat)
ppoi = sort(qchisq(ppoi, df = df, lower.tail = FALSE))
data_prop = stat[1:ntp]
ppoi_prop = ppoi[1:ntp]

lambda = median(data_prop, na.rm = TRUE)/qchisq(0.5, 1)

png(paste(file_path,outprefix,"_QQ_all.png",sep=""),width=800,height=800)
par(mai=c(1,1,0.05,0.05))
plot(ppoi, stat, pch=".",cex=1,col="black",
    ylab=expression(paste("Observed ",group("",chi^2,"")," values",sep="")),
    xlab=expression(paste("Expected ",group("",chi^2,"")," values",sep=""))) 
    #,xlim=c(0,60),ylim=c(0,60))
legend("bottomright",parse(text=paste("lambda==",round(lambda,4))))
abline(a=0,b=1,col="red")
dev.off()






print(estlambda(inDat$frequentist_add_pvalue, proportion = 0.9, method="median"))
print(estlambda(inDat$frequentist_add_pvalue, proportion = 0.9))