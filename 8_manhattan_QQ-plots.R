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

#build 38, X chr size: 156040895
chrSizes<-c(0,248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468)
#build 37, X chr size: 155270560
chrSizes<-c(0,249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566)
for(i in 2:length(chrSizes)){
  chrSizes[i]<-chrSizes[i-1]+chrSizes[i]
}
filtered_dat = filtered_dat %>% mutate(chromosome=as.numeric(chromosome))
filtered_dat = filtered_dat %>% mutate(gwPos=chrSizes[chr]+pos, chromcols=c("A", "B")[1+chr%%2])

p=ggplot(data=filtered_dat, aes(x=gwPos, y=-log10(P_value), color=chromcols)) +
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
stat[which(abs(stat) < 1e-08)] = NA

stat = sort(stat)
ppoi = ppoints(stat)
ppoi = sort(qchisq(ppoi, df = 1, lower.tail = FALSE))
data_prop = stat[1:ntp]
ppoi_prop = ppoi[1:ntp]

lambda_median = median(data_prop, na.rm = TRUE)/qchisq(0.5, 1)
s = summary(lm(data_prop ~ 0 + ppoi_prop))$coeff
lambda_regression = s[1, 1]
#lambda_se = s[1, 2]

#png(paste(file_path,outprefix,"_QQ_all.png",sep=""),width=800,height=800)
#par(mai=c(1,1,0.05,0.05))
#plot(ppoi, stat, pch=".",cex=1,col="black",
#    ylab=expression(paste("Observed ",group("",chi^2,"")," values",sep="")),
#    xlab=expression(paste("Expected ",group("",chi^2,"")," values",sep=""))) 
    #,xlim=c(0,60),ylim=c(0,60))
#legend("bottomright",parse(text=paste("lambda==",round(lambda,4))))
#abline(a=0,b=1,col="red")
#dev.off()

qq_data=tibble(expected=ppoi, observed=stat)
axis_limit=min(max(qq_data), 60) #don't want to plot anything beyond 60

p = ggplot(qq_data, aes(x=expected, y=observed))+geom_point(shape=".")+
  theme_classic() %+replace% theme(plot.title = element_text(hjust = 0.5, face="bold", size=12))+
  labs(x=expression("Expected "*chi^2*"values"),  y=expression("Observed "*chi^2*"values"))+
  geom_abline(intercept=0, slope=1, linetype=3, colour="red")+
  #xlim(0, axis_limit)+ylim(0, axis_limit)+
  annotate("text", x=0, y = axis_limit, vjust=1, hjust=0, label=paste0("lambda [median] == ", round(lambda_median, 2)), parse=T)+
  annotate("text", x=0, y = axis_limit-2, vjust=1, hjust=0, label=paste0("lambda [regression] == ", round(lambda_regression, 2)), parse=T)

ggsave(paste0(file_path,outprefix,"_QQ.png"), plot=p, dpi=1000, width=15, height=15, units="cm")


