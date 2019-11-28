path=""
expername=""

data=read.table(paste0(path, "/", expername, ".genome"),h=T)
pdf(path,"/IBD-hist.pdf")
hist(data$PI_HAT,ylim=c(0,100),col="RED",breaks=100,xlab="Estimated mean pairwise IBD",main="")
out=which(data$PI_HAT>18.5)
write.table("fail-IBD-check.txt",data[out,],col.names=F,row.names=F,sep="")
dev.off()
