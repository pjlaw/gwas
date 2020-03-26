path=""
expername=""

data=read.table(paste0(path,"/", expername, ".pca.evec"),h=F,skip=1)

lastcol=ncol(data)
cont=which(data[,lastcol]=="Control")
case=which(data[,lastcol]=="Case")
CEU=which(data[,lastcol]=="3")
CHB=which(data[,lastcol]=="4")
JPT=which(data[,lastcol]=="5")
YRI=which(data[,lastcol]=="6")
pdf(paste0(path,"/pca-ancestry-plot.pdf"))

plot(0,0,pch="", xlim=c(min(data$V2), max(data$V2)), ylim=c(min(data$V3), max(data$V3)), xlab="principal component 1", ylab="principal component 2")

points(data$V2[JPT],data$V3[JPT],pch=20,col="PURPLE")
points(data$V2[CHB],data$V3[CHB],pch=20,col="PURPLE")
points(data$V2[YRI],data$V3[YRI],pch=20,col="GREEN")
points(data$V2[CEU],data$V3[CEU],pch=20,col="RED")
points(data$V2[cont],data$V3[cont],pch=20,col="orange")
points(data$V2[case],data$V3[case],pch="+",col="BLACK")

legend("topleft", legend=c("CHB/JPT", "YRI", "CEU", "Controls", "Cases"), pch=c(20, 20, 20, 20, 3), col = c("purple", "green", "red", "orange", "black"), cex=0.8, y.intersp=0.8)

dev.off()
# 
# edit the limits as appropriate and run this to identify ca/co samples that are out of bounds
# exclude=as.character(data[which((data[,lastcol]=="Case" | data[,lastcol]=="Control") & (data$V2<=0.013 | data$V2>=0.019 | data$V3<=0.07)), 1])
# exclude=gsub(":", " ", exclude)
# writeLines(exclude, paste0(path, "/fail-pca.txt"))
