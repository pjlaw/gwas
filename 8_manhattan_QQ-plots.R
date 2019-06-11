#-----input arguments
	args<-commandArgs(TRUE)
	inFile<-args[1]
	outPrefix<-args[2]

#-----read the data
	inDat<-read.table(inFile,sep="\t",header=T)
	
#-----Manhattan plot
	chrSizes<-c(0,249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566)
	for(i in 2:length(chrSizes)){
		chrSizes[i]<-chrSizes[i-1]+chrSizes[i]
	}
	gwPos<-chrSizes[inDat[,3]]+inDat[,4]
	idx.plot<-which(inDat[,40]<5e-2)
	pdf(paste(outPrefix,"_manhattan.pdf",sep=""),width=25,height=7)
	par(mai=c(0.5,1.25,0.05,0.05))
  
	plot(gwPos[idx.plot],-log10(inDat[idx.plot,40]),pch=".",cex=2,col=c("turquoise4","steelblue4")[1+inDat[idx.plot,3]%%2],xlab="chromosome & position",ylab="-log10(P)",axes=F,ylim=c(0,20))
	axis(side=2)
	axis(side=1,at=(chrSizes[1:(length(chrSizes)-1)]+chrSizes[2:length(chrSizes)])/2,labels=paste("chr",1:22,sep=""),hadj=0.5)
	dev.off()

#-----inflation factor (based on a random subset of 1e6 variants)
	#idx.sample<-sample(1:nrow(inDat),size=min(1e6,nrow(inDat)),replace=F)
	#inDat<-inDat[idx.sample,]
	stat<-(inDat[,42]/inDat[,43])^2 #beta/stderr
	Per<-0.9
	len<-length(stat)
	lambda<-median(sort(stat)[1:(Per*len)])/median(qchisq(ppoints(len)[1:(Per*len)],1))
	pdf(paste(outPrefix,"_QQ.pdf",sep=""),width=8,height=8)
	par(mai=c(1,1,0.05,0.05))
	plot(qchisq(ppoints(len),df=1),sort(stat),pch=".",cex=1,col="black",ylab=expression(paste("observed ",group("",chi^2,"")," values",sep="")),xlab=expression(paste("expected ",group("",chi^2,"")," values",sep="")),xlim=c(0,60),ylim=c(0,60))
	legend("topleft",parse(text=paste("lambda==",round(lambda,4))))
	abline(a=0,b=1,col="red")
	dev.off()
