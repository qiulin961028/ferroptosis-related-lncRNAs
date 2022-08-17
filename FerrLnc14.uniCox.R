######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056

#install.packages("survival")


library(survival)      #���ð�
pFilter=0.05           #�����Թ�������
setwd("D:\\biowolf\\FerrLnc\\14.uniCox")     #���ù���Ŀ¼
rt=read.table("expTime.txt", header=T, sep="\t", check.names=F, row.names=1)     #��ȡ�����ļ�
rt$futime=rt$futime/365      #���浥λ�ĳ���

#������cox����
outTab=data.frame()
sigGenes=c("futime","fustat")
for(gene in colnames(rt[,3:ncol(rt)])){
	cox=coxph(Surv(futime, fustat) ~ rt[,gene], data = rt)
	coxSummary = summary(cox)
	coxP=coxSummary$coefficients[,"Pr(>|z|)"]
		
	if(coxP<pFilter){
		sigGenes=c(sigGenes,gene)
		outTab=rbind(outTab,
			       cbind(gene=gene,
			       HR=coxSummary$conf.int[,"exp(coef)"],
			       HR.95L=coxSummary$conf.int[,"lower .95"],
			       HR.95H=coxSummary$conf.int[,"upper .95"],
				   pvalue=coxP) )
	}
}

#��������ؽ��
write.table(outTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)
surSigExp=rt[,sigGenes]
surSigExp=cbind(id=row.names(surSigExp),surSigExp)
write.table(surSigExp,file="uniSigExp.txt",sep="\t",row.names=F,quote=F)


############����ɭ��ͼ����############
bioForest=function(coxFile=null, forestFile=null){
	#��ȡ�����ļ�
	rt <- read.table(coxFile,header=T,sep="\t",check.names=F,row.names=1)
	gene <- rownames(rt)
	hr <- sprintf("%.3f",rt$"HR")
	hrLow  <- sprintf("%.3f",rt$"HR.95L")
	hrLow[hrLow<0.001]=0.001
	hrHigh <- sprintf("%.3f",rt$"HR.95H")
	Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
	pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
		
	#���ͼ��
	height=(nrow(rt)/15)+5
	pdf(file=forestFile, width=8, height=height)
	n <- nrow(rt)
	nRow <- n+1
	ylim <- c(1,nRow)
	layout(matrix(c(1,2),nc=2),width=c(3,2.5))
		
	#����ɭ��ͼ��ߵ���Ϣ
	xlim = c(0,3)
	par(mar=c(4,2.5,2,1))
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
	text.cex=0.8
	text(0,n:1,gene,adj=0,cex=text.cex)
	text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,adj=1)
	text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,adj=1,)
		
	#����ɭ��ͼ
	par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
	LOGindex=10 
	hrLow = log(as.numeric(hrLow),LOGindex)
	hrHigh = log(as.numeric(hrHigh),LOGindex)
	hr = log(as.numeric(hr),LOGindex)
	xlim = c(floor(min(hrLow,hrHigh)),ceiling(max(hrLow,hrHigh)))
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
	arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
	abline(v=log(1,LOGindex),col="black",lty=2,lwd=2)
	boxcolor = ifelse(as.numeric(hr) > log(1,LOGindex), "red", "blue")
	points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
	a1 = axis(1,labels=F,tick=F)
	axis(1,a1,LOGindex^a1)
	dev.off()
}
############����ɭ��ͼ����############

#����ɭ��ͼ
bioForest(coxFile="uniCox.txt", forestFile="forest.pdf")


######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056