######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


#���ð�
library(limma)
logFCfilter=1       #logFC�ٽ�ֵ
fdrFilter=0.05      #fdr�ٽ�ֵ
setwd("D:\\biowolf\\FerrLnc\\10.diff")   #���ù���Ŀ¼

#�����������ĺ���
Diff=function(expFile=null, project=null){
	#��ȡ�����ļ�
	rt=read.table(expFile, header=T, sep="\t", check.names=F)
	rt=as.matrix(rt)
	rownames(rt)=rt[,1]
	exp=rt[,2:ncol(rt)]
	dimnames=list(rownames(exp), colnames(exp))
	data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
	data=avereps(data)
	data=data[rowMeans(data)>0,]
	
	#������������Ŀ
	group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
	group=sapply(strsplit(group,""), "[", 1)
	group=gsub("2", "1", group)
	conNum=length(group[group==1])       #��������Ʒ��Ŀ
	treatNum=length(group[group==0])     #��������Ʒ��Ŀ
	Type=c(rep(1,conNum), rep(2,treatNum))
	
	#�������
	outTab=data.frame()
	for(i in row.names(data)){
		rt=data.frame(expression=data[i,], Type=Type)
		wilcoxTest=wilcox.test(expression ~ Type, data=rt)
		conGeneMeans=mean(data[i,1:conNum])
		treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
		logFC=log2(treatGeneMeans)-log2(conGeneMeans)
		pvalue=wilcoxTest$p.value
		conMed=median(data[i,1:conNum])
		treatMed=median(data[i,(conNum+1):ncol(data)])
		diffMed=treatMed-conMed
		if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
			outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
		}
	}
	pValue=outTab[,"pValue"]
	fdr=p.adjust(as.numeric(as.vector(pValue)), method="fdr")
	outTab=cbind(outTab, fdr=fdr)
	
	#������л���Ĳ������
	write.table(outTab,file=paste0(project,".all.xls"),sep="\t",row.names=F,quote=F)
	
	#����������
	outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
	write.table(outDiff,file=paste0(project,".diff.xls"),sep="\t",row.names=F,quote=F)
	write.table(outDiff,file=paste0(project,".diff.txt"),sep="\t",row.names=F,quote=F)
	
	#��������ı����ļ�
	diffExp=rbind(ID=colnames(data[as.vector(outDiff[,1]),]),data[as.vector(outDiff[,1]),])
	write.table(diffExp,file=paste0(project,".diffExp.txt"),sep="\t",col.names=F,quote=F)
}

#���ú��������в������
Diff(expFile="FerrGeneExp.txt", project="gene")
Diff(expFile="FerrLncExp.txt", project="lncRNA")


######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056