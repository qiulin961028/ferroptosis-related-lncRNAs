######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library(limma)            #���ð�
expFile="symbol.txt"      #���������ļ�
geneFile="gene.txt"       #�����б��ļ�
setwd("D:\\biowolf\\FerrLnc\\08.FerrExp")    #���ù���Ŀ¼

#��ȡ�����ļ����������ݽ��д���
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#��ȡ����ı�����
gene=read.table(geneFile, header=F, check.names=F, sep="\t")
sameGene=intersect(as.vector(gene[,1]), rownames(data))
geneExp=data[sameGene,]

#������
out=rbind(ID=colnames(geneExp), geneExp)
write.table(out, file="FerrGeneExp.txt", sep="\t", quote=F, col.names=F)


######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056