######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library(limma)
corFilter=0.4            #���ϵ�����˱�׼
pvalueFilter=0.001       #pֵ���˱�׼
setwd("D:\\biowolf\\FerrLnc\\09.FerrLncExp")     #���ù���Ŀ¼

#��ȡlncRNA�����ļ�,�������ݽ��д���
rt = read.table("lncRNA.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.2,]

#ɾ��������Ʒ
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
lncRNA=data[,group==0]

#��ȡ��������������ļ�,�������ݽ��д���
rt = read.table("FerrGeneExp.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
ferrGene=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
ferrGene=avereps(ferrGene)
ferrGene=ferrGene[rowMeans(ferrGene)>0.2,]

#ɾ��������Ʒ
group=sapply(strsplit(colnames(ferrGene),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
ferrGene=ferrGene[,group==0]

#����Լ���
outTab=data.frame()
for(i in row.names(lncRNA)){
	if(sd(lncRNA[i,])>0.2){
		for(j in row.names(ferrGene)){
			x=as.numeric(lncRNA[i,])
			y=as.numeric(ferrGene[j,])
			corT=cor.test(x,y)
			cor=corT$estimate
			pvalue=corT$p.value
			if((cor>corFilter) & (pvalue<pvalueFilter)){
				outTab=rbind(outTab,cbind(ferrGene=j,lncRNA=i,cor,pvalue,Regulation="postive"))
			}
			if((cor< -corFilter) & (pvalue<pvalueFilter)){
				outTab=rbind(outTab,cbind(ferrGene=j,lncRNA=i,cor,pvalue,Regulation="negative"))
			}
		}
	}
}

#�������Խ��
write.table(file="corResult.txt",outTab,sep="\t",quote=F,row.names=F)

#���������lncRNA�ı�����
ferrLncRNA=unique(as.vector(outTab[,"lncRNA"]))
ferrLncRNAexp=data[ferrLncRNA,]
ferrLncRNAexp=rbind(ID=colnames(ferrLncRNAexp), ferrLncRNAexp)
write.table(ferrLncRNAexp,file="FerrLncExp.txt",sep="\t",quote=F,col.names=F)


######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056