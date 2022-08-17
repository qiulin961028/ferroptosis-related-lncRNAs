######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library(limma)
corFilter=0.4            #相关系数过滤标准
pvalueFilter=0.001       #p值过滤标准
setwd("D:\\biowolf\\FerrLnc\\09.FerrLncExp")     #设置工作目录

#读取lncRNA表达文件,并对数据进行处理
rt = read.table("lncRNA.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.2,]

#删掉正常样品
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
lncRNA=data[,group==0]

#读取铁死亡基因表达文件,并对数据进行处理
rt = read.table("FerrGeneExp.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
ferrGene=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
ferrGene=avereps(ferrGene)
ferrGene=ferrGene[rowMeans(ferrGene)>0.2,]

#删掉正常样品
group=sapply(strsplit(colnames(ferrGene),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
ferrGene=ferrGene[,group==0]

#相关性检验
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

#输出相关性结果
write.table(file="corResult.txt",outTab,sep="\t",quote=F,row.names=F)

#输出铁死亡lncRNA的表达量
ferrLncRNA=unique(as.vector(outTab[,"lncRNA"]))
ferrLncRNAexp=data[ferrLncRNA,]
ferrLncRNAexp=rbind(ID=colnames(ferrLncRNAexp), ferrLncRNAexp)
write.table(ferrLncRNAexp,file="FerrLncExp.txt",sep="\t",quote=F,col.names=F)


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056
