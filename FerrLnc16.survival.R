######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("survival")
#install.packages("survminer")


#引用包
library(survival)
library(survminer)
setwd("D:\\biowolf\\FerrLnc\\16.survival")     #设置工作目录

#定义生存曲线的函数
bioSurvival=function(inputFile=null, outFile=null){
	#读取输入文件
	rt=read.table(inputFile, header=T, sep="\t", check.names=F)
	#比较高低风险组生存差异，得到显著性p值
	diff=survdiff(Surv(futime, fustat) ~ risk, data=rt)
	pValue=1-pchisq(diff$chisq, df=1)
	if(pValue<0.001){
		pValue="p<0.001"
	}else{
		pValue=paste0("p=",sprintf("%.03f",pValue))
	}
	fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
	#print(surv_median(fit))
		
	#绘制生存曲线
	surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=T,
		           pval=pValue,
		           pval.size=6,
		           surv.median.line = "hv",
		           legend.title="Risk",
		           legend.labs=c("High risk", "Low risk"),
		           xlab="Time(years)",
		           break.time.by = 1,
		           palette=c("red", "blue"),
		           risk.table=TRUE,
		       	   risk.table.title="",
		           risk.table.col = "strata",
		           risk.table.height=.25)
	pdf(file=outFile, onefile=FALSE, width=6.5, height=5.5)
	print(surPlot)
	dev.off()
}

#调用函数，绘制生存曲线
bioSurvival(inputFile="risk.txt", outFile="survival.pdf")


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056
