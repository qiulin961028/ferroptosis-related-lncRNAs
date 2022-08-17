######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("survival")
#install.packages("ggDCA")


#引用包
library(survival)
library(ggDCA)
riskFile="risk.txt"         #风险输入文件
cliFile="clinical.txt"      #临床数据文件
setwd("D:\\biowolf\\FerrLnc\\20.DCA")     #修改工作目录

#读取风险输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk=risk[,c("futime", "fustat", "risk")]

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)
rt[,"Age"]=ifelse(rt[,"Age"]>65, 1, 0)

#DCA分析
predictTime=1    #预测时间
Risk<-coxph(Surv(futime,fustat)~risk,rt)
Age<-coxph(Surv(futime,fustat)~Age,rt)
Gender<-coxph(Surv(futime,fustat)~Gender,rt)
Grade<-coxph(Surv(futime,fustat)~Grade,rt)
Stage<-coxph(Surv(futime,fustat)~Stage,rt)

#绘制决策曲线
pdf(file="DCA.pdf", width=6.5, height=5.2)
d_train=dca(Risk,Age,Gender,Grade,Stage, times=predictTime)
ggplot(d_train, linetype=1)
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056
