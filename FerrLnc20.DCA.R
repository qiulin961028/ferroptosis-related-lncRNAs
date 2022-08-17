######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056

#install.packages("survival")
#install.packages("ggDCA")


#���ð�
library(survival)
library(ggDCA)
riskFile="risk.txt"         #���������ļ�
cliFile="clinical.txt"      #�ٴ������ļ�
setwd("D:\\biowolf\\FerrLnc\\20.DCA")     #�޸Ĺ���Ŀ¼

#��ȡ���������ļ�
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk=risk[,c("futime", "fustat", "risk")]

#��ȡ�ٴ������ļ�
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

#�ϲ�����
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)
rt[,"Age"]=ifelse(rt[,"Age"]>65, 1, 0)

#DCA����
predictTime=1    #Ԥ��ʱ��
Risk<-coxph(Surv(futime,fustat)~risk,rt)
Age<-coxph(Surv(futime,fustat)~Age,rt)
Gender<-coxph(Surv(futime,fustat)~Gender,rt)
Grade<-coxph(Surv(futime,fustat)~Grade,rt)
Stage<-coxph(Surv(futime,fustat)~Stage,rt)

#���ƾ�������
pdf(file="DCA.pdf", width=6.5, height=5.2)
d_train=dca(Risk,Age,Gender,Grade,Stage, times=predictTime)
ggplot(d_train, linetype=1)
dev.off()


######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056