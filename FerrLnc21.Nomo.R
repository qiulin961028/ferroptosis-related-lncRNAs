######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("survival")
#install.packages("regplot")


#引用包
library(survival)
library(regplot)
riskFile="risk.txt"       #风险输入文件
cliFile="clinical.txt"    #临床数据文件
setwd("D:\\biowolf\\FerrLnc\\21.Nomo")     #修改工作目录

#读取风险输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk=risk[,c("futime", "fustat", "risk")]

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
cli$Age=as.numeric(cli$Age)

#合并数据
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)

#绘制列线图
res.cox=coxph(Surv(futime, fustat) ~ . , data = rt)
nom1<-regplot(res.cox,
              plots = c("density", "boxes"),
              clickable=F,
              title="",
              points=TRUE,
              droplines=TRUE,
              observation=rt[1,],
              rank="sd",
              failtime = c(1,3,5),
              prfail = T)


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

  