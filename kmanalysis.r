library(dplyr)
library(survival)
library(survminer)

head(liver_km_R[, 1:4])
res.cut <- surv_cutpoint(liver_km_R, time = "months", event = "vital",
                         variables = c("FRG1"))
summary(res.cut)
plot(res.cut, "FRG1", palette = "npg",font.x=14, font.y=14,
     font.tickslab=14, font.legend=15)
res.cat <- surv_categorize(res.cut)
head(res.cat)
library("survival")
fit <- survfit(Surv(months, vital) ~FRG1, data = res.cat)

ggsurvplot(fit, risk.table = TRUE,  pval=TRUE, font.x=18, 
           font.y=18,font.tickslab=18, font.legend=15,palette=c("red", "blue") ,
           xlab = "Time (Months)")
dev.off()



 ggsurvplot(fit, conf.int=TRUE, pval=TRUE, risk.table=TRUE, 
           legend.title="Sex",  
           palette=c("dodgerblue2", "orchid2"),
           title="Kaplan-Meier Curve for Lung Cancer Survival", 
           risk.table.height=.15)
View(gastric_km_R)
# conf.int = TRUE palette=c("red", "blue")
pdf("pract12.pdf", width = 6,height = 6)
jpeg("rplot2.jpg", width = 3500, height = 3500)
