library(data.table)

pheno=fread("../../GWAS/phenotype_simul_noedu.txt")
bmi.prs=fread("BMI_all_prs.tsv")
phen=merge(pheno,bmi.prs,by="iid")
bmi.prs=fread("BMI_male_prs.tsv")
phen=merge(phen,bmi.prs,by="iid")
bmi.prs=fread("BMI_female_prs.tsv")
phen=merge(phen,bmi.prs,by="iid")

phen=phen[phen$ethnicity=="g_bri" & phen$unrelated==TRUE,]


#### Baseline results

### set baseline cov
std.cov=paste0(paste0("pc",1:40,collapse="+"),"+batch+array")


general.t2d.phen=names(phen)[grep("t2d_bias",names(phen))]
general.t2d.all=general.t2d.phen[-grep("_m|_f",general.t2d.phen)]
general.t2d.f=general.t2d.phen[grep("_f",general.t2d.phen)]
general.t2d.m=general.t2d.phen[grep("_m",general.t2d.phen)]


####### Both sexes
unbiased=matrix(nrow=3,ncol=4)

phen$bmi=phen$body_mass_index_bmi
phen$bmi[phen$bmi<15]=NA
phen$bmi=(phen$bmi-mean(phen$bmi,na.rm=T))/sd(phen$bmi,na.rm=T)

#phen$bmi[phen$bmi>4]=NA

step1=summary(lm(as.formula(paste("bmi~bmi_all+age+sex+I(age^2)+",std.cov)),data=phen))
step2=try(summary(glm(as.formula(paste("fc1_poss_t2dm~bmi_all+age+sex+I(age^2)+",std.cov)),family=binomial(),data=phen)))
beta.zy=step2$coefficients["bmi_all","Estimate"]
beta.zx=step1$coefficients["bmi_all","Estimate"]
se.zy=step2$coefficients["bmi_all","Std. Error"]
se.zx=step1$coefficients["bmi_all","Std. Error"]
cov.zx.zy=0
beta=beta.zy/beta.zx
se=((se.zy^2) / beta.zx^2) + (((beta.zy^2)/(beta.zx^4))*(se.zx^2)) - (2*(beta.zy/(beta.zx^3)) *cov.zx.zy)
se=sqrt(se)
unbiased[1,]=c("Unbiased",beta,se,"Combined")


idx.samp=which(phen[,sex==1])
step1=summary(lm(as.formula(paste("bmi~bmi_male+age+sex+I(age^2)+",std.cov)),data=phen[idx.samp,]))
step2=try(summary(glm(as.formula(paste("fc1_poss_t2dm~bmi_male+age+sex+I(age^2)+",std.cov)),family=binomial(),data=phen[idx.samp,])))
beta.zy=step2$coefficients["bmi_male","Estimate"]
beta.zx=step1$coefficients["bmi_male","Estimate"]
se.zy=step2$coefficients["bmi_male","Std. Error"]
se.zx=step1$coefficients["bmi_male","Std. Error"]
cov.zx.zy=0
beta=beta.zy/beta.zx
se=((se.zy^2) / beta.zx^2) + (((beta.zy^2)/(beta.zx^4))*(se.zx^2)) - (2*(beta.zy/(beta.zx^3)) *cov.zx.zy)
se=sqrt(se)
unbiased[2,]=c("Unbiased",beta,se,"Men")

idx.samp=which(phen[,sex==0])
step1=summary(lm(as.formula(paste("bmi~bmi_female+age+sex+I(age^2)+",std.cov)),data=phen[idx.samp,]))
step2=try(summary(glm(as.formula(paste("fc1_poss_t2dm~bmi_female+age+sex+I(age^2)+",std.cov)),family=binomial(),data=phen[idx.samp,])))
beta.zy=step2$coefficients["bmi_female","Estimate"]
beta.zx=step1$coefficients["bmi_female","Estimate"]
se.zy=step2$coefficients["bmi_female","Std. Error"]
se.zx=step1$coefficients["bmi_female","Std. Error"]
cov.zx.zy=0
beta=beta.zy/beta.zx
se=((se.zy^2) / beta.zx^2) + (((beta.zy^2)/(beta.zx^4))*(se.zx^2)) - (2*(beta.zy/(beta.zx^3)) *cov.zx.zy)
se=sqrt(se)
unbiased[3,]=c("Unbiased",beta,se,"Women")

colnames(unbiased)=c("BMI_OR","beta","se","Type")
unbiased=data.frame(unbiased)
unbiased$beta=as.numeric(as.character(unbiased$beta))
unbiased$se=as.numeric(as.character(unbiased$se))





parameter.values=c(exp(-log(rev(c(1,1.5,1.8,2,4)))),exp(log(c(1.5,1.8,2,4))))
set.seed(123456)
res.all=c()
res.men=c()
res.women=c()
k=1
for(i in parameter.values){
  for(j in 1:100){
  
  z.pmen=(phen$bmi*log(i))
  z.pmen[phen$sex==0]=z.pmen[phen$sex==0]*-1
  prob.men=1/(1+exp(-z.pmen))
  phen$select=rbinom(n=nrow(phen),size = 1,prob = prob.men)
  
  
  
  step1=summary(lm(as.formula(paste("bmi~bmi_all+sex+age+I(age^2)+",std.cov)),data=phen[which(phen[,select]==1),]))
  step2=try(summary(glm(as.formula(paste("fc1_poss_t2dm~bmi_all+age+sex+I(age^2)+",std.cov)),family=binomial(),data=phen[which(phen[,select]==1),])))
  beta.zy=step2$coefficients["bmi_all","Estimate"]
  beta.zx=step1$coefficients["bmi_all","Estimate"]
  se.zy=step2$coefficients["bmi_all","Std. Error"]
  se.zx=step1$coefficients["bmi_all","Std. Error"]
  cov.zx.zy=0
  beta=beta.zy/beta.zx
  se=((se.zy^2) / beta.zx^2) + (((beta.zy^2)/(beta.zx^4))*(se.zx^2)) - (2*(beta.zy/(beta.zx^3)) *cov.zx.zy)
  se=sqrt(se)
  res.all=rbind(res.all,c(i,beta,se))
  
  step1=summary(lm(as.formula(paste("bmi~bmi_male+age+I(age^2)+",std.cov)),data=phen[which(phen[,select]==1 & phen[,sex]==1),]))
  step2=try(summary(glm(as.formula(paste("fc1_poss_t2dm~bmi_male+age+I(age^2)+",std.cov)),family=binomial(),data=phen[which(phen[,select]==1 & phen[,sex]==1),])))
  beta.zy=step2$coefficients["bmi_male","Estimate"]
  beta.zx=step1$coefficients["bmi_male","Estimate"]
  se.zy=step2$coefficients["bmi_male","Std. Error"]
  se.zx=step1$coefficients["bmi_male","Std. Error"]
  cov.zx.zy=0
  beta=beta.zy/beta.zx
  se=((se.zy^2) / beta.zx^2) + (((beta.zy^2)/(beta.zx^4))*(se.zx^2)) - (2*(beta.zy/(beta.zx^3)) *cov.zx.zy)
  se=sqrt(se)
  res.men=rbind(res.men,c(i,beta,se))
  
  step1=summary(lm(as.formula(paste("bmi~bmi_female+age+I(age^2)+",std.cov)),data=phen[which(phen[,select]==1 & phen[,sex]==0),]))
  step2=try(summary(glm(as.formula(paste("fc1_poss_t2dm~bmi_female+age+I(age^2)+",std.cov)),family=binomial(),data=phen[which(phen[,select]==1 & phen[,sex]==0),])))
  beta.zy=step2$coefficients["bmi_female","Estimate"]
  beta.zx=step1$coefficients["bmi_female","Estimate"]
  se.zy=step2$coefficients["bmi_female","Std. Error"]
  se.zx=step1$coefficients["bmi_female","Std. Error"]
  cov.zx.zy=0
  beta=beta.zy/beta.zx
  se=((se.zy^2) / beta.zx^2) + (((beta.zy^2)/(beta.zx^4))*(se.zx^2)) - (2*(beta.zy/(beta.zx^3)) *cov.zx.zy)
  se=sqrt(se)
  res.women=rbind(res.women,c(i,beta,se))
  
  print(j)
  }

}

#### Summarise results



combined=cbind(unique(res.all[,1]),by(res.all[,2],res.all[,1],mean),by(res.all[,2],res.all[,1],sd),"Combined")
men=cbind(unique(res.men[,1]),by(res.men[,2],res.men[,1],mean),by(res.men[,2],res.men[,1],sd),"Men")
women=cbind(unique(res.women[,1]),by(res.women[,2],res.women[,1],mean),by(res.women[,2],res.all[,1],sd),"Women")

results=rbind(combined,men,women)

results=as.data.frame(results)
names(results)=c("BMI_OR","beta","se","Type")
results$beta=as.numeric(as.character(results$beta))
results$se=as.numeric(as.character(results$se))

results$high=results$beta+qnorm(p = 1-0.025,lower.tail = T)*results$se
results$low=results$beta+qnorm(p = 1-0.025,lower.tail = F)*results$se





library(ggplot2)
unbiased$high=unbiased$beta+qnorm(p = 1-0.025,lower.tail = T)*unbiased$se
unbiased$low=unbiased$beta+qnorm(p = 1-0.025,lower.tail = F)*unbiased$se

pdf("original_forestplot.pdf")
ggplot(unbiased,aes(x=Type,y=exp(beta),ymin=exp(low),ymax=exp(high)))+geom_pointrange()+theme_minimal()+ylim(0,4)+geom_hline(yintercept=1, linetype="dashed", color="red")
dev.off()
results$BMI_OR=as.factor(round(as.numeric(as.character(results$BMI_OR)),2))
  
massimi=by(exp(results$high),INDICES = results$BMI_OR,FUN = max)

df1=data.frame(x1 =c(1,1.2),y1 =c(massimi[1]+0.1,massimi[1]+0.1))  ## x=y
df2=data.frame(x1 =c(2,2.2),y1 =c(massimi[2]+0.1,massimi[2]+0.1))  ## x=y
df3=data.frame(x1 =c(3,3.2),y1 =c(massimi[3]+0.1,massimi[3]+0.1))  ## x=y
df4=data.frame(x1 =c(4,4.2),y1 =c(massimi[4]+0.1,massimi[4]+0.1))  ## x=y
df5=data.frame(x1 =c(5,5.2),y1 =c(massimi[5]+0.1,massimi[5]+0.1))  ## x=y
df6=data.frame(x1 =c(6,6.2),y1 =c(massimi[6]+0.1,massimi[6]+0.1))  ## x=y
df7=data.frame(x1 =c(7,7.2),y1 =c(massimi[7]+0.1,massimi[7]+0.1))  ## x=y
df8=data.frame(x1 =c(8,8.2),y1 =c(massimi[8]+0.1,massimi[8]+0.1))  ## x=y
df9=data.frame(x1 =c(9,9.2),y1 =c(massimi[9]+0.1,massimi[9]+0.1))  ## x=y

pvals=pchisq(((as.numeric(men[,2])-as.numeric(women[,2]))/sqrt(as.numeric(men[,3])^2+as.numeric(women[,3])^2))^2,df=1,lower=F)
pvals=formatC(pvals, format = "e", digits = 1)

pdf("Results_simulations.pdf",width=21,height=7)
p1=ggplot(results,aes(x=BMI_OR,y=exp(beta),ymin=exp(low),ymax=exp(high),color=as.factor(Type)))+geom_pointrange(position=position_dodge(0.5),size=1.5)+
  theme_minimal()+ylim(0,6)+geom_hline(yintercept=1, linetype="dashed", color="red")+
  geom_hline(yintercept=exp(unbiased$beta[1]), linetype="solid", color="#F8766D")+geom_hline(yintercept=exp(unbiased$beta[2]), linetype="solid", color="#00BA38")+
  geom_hline(yintercept=exp(unbiased$beta[3]), linetype="solid", color="#619CFF")+theme(axis.text.x = element_text(size=18),
                                                                                        axis.text.y = element_text(size=18),
                                                                                        axis.title.y =element_text(size=20),
                                                                                        axis.title.x =element_text(size=20),
                                                                                        legend.title = element_blank(),
                                                                                        legend.text = element_text(size=20))+
  
  xlab("Participation Bias")+ylab("Causal Effect \n BMI->T2D (OR x SD)")

  p1+geom_line(aes(x =x1, y =y1 ),data=df1,inherit.aes = FALSE) + annotate("text", x = mean(df1$x1), y = (massimi[1]+0.2), label = pvals[1], size = 3)+
    geom_line(aes(x =x1, y =y1 ),data=df2,inherit.aes = FALSE) + annotate("text", x = mean(df2$x1), y = (massimi[2]+0.2), label = pvals[2], size = 3)+
    geom_line(aes(x =x1, y =y1 ),data=df3,inherit.aes = FALSE) + annotate("text", x = mean(df3$x1), y = (massimi[3]+0.2), label = pvals[3], size = 3)+
    geom_line(aes(x =x1, y =y1 ),data=df4,inherit.aes = FALSE) + annotate("text", x = mean(df4$x1), y = (massimi[4]+0.2), label = pvals[4], size = 3)+
    geom_line(aes(x =x1, y =y1 ),data=df5,inherit.aes = FALSE) + annotate("text", x = mean(df5$x1), y = (massimi[5]+0.2), label = pvals[5], size = 3)+
    geom_line(aes(x =x1, y =y1 ),data=df6,inherit.aes = FALSE) + annotate("text", x = mean(df6$x1), y = (massimi[6]+0.2), label = pvals[6], size = 3)+
    geom_line(aes(x =x1, y =y1 ),data=df7,inherit.aes = FALSE) + annotate("text", x = mean(df7$x1), y = (massimi[7]+0.2), label = pvals[7], size = 3)+
    geom_line(aes(x =x1, y =y1 ),data=df8,inherit.aes = FALSE) + annotate("text", x = mean(df8$x1), y = (massimi[8]+0.2), label = pvals[8], size = 3)+
    geom_line(aes(x =x1, y =y1 ),data=df9,inherit.aes = FALSE) + annotate("text", x = mean(df9$x1), y = (massimi[9]+0.2), label = pvals[9], size = 3)
dev.off()

  













