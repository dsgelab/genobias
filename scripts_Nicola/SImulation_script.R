library(data.table)

pheno=fread("../../GWAS/phenotype_simul_noedu.txt")
bmi.prs=fread("BMI_all_prs.tsv")
phen=merge(pheno,bmi.prs,by="iid")
bmi.prs=fread("BMI_male_prs.tsv")
phen=merge(phen,bmi.prs,by="iid")
bmi.prs=fread("BMI_female_prs.tsv")
phen=merge(phen,bmi.prs,by="iid")



#### Baseline results

### set baseline cov
std.cov=paste0(paste0("pc",1:40,collapse="+"),"+g_bri+unrelated+batch+array")


general.t2d.phen=names(phen)[grep("t2d_bias",names(phen))]
general.t2d.all=general.t2d.phen[-grep("_m|_f",general.t2d.phen)]
general.t2d.f=general.t2d.phen[grep("_f",general.t2d.phen)]
general.t2d.m=general.t2d.phen[grep("_m",general.t2d.phen)]


####### Both sexes
res.tot=matrix(nrow=length(general.t2d.all)+1,ncol=7)

k=1
phen$all=1
a="all"
i="t2d_or"

idx.samp=which(phen[,..a]==1)

step1=summary(lm(as.formula(paste("bmi~bmi_all+age+sex+I(age^2)+",std.cov)),data=phen[idx.samp,]))
step2=try(summary(glm(as.formula(paste(i,"~bmi_all+age+sex+I(age^2)+",std.cov)),family=binomial(),data=phen[idx.samp,])))
beta.zy=step2$coefficients["bmi_all","Estimate"]
beta.zx=step1$coefficients["bmi_all","Estimate"]
se.zy=step2$coefficients["bmi_all","Std. Error"]
se.zx=step1$coefficients["bmi_all","Std. Error"]
cov.zx.zy=0
beta=beta.zy/beta.zx
se=((se.zy^2) / beta.zx^2) + (((beta.zy^2)/(beta.zx^4))*(se.zx^2)) - (2*(beta.zy/(beta.zx^3)) *cov.zx.zy)
se=sqrt(se)
c.c=table(model.frame(as.formula(paste(i,"~bmi_all+age+sex+I(age^2)+",std.cov)),data=phen[which(phen[,..a]==1),])[,c(i,"sex")])

res.tot[k,]=c(i,beta,se,c.c[,1],c.c[,2])
k=k+1


parameter.values=c(exp(-log(rev(c(1,1.5,1.8,2,4)))),exp(log(c(1.5,1.8,2,4))))
set.seed(123456)
res.all=c()
res.man=c()
res.women=c()
k=1
for(i in parameter.values){
  for(j in 1:100){
  
  z.pmen=(phen$bmi*log(i))
  z.pmen[phen$sex==0]=z.pmen[phen$sex==0]*-1
  prob.men=1/(1+exp(-z.pmen))
  phen$select=rbinom(n=nrow(phen),size = 1,prob = prob.men)
  
  
  step1=summary(lm(as.formula(paste("bmi~bmi_all+sex+age+I(age^2)+",std.cov)),data=phen[which(phen[,select]==1),]))
  step2=try(summary(glm(as.formula(paste("t2d_or~bmi_all+age+sex+I(age^2)+",std.cov)),family=binomial(),data=phen[which(phen[,select]==1),])))
  beta.zy=step2$coefficients["bmi_all","Estimate"]
  beta.zx=step1$coefficients["bmi_all","Estimate"]
  se.zy=step2$coefficients["bmi_all","Std. Error"]
  se.zx=step1$coefficients["bmi_all","Std. Error"]
  cov.zx.zy=0
  beta=beta.zy/beta.zx
  se=((se.zy^2) / beta.zx^2) + (((beta.zy^2)/(beta.zx^4))*(se.zx^2)) - (2*(beta.zy/(beta.zx^3)) *cov.zx.zy)
  se=sqrt(se)
  res.all=rbind(res.all,c(i,beta,se))
  
  step1=summary(lm(as.formula(paste("bmi~bmi_all+age+I(age^2)+",std.cov)),data=phen[which(phen[,select]==1 & phen[,sex]==1),]))
  step2=try(summary(glm(as.formula(paste("t2d_or~bmi_all+age+I(age^2)+",std.cov)),family=binomial(),data=phen[which(phen[,select]==1 & phen[,sex]==1),])))
  beta.zy=step2$coefficients["bmi_all","Estimate"]
  beta.zx=step1$coefficients["bmi_all","Estimate"]
  se.zy=step2$coefficients["bmi_all","Std. Error"]
  se.zx=step1$coefficients["bmi_all","Std. Error"]
  cov.zx.zy=0
  beta=beta.zy/beta.zx
  se=((se.zy^2) / beta.zx^2) + (((beta.zy^2)/(beta.zx^4))*(se.zx^2)) - (2*(beta.zy/(beta.zx^3)) *cov.zx.zy)
  se=sqrt(se)
  res.men=rbind(res.men,c(i,beta,se))
  
  step1=summary(lm(as.formula(paste("bmi~bmi_all+age+I(age^2)+",std.cov)),data=phen[which(phen[,select]==1 & phen[,sex]==0),]))
  step2=try(summary(glm(as.formula(paste("t2d_or~bmi_all+age+I(age^2)+",std.cov)),family=binomial(),data=phen[which(phen[,select]==1 & phen[,sex]==0),])))
  beta.zy=step2$coefficients["bmi_all","Estimate"]
  beta.zx=step1$coefficients["bmi_all","Estimate"]
  se.zy=step2$coefficients["bmi_all","Std. Error"]
  se.zx=step1$coefficients["bmi_all","Std. Error"]
  cov.zx.zy=0
  beta=beta.zy/beta.zx
  se=((se.zy^2) / beta.zx^2) + (((beta.zy^2)/(beta.zx^4))*(se.zx^2)) - (2*(beta.zy/(beta.zx^3)) *cov.zx.zy)
  se=sqrt(se)
  res.women=rbind(res.women,c(i,beta,se))
  
  print(j)
  }

}


colnames(res.tot)=c("Trait","beta","se","controls_F","Cases_F","Controls_M","Cases_M")
write.table(res.tot,file="BMI_T2D_alla_samples.tsv",row.names=F,quote=F,sep="\t")










