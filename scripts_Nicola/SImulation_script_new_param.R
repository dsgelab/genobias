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

write.table(unbiased,file="unbiased_T2D_BMI_MR_gen_briti.tsv",row.names=F,quote=F,sep="\t")



parameter.values=c(exp(-log(rev(c(1,1.2,1.5,1.8,2,3)))),exp(log(c(1.2,1.5,1.8,2,3))))
set.seed(123456)
res.all=c()
res.men=c()
res.women=c()
ks=c(-0.5, -0.3, 0, 0.3, 0.7, 1, 1.5)
for(i in parameter.values){
  for(k in ks){
  for(j in 1:100){
  
  z.pmen=(phen$bmi*log(i))+rnorm(n = length(phen$bmi),mean = 0,sd = 0.1)
  z.pmen[phen$sex==0]=z.pmen[phen$sex==0]*k
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
  res.all=rbind(res.all,c(i,k,beta,se))
  
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
  res.men=rbind(res.men,c(i,k,beta,se))
  
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
  res.women=rbind(res.women,c(i,k,beta,se))
  
  print(j)
  }

 }
}
#### Summarise results

write.table(res.all,file="Full_results_simulations_combined_new_param.tsv",row.names=F,quote=F,sep="\t")
write.table(res.men,file="Full_results_simulations_men_new_param.tsv",row.names=F,quote=F,sep="\t")
write.table(res.women,file="Full_results_simulations_women_new_param.tsv",row.names=F,quote=F,sep="\t")

res.all=fread("Full_results_simulations_combined_new_param.tsv")
res.men=fread("Full_results_simulations_men_new_param.tsv")
res.women=fread("Full_results_simulations_women_new_param.tsv")
res.all=data.frame(res.all)
res.men=data.frame(res.men)
res.women=data.frame(res.women)

names(res.all)=c("BMI_OR","k","beta","se")
names(res.men)=c("BMI_OR","k","beta","se")
names(res.women)=c("BMI_OR","k","beta","se")

res.all$type="Combined"
res.men$type="Men"
res.women$type="Women"

res.all.beta=merge(aggregate(beta~BMI_OR+k,res.all,FUN = mean),aggregate(beta~BMI_OR+k,res.all,FUN = sd),by=c("BMI_OR","k"))
names(res.all.beta)[3:4]=c("beta","se")
res.all.beta$Type="Combined"

res.men.beta=merge(aggregate(beta~BMI_OR+k,res.men,FUN = mean),aggregate(beta~BMI_OR+k,res.men,FUN = sd),by=c("BMI_OR","k"))
names(res.men.beta)[3:4]=c("beta","se")
res.men.beta$Type="Men"

res.women.beta=merge(aggregate(beta~BMI_OR+k,res.women,FUN = mean),aggregate(beta~BMI_OR+k,res.women,FUN = sd),by=c("BMI_OR","k"))
names(res.women.beta)[3:4]=c("beta","se")
res.women.beta$Type="Women"



results=rbind(res.all.beta,res.men.beta,res.women.beta)

results$beta=as.numeric(as.character(results$beta))
results$se=as.numeric(as.character(results$se))

results$high=results$beta+qnorm(p = 1-0.025,lower.tail = T)*results$se
results$low=results$beta+qnorm(p = 1-0.025,lower.tail = F)*results$se
results$BMI_OR=as.factor(round(as.numeric(as.character(results$BMI_OR)),2))

write.table(results,file="Results_simulations_k.tsv",row.names=F,quote=F,sep="\t")



