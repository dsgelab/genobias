### Create rg subset

file.list.bgen=system("ls *.bgen",intern=T)

file.list.bgi=system("ls *.bgi",intern=T)


for( i in file.list.bgen){
  
  if(!(paste0(i,".bgi")%in%file.list.bgi)) system(paste("bgenix -g",i,"-index"))
    
}


for(i in 1:241){
  
  
   system(paste0("bgenix -g /exports/igmm/eddie/wilson-lab/data/base_data/ukbb_19655/genotypes/hrc/bgen/chunks/ukbb_19655_chunk_",i,".bgen -incl-rsids SNP_4_MR.tsv > ukbb_19655_MR_chunk_",i,".bgen"))
  
}

library(data.table)
library(bigsnpr)
#system("/exports/igmm/eddie/wilson-lab/app_full_stuff/bgen/gavinband-bgen-4e33223a8dc4/build/apps/cat-bgen -g ukbb_19655_MR_chunk_*.bgen  -og MR_genotypes_BMI.bgen")
#system("bgenix -index MR_genotypes_BMI.bgen")
info=fread("info.tsv")
info$alternate_ids=paste0(info$chromosome,"_",info$position,"_",info$first_allele,"_",info$alternative_alleles)
snp.names=list()
snp.names[[1]]=info$alternate_ids
rds=snp_readBGEN(bgenfiles="MR_genotypes_BMI.bgen"
                   , backingfile="./rds"
                   , ind_row = NULL
                 ,list_snp_id=snp.names
                 , read_as = "dosage"
                 ,ncores = 1
                 )

obj.bigSNP=snp_attach("/gpfs/igmmfs01/eddie/wilson-lab/projects/prj_109_sex_gwas/MR_simul/ldscore/rds.rds")
dos=obj.bigSNP$genotypes[,1:ncol(obj.bigSNP$genotypes)]

sample=fread( "/exports/igmm/eddie/wilson-lab/data/base_data/ukbb_19655/genotypes/hrc/bgen/ukbb_19655.sample")
# Create all PRS
bmi.prs=fread("../bmi_all_intrument.txt")
bmi.prs$SNP=apply(t(bmi.prs$SNP),2,function(x)unlist(strsplit(x,split="_"))[1])
bmi.prs=bmi.prs[,c("SNP","Chr","Pos","NEA","EA","beta.combined")]
names(bmi.prs)=c("rsid","chr", "pos", "a0", "a1" , "beta")

info.prs=info[,c("rsid","chromosome", "position","first_allele","alternative_alleles")]
names(info.prs)=c("rsid","chr", "pos", "a0", "a1" )


flipp=snp_match(sumstats=bmi.prs, info_snp=info.prs, strand_flip = FALSE, join_by_pos = FALSE,
          match.min.prop = 0.5)

prs=snp_PRS(obj.bigSNP$genotypes, betas.keep=flipp$beta)
prs.all=data.frame(iid=sample$ID_1[-1],bmi_all=prs,stringsAsFactors = F)
write.table(prs.all,file="BMI_all_prs.tsv",row.names=F,quote=F,sep="\t")


### BMI men
bmi.prs.male=fread("../bmi_men_intrument.txt")
bmi.prs.male$SNP=apply(t(bmi.prs.male$SNP),2,function(x)unlist(strsplit(x,split="_"))[1])
bmi.prs.male=bmi.prs.male[,c("SNP","Chr","Pos","NEA","EA","beta.males")]
names(bmi.prs.male)=c("rsid","chr", "pos", "a0", "a1" , "beta")

flipp=snp_match(sumstats=bmi.prs.male, info_snp=info.prs, strand_flip = FALSE, join_by_pos = FALSE,
                match.min.prop = 0.5)
which(flipp$rsid!=info$rsid)
prs=snp_PRS(obj.bigSNP$genotypes, betas.keep=flipp$beta)
prs.male=data.frame(iid=sample$ID_1[-1],bmi_male=prs,stringsAsFactors = F)
write.table(prs.male,file="BMI_male_prs.tsv",row.names=F,quote=F,sep="\t")

### PRS female

bmi.prs.female=fread("../bmi_women_intrument.txt")
bmi.prs.female$SNP=apply(t(bmi.prs.female$SNP),2,function(x)unlist(strsplit(x,split="_"))[1])
bmi.prs.female=bmi.prs.female[,c("SNP","Chr","Pos","NEA","EA","beta.females")]
names(bmi.prs.female)=c("rsid","chr", "pos", "a0", "a1" , "beta")

flipp=snp_match(sumstats=bmi.prs.female, info_snp=info.prs, strand_flip = FALSE, join_by_pos = FALSE,
                match.min.prop = 0.5)
which(flipp$rsid!=info$rsid)
prs=snp_PRS(obj.bigSNP$genotypes, betas.keep=flipp$beta)
prs.female=data.frame(iid=sample$ID_1[-1],bmi_female=prs,stringsAsFactors = F)
write.table(prs.female,file="BMI_female_prs.tsv",row.names=F,quote=F,sep="\t")


######## MR

library(data.table)

pheno=fread("../../GWAS/phenotype_simul.txt")
bmi.prs=fread("BMI_all_prs.tsv")
phen=merge(pheno,bmi.prs,by="iid")
bmi.prs=fread("BMI_male_prs.tsv")
phen=merge(phen,bmi.prs,by="iid")
bmi.prs=fread("BMI_female_prs.tsv")
phen=merge(phen,bmi.prs,by="iid")


std.cov=paste0(paste0("pc",1:40,collapse="+"),"+g_bri+unrelated+batch+array")


library(estimatr)
library(speedglm)
### BMI all-> T2d

Briti+g_bri+Irish+White+Any o

general.t2d.phen=names(phen)[grep("t2d_bias",names(phen))]
#general.t2d.phen=general.t2d.phen[grep("_edu_1_",general.t2d.phen)]

general.t2d.all=general.t2d.phen[-grep("_m|_f",general.t2d.phen)]
general.t2d.f=general.t2d.phen[grep("_f",general.t2d.phen)]
general.t2d.m=general.t2d.phen[grep("_m",general.t2d.phen)]


####### Both sexes
res.tot=matrix(nrow=length(general.t2d.all)+1,ncol=7)

k=1
phen$all=1
a="all"
i="t2d_or"
step1=summary(lm(as.formula(paste("bmi~bmi_all+age+sex+I(age^2)+",std.cov)),data=phen[which(phen[,..a]==1),]))
step2=try(summary(glm(as.formula(paste(i,"~bmi_all+age+sex+I(age^2)+",std.cov)),family=binomial(),data=phen[which(phen[,..a]==1),])))
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

for(i in general.t2d.all){
  
    a=gsub("t2d_bias","select",i)
    step1=summary(lm(as.formula(paste("bmi~bmi_all+sex+age+I(age^2)+",std.cov)),data=phen[which(phen[,..a]==1),]))
    step2=try(summary(glm(as.formula(paste(i,"~bmi_all+age+sex+I(age^2)+",std.cov)),family=binomial(),data=phen[which(phen[,..a]==1),])))
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
    print(res.tot[k,])
    k=k+1
}
colnames(res.tot)=c("Trait","beta","se","controls_F","Cases_F","Controls_M","Cases_M")
write.table(res.tot,file="BMI_T2D_alla_samples.tsv",row.names=F,quote=F,sep="\t")

###### Men

res.tot.m=matrix(nrow=length(general.t2d.m)+1,ncol=7)

k=1
phen$all=1
a="all"
i="t2d_or_m"

tmp.d=phen[which(phen[,..a]==1 & phen$sex==1),]
step1=summary(lm(as.formula(paste("bmi~bmi_male+age+I(age^2)+",std.cov)),data=tmp.d))
step2=try(summary(glm(as.formula(paste(i,"~bmi_male+age+I(age^2)+",std.cov)),family=binomial(),data=tmp.d)))
beta.zy=step2$coefficients["bmi_male","Estimate"]
beta.zx=step1$coefficients["bmi_male","Estimate"]
se.zy=step2$coefficients["bmi_male","Std. Error"]
se.zx=step1$coefficients["bmi_male","Std. Error"]
cov.zx.zy=0
beta=beta.zy/beta.zx
se=((se.zy^2) / beta.zx^2) + (((beta.zy^2)/(beta.zx^4))*(se.zx^2)) - (2*(beta.zy/(beta.zx^3)) *cov.zx.zy)
se=sqrt(se)
c.c=table(model.frame(as.formula(paste(i,"~bmi_male+age+I(age^2)+",std.cov)),data=tmp.d)[,i])

res.tot.m[k,]=c(i,beta,se,NA,NA,c.c)
k=k+1

for(i in general.t2d.m){
  a=gsub("t2d_bias_m","select",i)
  tmp.d=phen[which(phen[,..a]==1 & phen$sex==1),]
  step1=summary(lm(as.formula(paste("bmi~bmi_male+age+I(age^2)+",std.cov)),data=tmp.d))
  step2=try(summary(glm(as.formula(paste(i,"~bmi_male+age+I(age^2)+",std.cov)),family=binomial(),data=tmp.d)))
  beta.zy=step2$coefficients["bmi_male","Estimate"]
  beta.zx=step1$coefficients["bmi_male","Estimate"]
  se.zy=step2$coefficients["bmi_male","Std. Error"]
  se.zx=step1$coefficients["bmi_male","Std. Error"]
  cov.zx.zy=0
  beta=beta.zy/beta.zx
  se=((se.zy^2) / beta.zx^2) + (((beta.zy^2)/(beta.zx^4))*(se.zx^2)) - (2*(beta.zy/(beta.zx^3)) *cov.zx.zy)
  se=sqrt(se)
  c.c=table(model.frame(as.formula(paste(i,"~bmi_male+age+I(age^2)+",std.cov)),data=tmp.d)[,i])
  
  res.tot.m[k,]=c(i,beta,se,NA,NA,c.c)
  k=k+1

}



colnames(res.tot.m)=c("Trait","beta","se","controls_F","Cases_F","Controls_M","Cases_M")
write.table(res.tot.m,file="BMI_T2D_Males_samples.tsv",row.names=F,quote=F,sep="\t")



###### Women

res.tot.f=matrix(nrow=length(general.t2d.f)+1,ncol=7)

k=1
phen$all=1
a="all"
i="t2d_or_f"

tmp.d=phen[which(phen[,..a]==1 & phen$sex==0),]
step1=summary(lm(as.formula(paste("bmi~bmi_female+",std.cov)),data=tmp.d))
step2=try(summary(glm(as.formula(paste(i,"~bmi_female+age+I(age^2)+",std.cov)),family=binomial(),data=tmp.d)))
beta.zy=step2$coefficients["bmi_female","Estimate"]
beta.zx=step1$coefficients["bmi_female","Estimate"]
se.zy=step2$coefficients["bmi_female","Std. Error"]
se.zx=step1$coefficients["bmi_female","Std. Error"]
cov.zx.zy=0
beta=beta.zy/beta.zx
se=((se.zy^2) / beta.zx^2) + (((beta.zy^2)/(beta.zx^4))*(se.zx^2)) - (2*(beta.zy/(beta.zx^3)) *cov.zx.zy)
se=sqrt(se)
c.c=table(model.frame(as.formula(paste(i,"~bmi_female+age+I(age^2)+",std.cov)),data=tmp.d)[,i])

res.tot.f[k,]=c(i,beta,se,NA,NA,c.c)
k=k+1

for(i in general.t2d.f){
  a=gsub("t2d_bias_f","select",i)
  tmp.d=phen[which(phen[,..a]==1 & phen$sex==0),]
  step1=summary(lm(as.formula(paste("bmi~bmi_female+age+I(age^2)+",std.cov)),data=tmp.d))
  step2=try(summary(glm(as.formula(paste(i,"~bmi_female+age+I(age^2)+",std.cov)),family=binomial(),data=tmp.d)))
  beta.zy=step2$coefficients["bmi_female","Estimate"]
  beta.zx=step1$coefficients["bmi_female","Estimate"]
  se.zy=step2$coefficients["bmi_female","Std. Error"]
  se.zx=step1$coefficients["bmi_female","Std. Error"]
  cov.zx.zy=0
  beta=beta.zy/beta.zx
  se=((se.zy^2) / beta.zx^2) + (((beta.zy^2)/(beta.zx^4))*(se.zx^2)) - (2*(beta.zy/(beta.zx^3)) *cov.zx.zy)
  se=sqrt(se)
  c.c=table(model.frame(as.formula(paste(i,"~bmi_female+age+I(age^2)+",std.cov)),data=tmp.d)[,i])
  
  res.tot.f[k,]=c(i,beta,se,c.c,NA,NA)
  k=k+1
  
}


colnames(res.tot.f)=c("Trait","beta","se","controls_F","Cases_F","Controls_M","Cases_M")
write.table(res.tot.f,file="BMI_T2D_Females_samples.tsv",row.names=F,quote=F,sep="\t")

#### gather data and plot

library(data.table)

unbiased=c()

combined=fread("BMI_T2D_alla_samples.tsv")
combined$type="combined"
unbiased=rbind(unbiased,combined[1,])
combined=combined[-1,]
men=fread("BMI_T2D_Males_samples.tsv")
men$type="men"
unbiased=rbind(unbiased,men[1,])
men=men[-1,]

women=fread("BMI_T2D_Females_samples.tsv")
women$type="women"
unbiased=rbind(unbiased,women[1,])
women=women[-1,]

library(ggplot2)
unbiased$high=unbiased$beta+qnorm(p = 1-0.025,lower.tail = T)*unbiased$se
unbiased$low=unbiased$beta+qnorm(p = 1-0.025,lower.tail = F)*unbiased$se

pdf("original_forestplot.pdf")
ggplot(unbiased,aes(x=Trait,y=exp(beta),ymin=exp(low),ymax=exp(high)))+geom_pointrange()+theme_minimal()+ylim(0,4)+geom_hline(yintercept=1, linetype="dashed", color="red")
dev.off()


all.res=rbind(combined,men)
all.res=rbind(all.res,women)

splitter=function(x)unlist(strsplit(x,split="_"))[c(1,3,5)]
param=t(apply(t(all.res$Trait),2,splitter))
all.res=cbind(all.res,param)
names(all.res)[9:11]=c("Scenario","BMI_OR","Edu_OR")
all.res$high=all.res$beta+qnorm(p = 1-0.025,lower.tail = T)*all.res$se
all.res$low=all.res$beta+qnorm(p = 1-0.025,lower.tail = F)*all.res$se


pdf("test.pdf",width=21,height=21)
ggplot(all.res[all.res$Scenario=="sc2",],aes(x=BMI_OR,y=exp(beta),ymin=exp(low),ymax=exp(high),color=as.factor(type)))+geom_pointrange(position=position_dodge(0.5))+
  theme_minimal()+ylim(0,6)+geom_hline(yintercept=1, linetype="dashed", color="red")+
  facet_grid(Edu_OR~.)+geom_hline(yintercept=3.176839, linetype="solid", color="#F8766D")+geom_hline(yintercept=3.186489, linetype="solid", color="#00BA38")+
  geom_hline(yintercept=3.273303, linetype="solid", color="#619CFF")
dev.off()


all.res

tmp=all.res[which(all.res$Edu_OR==1 & all.res$Scenario=='sc1'),]

tmp=all.res[which( all.res$Scenario=='sc2'),]

pdf("MR_simul_sc2.pdf",width=21,height=7)
ggplot(tmp[tmp$Edu_OR==1,],aes(x=BMI_OR,y=exp(beta),ymin=exp(low),ymax=exp(high),color=type))+geom_pointrange(position=position_dodge(0.5))+
  theme_minimal()+ylim(0,6)+geom_hline(yintercept=1, linetype="dashed", color="red")+
  geom_hline(yintercept=3.176839, linetype="solid", color="#F8766D")+geom_hline(yintercept=3.186489, linetype="solid", color="#00BA38")+
  geom_hline(yintercept=3.273303, linetype="solid", color="#619CFF")
dev.off()


tmp=all.res[which( all.res$Scenario=='sc1'),]
pdf("MR_simul_sc1.pdf",width=21,height=7)
ggplot(tmp[tmp$Edu_OR==1,],aes(x=BMI_OR,y=exp(beta),ymin=exp(low),ymax=exp(high),color=type))+geom_pointrange(position=position_dodge(0.5))+
  theme_minimal()+ylim(0,6)+geom_hline(yintercept=1, linetype="dashed", color="red")+
  geom_hline(yintercept=3.176839, linetype="solid", color="#F8766D")+geom_hline(yintercept=3.186489, linetype="solid", color="#00BA38")+
  geom_hline(yintercept=3.273303, linetype="solid", color="#619CFF")
dev.off()




pheno=fread("../T2D_phenotype.tsv")




qctool_v2.0-rc5 -g $temp1.bgen -incl-rsids $dose_snps -og $temp2.gen


res.test=matrix(nrow=3,ncol=7)

std.cov=paste0(paste0("pc",1:40,collapse="+"),"+g_bri+unrelated+batch+array")

k=1
phen$all=1
a="all"
i="t2d_or"
step1=summary(lm(as.formula(paste("bmi~bmi_all+sex+age+I(age^2)+",std.cov)),data=phen[which(phen[,..a]==1),]))
step2=try(summary(glm(as.formula(paste(i,"~bmi_all+age+sex+I(age^2)+",std.cov)),family=binomial(),data=phen[which(phen[,..a]==1),])))
beta.zy=step2$coefficients["bmi_all","Estimate"]
beta.zx=step1$coefficients["bmi_all","Estimate"]
se.zy=step2$coefficients["bmi_all","Std. Error"]
se.zx=step1$coefficients["bmi_all","Std. Error"]
cov.zx.zy=0
beta=beta.zy/beta.zx
se=((se.zy^2) / beta.zx^2) + (((beta.zy^2)/(beta.zx^4))*(se.zx^2)) - (2*(beta.zy/(beta.zx^3)) *cov.zx.zy)
se=sqrt(se)
c.c=table(model.frame(as.formula(paste(i,"~bmi_all+age+sex+I(age^2)+",std.cov)),data=phen[which(phen[,..a]==1),])[,c(i,"sex")])

res.test[1,]=c(i,beta,se,c.c[,1],c.c[,2])




phen$all=1
a="all"
i="t2d_or_f"

tmp.d=phen[which(phen[,..a]==1 & phen$sex==0),]
step1=summary(lm(as.formula(paste("bmi~bmi_female+age+I(age^2)+",std.cov)),data=tmp.d))
step2=try(summary(glm(as.formula(paste(i,"~bmi_female+age+I(age^2)+",std.cov)),family=binomial(),data=tmp.d)))
beta.zy=step2$coefficients["bmi_female","Estimate"]
beta.zx=step1$coefficients["bmi_female","Estimate"]
se.zy=step2$coefficients["bmi_female","Std. Error"]
se.zx=step1$coefficients["bmi_female","Std. Error"]
cov.zx.zy=0
beta=beta.zy/beta.zx
se=((se.zy^2) / beta.zx^2) + (((beta.zy^2)/(beta.zx^4))*(se.zx^2)) - (2*(beta.zy/(beta.zx^3)) *cov.zx.zy)
se=sqrt(se)
c.c=table(model.frame(as.formula(paste(i,"~bmi_female+age+I(age^2)+",std.cov)),data=tmp.d)[,i])

res.test[2,]=c(i,beta,se,c.c,NA,NA)



k=1
phen$all=1
a="all"
i="t2d_or_m"

tmp.d=phen[which(phen[,..a]==1 & phen$sex==1),]
step1=summary(lm(as.formula(paste("bmi~bmi_male+age+I(age^2)+",std.cov)),data=tmp.d))
step2=try(summary(glm(as.formula(paste(i,"~bmi_male+age+I(age^2)+",std.cov)),family=binomial(),data=tmp.d)))
beta.zy=step2$coefficients["bmi_male","Estimate"]
beta.zx=step1$coefficients["bmi_male","Estimate"]
se.zy=step2$coefficients["bmi_male","Std. Error"]
se.zx=step1$coefficients["bmi_male","Std. Error"]
cov.zx.zy=0
beta=beta.zy/beta.zx
se=((se.zy^2) / beta.zx^2) + (((beta.zy^2)/(beta.zx^4))*(se.zx^2)) - (2*(beta.zy/(beta.zx^3)) *cov.zx.zy)
se=sqrt(se)
c.c=table(model.frame(as.formula(paste(i,"~bmi_male+age+I(age^2)+",std.cov)),data=tmp.d)[,i])

res.test[3,]=c(i,beta,se,NA,NA,c.c)






