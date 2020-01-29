source("MR_functions_Metabolites_MR.R")
library(data.table)
library(TwoSampleMR)

### Create instrument object from Lindgren paper
all.instrument=fread("bmi_all_intrument.txt")
all.instr=format_data(dat=all.instrument,type ="exposure",snp_col = "SNP", beta_col = "beta.combined", 
            se_col = "se.combined", eaf_col = "frqA1.combined", effect_allele_col = "EA", 
            other_allele_col = "NEA", pval_col = "pval.combined", 
            samplesize_col = "nmeta.combined",
            info_col = "info.combined", chr_col = "Chr", pos_col = "Pos")

all.instr$exposure="BMI_both_sex"
all.instr$SNP=apply(t(all.instr$SNP),2,function(x)unlist(strsplit(x,split="_"))[1])



men.instrument=fread("bmi_men_intrument.txt")
men.instr=format_data(dat=men.instrument,type ="exposure",snp_col = "SNP", beta_col = "beta.males", 
                      se_col = "se.males", eaf_col = "frqA1.males", effect_allele_col = "EA", 
                      other_allele_col = "NEA", pval_col = "pval.males", 
                      samplesize_col = "nmeta.males",
                      info_col = "info.males", chr_col = "Chr", pos_col = "Pos")

men.instr$exposure="BMI_male"
men.instr$SNP=apply(t(men.instr$SNP),2,function(x)unlist(strsplit(x,split="_"))[1])


women.instrument=fread("bmi_women_intrument.txt")
women.instr=format_data(dat=women.instrument,type ="exposure",snp_col = "SNP", beta_col = "beta.females", 
                      se_col = "se.females", eaf_col = "frqA1.females", effect_allele_col = "EA", 
                      other_allele_col = "NEA", pval_col = "pval.females", 
                      samplesize_col = "nmeta.females",
                      info_col = "info.females", chr_col = "Chr", pos_col = "Pos")

women.instr$exposure="BMI_female"
women.instr$SNP=apply(t(women.instr$SNP),2,function(x)unlist(strsplit(x,split="_"))[1])

phenos=fread("../GWAS/T2D_sim_BMI/phenotype_simul.txt")


### Create outcome object from UKB
## T2D no sex
t2d.all=fread("../GWAS/T2D_sim_BMI/p05_comb_chr/t2d_or.tsv.gz")
prev=mean(phenos$fc1_poss_t2dm,na.rm=T)
t2d.all.data=format_data(snps =unique(c(women.instr$SNP,men.instr$SNP,all.instr$SNP)) ,dat=t2d.all,type ="outcome",snp_col = "rsid", beta_col = "beta1", 
                         se_col = "se", eaf_col = "freq1", effect_allele_col = "a1", 
                         other_allele_col = "a0", pval_col = "p", 
                         samplesize_col = "n",
                         info_col = "info_pop", chr_col = "chr", pos_col = "pos")
t2d.all.data$outcome="T2D_OR"
t2d.all.data$beta.outcome=t2d.all.data$beta.outcome/(prev*(1-prev))
t2d.all.data$se.outcome=t2d.all.data$se.outcome/(prev*(1-prev))

## T2D  sex

t2d.all=fread("../GWAS/T2D_sim_BMI/p05_comb_chr/t2d_sex.tsv.gz")
prev=mean(phenos$fc1_poss_t2dm,na.rm=T)

t2d.all.sex.data=format_data(snps =unique(c(women.instr$SNP,men.instr$SNP,all.instr$SNP)) ,dat=t2d.all,type ="outcome",snp_col = "rsid", beta_col = "beta1", 
                         se_col = "se", eaf_col = "freq1", effect_allele_col = "a1", 
                         other_allele_col = "a0", pval_col = "p", 
                         samplesize_col = "n",
                         info_col = "info_pop", chr_col = "chr", pos_col = "pos")
t2d.all.sex.data$outcome="T2D_OR_sex"
t2d.all.sex.data$beta.outcome=t2d.all.sex.data$beta.outcome/(prev*(1-prev))
t2d.all.sex.data$se.outcome=t2d.all.sex.data$se.outcome/(prev*(1-prev))



## T2D original Male
t2d.all=fread("../GWAS/T2D_sim_BMI/p05_comb_chr/t2d_m.tsv.gz")

t2d.male.data=format_data(snps =unique(c(women.instr$SNP,men.instr$SNP,all.instr$SNP)) ,dat=t2d.all,type ="outcome",snp_col = "rsid", beta_col = "beta1", 
                             se_col = "se", eaf_col = "freq1", effect_allele_col = "a1", 
                             other_allele_col = "a0", pval_col = "p", 
                             samplesize_col = "n",
                             info_col = "info_pop", chr_col = "chr", pos_col = "pos")
t2d.male.data$outcome="T2D_OR_male"

prev=mean(phenos$t2d.m,na.rm=T)
t2d.male.data$beta.outcome=t2d.male.data$beta.outcome/(prev*(1-prev))
t2d.male.data$se.outcome=t2d.male.data$se.outcome/(prev*(1-prev))

## T2D original Female
t2d.all=fread("../GWAS/T2D_sim_BMI/p05_comb_chr/t2d_f.tsv.gz")

t2d.female.data=format_data(snps =unique(c(women.instr$SNP,men.instr$SNP,all.instr$SNP)) ,dat=t2d.all,type ="outcome",snp_col = "rsid", beta_col = "beta1", 
                          se_col = "se", eaf_col = "freq1", effect_allele_col = "a1", 
                          other_allele_col = "a0", pval_col = "p", 
                          samplesize_col = "n",
                          info_col = "info_pop", chr_col = "chr", pos_col = "pos")
t2d.female.data$outcome="T2D_OR_female"

prev=mean(phenos$t2d.f,na.rm=T)
t2d.female.data$beta.outcome=t2d.female.data$beta.outcome/(prev*(1-prev))
t2d.female.data$se.outcome=t2d.female.data$se.outcome/(prev*(1-prev))


### T2D bias no sex
t2d.all=fread("../GWAS/T2D_sim_BMI/p05_comb_chr/t2d_new.tsv.gz")

t2d.bias.data=format_data(snps =unique(c(women.instr$SNP,men.instr$SNP,all.instr$SNP)) ,dat=t2d.all,type ="outcome",snp_col = "rsid", beta_col = "beta1", 
                            se_col = "se", eaf_col = "freq1", effect_allele_col = "a1", 
                            other_allele_col = "a0", pval_col = "p", 
                            samplesize_col = "n",
                            info_col = "info_pop", chr_col = "chr", pos_col = "pos")
t2d.bias.data$outcome="T2D_bias"

prev=mean(phenos$t2d.new,na.rm=T)
t2d.bias.data$beta.outcome=t2d.bias.data$beta.outcome/(prev*(1-prev))
t2d.bias.data$se.outcome=t2d.bias.data$se.outcome/(prev*(1-prev))



### T2D bias sex
t2d.all=fread("../GWAS/T2D_sim_BMI/p05_comb_chr/t2d_new_sex.tsv.gz")

t2d.bias.sex.data=format_data(snps =unique(c(women.instr$SNP,men.instr$SNP,all.instr$SNP)) ,dat=t2d.all,type ="outcome",snp_col = "rsid", beta_col = "beta1", 
                          se_col = "se", eaf_col = "freq1", effect_allele_col = "a1", 
                          other_allele_col = "a0", pval_col = "p", 
                          samplesize_col = "n",
                          info_col = "info_pop", chr_col = "chr", pos_col = "pos")
t2d.bias.sex.data$outcome="T2D_bias_sex"
prev=mean(phenos$t2d.new,na.rm=T)
t2d.bias.sex.data$beta.outcome=t2d.bias.sex.data$beta.outcome/(prev*(1-prev))
t2d.bias.sex.data$se.outcome=t2d.bias.sex.data$se.outcome/(prev*(1-prev))


### T2D bias Male
t2d.all=fread("../GWAS/T2D_sim_BMI/p05_comb_chr/t2d_m_new.tsv.gz")

t2d.bias.male.data=format_data(snps =unique(c(women.instr$SNP,men.instr$SNP,all.instr$SNP)) ,dat=t2d.all,type ="outcome",snp_col = "rsid", beta_col = "beta1", 
                              se_col = "se", eaf_col = "freq1", effect_allele_col = "a1", 
                              other_allele_col = "a0", pval_col = "p", 
                              samplesize_col = "n",
                              info_col = "info_pop", chr_col = "chr", pos_col = "pos")
t2d.bias.male.data$outcome="T2D_bias_Male"

prev=mean(phenos$t2d.new.m,na.rm=T)
t2d.bias.male.data$beta.outcome=t2d.bias.male.data$beta.outcome/(prev*(1-prev))
t2d.bias.male.data$se.outcome=t2d.bias.male.data$se.outcome/(prev*(1-prev))



### T2D bias Male
t2d.all=fread("../GWAS/T2D_sim_BMI/p05_comb_chr/t2d_f_new.tsv.gz")

t2d.bias.female.data=format_data(snps =unique(c(women.instr$SNP,men.instr$SNP,all.instr$SNP)) ,dat=t2d.all,type ="outcome",snp_col = "rsid", beta_col = "beta1", 
                               se_col = "se", eaf_col = "freq1", effect_allele_col = "a1", 
                               other_allele_col = "a0", pval_col = "p", 
                               samplesize_col = "n",
                               info_col = "info_pop", chr_col = "chr", pos_col = "pos")
t2d.bias.female.data$outcome="T2D_bias_Female"
prev=mean(phenos$t2d.new.f,na.rm=T)
t2d.bias.female.data$beta.outcome=t2d.bias.female.data$beta.outcome/(prev*(1-prev))
t2d.bias.female.data$se.outcome=t2d.bias.female.data$se.outcome/(prev*(1-prev))



### BMI -> T2D origin
res.tot=c()
harm=harmonise_data(all.instr,t2d.all.data )
res=mr.run(harm)
res.tot=rbind(res.tot,res)

### BMI -> T2D origin sex
harm=harmonise_data(all.instr,t2d.all.sex.data )
res=mr.run(harm)
res.tot=rbind(res.tot,res)

### BMI -> T2D bias

harm=harmonise_data(all.instr,t2d.bias.data )
res=mr.run(harm)
res.tot=rbind(res.tot,res)

t2d.bias.data

### BMI -> T2D bias sex

harm=harmonise_data(all.instr,t2d.bias.sex.data )
res=mr.run(harm)
res.tot=rbind(res.tot,res)

####Sex separated analyses: Male
# Original
harm=harmonise_data(men.instr,t2d.male.data )
res=mr.run(harm)
res.tot=rbind(res.tot,res)

# Bias
harm=harmonise_data(men.instr,t2d.bias.male.data )
res=mr.run(harm)
res.tot=rbind(res.tot,res)

# Original
harm=harmonise_data(women.instr,t2d.female.data )
res=mr.run(harm)
res.tot=rbind(res.tot,res)

# Bias
harm=harmonise_data(women.instr,t2d.bias.female.data )
res=mr.run(harm)
res.tot=rbind(res.tot,res)

res.tot$cat=c("All","All","All","All","Male","Male","Female","Female")


library(ggplot2)


ext.t2d=extract_outcome_data(snps=unique(c(women.instr$SNP,men.instr$SNP,all.instr$SNP)),outcomes =26,access_token = ".httr-oauth")
harm=harmonise_data(all.instr,ext.t2d)
ext.res=mr.run(harm)
ext.res$cat="External"
res.tot=rbind(res.tot,ext.res)

res.tot$low=res.tot$b+(qnorm(0.025)*res.tot$se)
res.tot$high=res.tot$b+(qnorm(0.025,lower=F)*res.tot$se)

library(ggplot2)
pdf("MR_results_BMI_bias.pdf",height=7)
ggplot(res.tot,aes(y=exp(b),ymin=exp(low),ymax=exp(high),x=outcome,fill=cat))+geom_pointrange(pch=21,size=1)+
  coord_flip()+theme_minimal()+ylim(0,4)+geom_hline(yintercept=1)+geom_text(aes(label = round(exp(b),digits=2),vjust=-1))
dev.off()
















t2d.bias.sex.data


