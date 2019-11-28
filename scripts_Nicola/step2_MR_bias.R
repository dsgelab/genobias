library(data.table)
height.sex=fread("simul/gwas/height_or_sex.tsv.gz",data.table=F)
height8.dummy.sex=fread("simul/gwas/height_or8_dummy.tsv.gz",data.table=F)
height8.dummy.plain=fread("simul/gwas/height_or8_plain.tsv.gz",data.table=F)
mcv.dummy8.sex=fread("simul/gwas/mcv_or8_dummy.tsv.gz",data.table=F)
mcv.dummy8.plain=fread("simul/gwas/mcv_or8_plain.tsv.gz",data.table=F)
height4.dummy.sex=fread("simul/gwas/height_dummy.tsv.gz",data.table=F)
height4.dummy.plain=fread("simul/gwas/height_plain.tsv.gz",data.table=F)
mcv.dummy4.sex=fread("simul/gwas/mcv_dummy.tsv.gz",data.table=F)
mcv.dummy4.plain=fread("simul/gwas/mcv_plain.tsv.gz",data.table=F)
mcv.sex=fread("simul/gwas/mcv_or_sex.tsv.gz",data.table=F)

  
pval.table=data.frame(rsid=height.sex$rsid,p.height.or=height.sex$p,p.mcv.or=mcv.sex$p,
                      height8.dummy.sex=height8.dummy.sex$p,height8.dummy.plain=height8.dummy.plain$p,
                      mcv.dummy8.sex=mcv.dummy8.sex$p,mcv.dummy8.plain=mcv.dummy8.plain$p,
                      height4.dummy.sex=height4.dummy.sex$p,height4.dummy.plain=height4.dummy.plain$p,
                      mcv.dummy4.sex=mcv.dummy4.sex$p,mcv.dummy4.plain=mcv.dummy4.plain$p,stringsAsFactors=F)



idx=which(pval.table$p.height.or>0.01& pval.table$height4.dummy.sex<5e-5 & pval.table$mcv.dummy4.plain<5e-8)

instrument.create=function(data,phen.name="test",p=5e-8){
  require(TwoSampleMR)
  data=data[which(data$p<p),]
  data$phen=phen.name
  data.instr=format_data(data,type = "exposure",phenotype_col = "phen"
                         ,snp_col = "rsid",beta_col = "beta1",se_col = "se"
                         ,eaf_col = "freq1",effect_allele_col = "a1",other_allele_col = "a0",pval_col = "p",samplesize_col = "n"
                         ,chr_col = "chr",pos_col = "pos")
  data.instr=clump_data(data.instr)
  
}



### MCV 4 cause sex
## Creat MCV 4 instrument
library(TwoSampleMR)
height.or=instrument.create(height.sex,phen.name="height.or",p=1e-15)
mcv.or=instrument.create(mcv.sex,phen.name="mcv.or")
mcv.4.dsex.instr=instrument.create(mcv.dummy4.sex,phen.name="mcv4.dsex")    ##mcv or 4 ~ dsex
mcv.8.dsex.instr=instrument.create(mcv.dummy8.sex,phen.name="mcv8.dsex")    ##mcv or 8 ~ dsex
height.8.dsex.instr=instrument.create(height8.dummy.sex,phen.name="height8.dsex") ##height or 4 ~ dsex
height.4.dsex.instr=instrument.create(height4.dummy.sex,phen.name="height4.dsex") ##height or 8 ~ dsex
mcv.4.instr=instrument.create(mcv.dummy4.plain,phen.name="mcv4.plain")  ##mcv or 4 
mcv.8.instr=instrument.create(mcv.dummy8.plain,phen.name="mcv8.plain")  ##mcv or 8 
height.8.instr=instrument.create(height8.dummy.plain,phen.name="height8.plain")  ##height or 4 
height.4.instr=instrument.create(height4.dummy.plain,phen.name="height4.plain")  ##height or 8 ~ dsex
ext.instr=extract_instruments(outcomes = c(89,1250),access_token = "../../../f001_ciara_supervisormeeting/MR/Two_way_MR/Part_B/Z_scores_based/CAD/mrbase.oauth",p1 = 6e-8) ##Wood et al height




overall.snp.list=unique(c(mcv.4.dsex.instr$SNP,mcv.8.dsex.instr$SNP,height.8.dsex.instr$SNP,height.4.dsex.instr$SNP,
                          mcv.4.instr$SNP,mcv.8.instr$SNP,height.8.instr$SNP,height.4.instr$SNP,ext.instr$SNP))


#### Outcomes
height.or.outcome=format_data(height.sex,type = "outcome",phenotype_col = "phen"
                              ,snp_col = "rsid",beta_col = "beta1",se_col = "se"
                              ,eaf_col = "freq1",effect_allele_col = "a1",other_allele_col = "a0",pval_col = "p",samplesize_col = "n"
                              ,chr_col = "chr",pos_col = "pos",snps =overall.snp.list )
height.or.outcome$outcome="height.or"

mcv.or.outcome=format_data(mcv.sex,type = "outcome",phenotype_col = "phen"
                              ,snp_col = "rsid",beta_col = "beta1",se_col = "se"
                              ,eaf_col = "freq1",effect_allele_col = "a1",other_allele_col = "a0",pval_col = "p",samplesize_col = "n"
                              ,chr_col = "chr",pos_col = "pos",snps =overall.snp.list )
mcv.or.outcome$outcome="mcv.or"


mcv.4.dsex.outcome=format_data(mcv.dummy4.sex,type = "outcome"
                              ,snp_col = "rsid",beta_col = "beta1",se_col = "se"
                              ,eaf_col = "freq1",effect_allele_col = "a1",other_allele_col = "a0",pval_col = "p",samplesize_col = "n"
                              ,chr_col = "chr",pos_col = "pos",snps =overall.snp.list )
mcv.4.dsex.outcome$outcome="mcv4.dsex"

mcv.8.dsex.outcome=format_data(mcv.dummy8.sex,type = "outcome"
                          ,snp_col = "rsid",beta_col = "beta1",se_col = "se"
                          ,eaf_col = "freq1",effect_allele_col = "a1",other_allele_col = "a0",pval_col = "p",samplesize_col = "n"
                          ,chr_col = "chr",pos_col = "pos",snps =overall.snp.list )
mcv.8.dsex.outcome$outcome="mcv8.dsex"

height.4.dsex.outcome=format_data(height4.dummy.sex,type = "outcome"
                          ,snp_col = "rsid",beta_col = "beta1",se_col = "se"
                          ,eaf_col = "freq1",effect_allele_col = "a1",other_allele_col = "a0",pval_col = "p",samplesize_col = "n"
                          ,chr_col = "chr",pos_col = "pos",snps =overall.snp.list )
height.4.dsex.outcome$outcome="height4.dsex"

height.8.dsex.outcome=format_data(height8.dummy.sex,type = "outcome"
                             ,snp_col = "rsid",beta_col = "beta1",se_col = "se"
                             ,eaf_col = "freq1",effect_allele_col = "a1",other_allele_col = "a0",pval_col = "p",samplesize_col = "n"
                             ,chr_col = "chr",pos_col = "pos",snps =overall.snp.list )
height.8.dsex.outcome$outcome="height8.dsex"

mcv.4.plain.outcome=format_data(mcv.dummy4.plain,type = "outcome"
                          ,snp_col = "rsid",beta_col = "beta1",se_col = "se"
                          ,eaf_col = "freq1",effect_allele_col = "a1",other_allele_col = "a0",pval_col = "p",samplesize_col = "n"
                          ,chr_col = "chr",pos_col = "pos",snps =overall.snp.list )
mcv.4.plain.outcome$outcome="mcv4.plain"

mcv.8.plain.outcome=format_data(mcv.dummy8.plain,type = "outcome"
                          ,snp_col = "rsid",beta_col = "beta1",se_col = "se"
                          ,eaf_col = "freq1",effect_allele_col = "a1",other_allele_col = "a0",pval_col = "p",samplesize_col = "n"
                          ,chr_col = "chr",pos_col = "pos",snps =overall.snp.list )
mcv.8.plain.outcome$outcome="mcv8.plain"

height.4.plain.outcome=format_data(height4.dummy.plain,type = "outcome"
                             ,snp_col = "rsid",beta_col = "beta1",se_col = "se"
                             ,eaf_col = "freq1",effect_allele_col = "a1",other_allele_col = "a0",pval_col = "p",samplesize_col = "n"
                             ,chr_col = "chr",pos_col = "pos",snps =overall.snp.list )
height.4.plain.outcome$outcome="height4.plain"

height.8.plain.outcome=format_data(height8.dummy.plain,type = "outcome"
                             ,snp_col = "rsid",beta_col = "beta1",se_col = "se"
                             ,eaf_col = "freq1",effect_allele_col = "a1",other_allele_col = "a0",pval_col = "p",samplesize_col = "n"
                             ,chr_col = "chr",pos_col = "pos",snps =overall.snp.list )
height.8.plain.outcome$outcome="height8.plain"

### Extract external gwas
#ao=available_outcomes(access_token ="../../../f001_ciara_supervisormeeting/MR/Two_way_MR/Part_B/Z_scores_based/CAD/mrbase.oauth")
#ao[grep("height",tolower(ao$trait)),]

#id 89 height
#id 1250 MCV

ext.outcome=extract_outcome_data(snps=overall.snp.list,outcomes=c(89,1250),access_token = "../../../f001_ciara_supervisormeeting/MR/Two_way_MR/Part_B/Z_scores_based/CAD/mrbase.oauth")

### Baseline MR

harm=harmonise_data(ext.instr,ext.outcome)
harm=harm[-which(harm$exposure==harm$outcome),]
st.filtered=steiger_filtering(harm)
library(mr.radial)
baseline.mr=mr(st.filtered)

### merge outcomes
colo=intersect(names(mcv.4.dsex.outcome),names(ext.outcome))
overall.outcome=rbind(height.or.outcome[,colo]
                      ,mcv.or.outcome[,colo]
                      ,mcv.4.dsex.outcome[,colo]
                      ,mcv.8.dsex.outcome[,colo]
                      ,height.4.dsex.outcome[,colo]
                      ,height.8.dsex.outcome[,colo]
                      ,mcv.4.plain.outcome[,colo]
                      ,mcv.8.plain.outcome[,colo]
                      ,height.4.plain.outcome[,colo]
                      ,height.8.plain.outcome[,colo]
                      ,ext.outcome[,colo])
overall.outcome$units.outcome="SD"

colo=intersect(names(mcv.4.dsex.instr),names(ext.instr))
overall.instrument=rbind(height.or[,colo]
                         ,mcv.or[,colo]
                         ,mcv.4.dsex.instr[,colo]
                         ,mcv.8.dsex.instr[,colo]
                         ,height.8.dsex.instr[,colo]
                         ,height.4.dsex.instr[,colo]
                         ,mcv.4.instr[,colo]
                         ,mcv.8.instr[,colo]
                         ,height.8.instr[,colo]
                         ,height.4.instr[,colo]
                         ,ext.instr[,colo])
overall.instrument$units.exposure="SD"



harm.overall=harmonise_data(overall.instrument,overall.outcome)
st.filtered.ov=steiger_filtering(harm.overall)
overall.mr=mr(st.filtered.ov)



#### RG edu codes
load("simul/MR/overall_MR.Rdata")
mr.mcv=overall.mr[intersect(which(overall.mr$method=="Inverse variance weighted" & 
                   overall.mr$exposure=="Mean corpuscular volume || id:1250"),
                   grep("height|Height",overall.mr$outcome)),]

mr.mcv4.plain=overall.mr[intersect(which(overall.mr$method=="Inverse variance weighted" & 
                                    overall.mr$exposure=="mcv4.plain"),
                            grep("height|Height",overall.mr$outcome)),]

mr.mcv8.plain=overall.mr[intersect(which(overall.mr$method=="Inverse variance weighted" & 
                                           overall.mr$exposure=="mcv8.plain"),
                                   grep("height|Height",overall.mr$outcome)),]

mr.mcv4.dsex=overall.mr[intersect(which(overall.mr$method=="Inverse variance weighted" & 
                                           overall.mr$exposure=="mcv4.dsex"),
                                   grep("height|Height",overall.mr$outcome)),]


mr.mcv8.dsex=overall.mr[intersect(which(overall.mr$method=="Inverse variance weighted" & 
                                          overall.mr$exposure=="mcv8.dsex"),
                                  grep("height|Height",overall.mr$outcome)),]
mr.mcv.all=rbind(mr.mcv,mr.mcv4.plain,mr.mcv4.dsex,mr.mcv8.plain,mr.mcv8.dsex)

mr.mcv.all$hig=mr.mcv.all$b+mr.mcv.all$se*qnorm(1-0.025)
mr.mcv.all$low=mr.mcv.all$b+mr.mcv.all$se*qnorm(0.025)
mr.mcv.all$exposure=recode(mr.mcv.all$exposure,"Mean corpuscular volume || id:1250"="MCV Astle",
                           mcv4.plain="MCV OR 4",mcv4.dsex="MCV OR 4~sex",mcv8.plain="MCV OR 8",mcv8.dsex="MCV OR 8 ~ sex")

mr.mcv.all$exposure=factor(mr.mcv.all$exposure,levels=rev(levels(mr.mcv.all$exposure)))
mr.mcv.all$outcome=recode(mr.mcv.all$outcome,"Height || id:89"="Height Woode",height.or="Height Ukb",height4.plain="Height OR 4",
                              height4.dsex="Height OR 4 ~ sex",height8.plain="Height OR 8",height8.dsex="Height OR 8 ~  sex")
mr.mcv.all$outcome=factor(mr.mcv.all$outcome,levels=rev(c("Height Woode","Height Ukb","Height OR 4",
                             "Height OR 4 ~ sex","Height OR 8","Height OR 8 ~  sex")))


library(ggplot2)
pdf("simul/MR/MCV_exposure_forestplot.pdf")
ggplot(data=mr.mcv.all,aes(x=outcome,y=b,ymin=low,ymax=hig,color=exposure))+geom_pointrange(position = position_dodge(width = 0.50))+
  coord_flip()+ theme_minimal()+ggtitle("MCV exposure")+geom_hline(yintercept = 0,color="darkgrey") 
dev.off()





##height

mr.height=overall.mr[intersect(which(overall.mr$method=="Inverse variance weighted" & 
                                       overall.mr$exposure=="Height || id:89"),
                               grep("mcv|corpuscular",overall.mr$outcome)),]

mr.height$hig=mr.height$b+mr.height$se*qnorm(1-0.025)
mr.height$low=mr.height$b+mr.height$se*qnorm(0.025)
library(ggplot2)
mr.height$outcome=factor(mr.height$outcome,levels=rev(c("Mean corpuscular volume || id:1250","mcv.or","mcv4.plain","mcv4.dsex",
                                                        "mcv8.plain","mcv8.dsex")))



mr.height4.plain=overall.mr[intersect(which(overall.mr$method=="Inverse variance weighted" & 
                                       overall.mr$exposure=="height4.plain"),
                               grep("mcv|corpuscular",overall.mr$outcome)),]

mr.height4.plain$hig=mr.height4.plain$b+mr.height4.plain$se*qnorm(1-0.025)
mr.height4.plain$low=mr.height4.plain$b+mr.height4.plain$se*qnorm(0.025)
library(ggplot2)
mr.height4.plain$outcome=factor(mr.height4.plain$outcome,levels=rev(c("Mean corpuscular volume || id:1250","mcv.or","mcv4.plain","mcv4.dsex",
                                                        "mcv8.plain","mcv8.dsex")))

mr.height4.dsex=overall.mr[intersect(which(overall.mr$method=="Inverse variance weighted" & 
                                              overall.mr$exposure=="height4.dsex"),
                                      grep("mcv|corpuscular",overall.mr$outcome)),]
mr.height4.dsex$hig=mr.height4.dsex$b+mr.height4.dsex$se*qnorm(1-0.025)
mr.height4.dsex$low=mr.height4.dsex$b+mr.height4.dsex$se*qnorm(0.025)
library(ggplot2)
mr.height4.dsex$outcome=factor(mr.height4.dsex$outcome,levels=rev(c("Mean corpuscular volume || id:1250","mcv.or","mcv4.plain","mcv4.dsex",
                                                                      "mcv8.plain","mcv8.dsex")))

ggplot(data=mr.height4.dsex,aes(x=outcome,y=b,ymin=low,ymax=hig))+geom_pointrange()+
  coord_flip()+ theme_minimal()+geom_hline(yintercept = 0,color="darkgrey") 





mr.height8.plain=overall.mr[intersect(which(overall.mr$method=="Inverse variance weighted" & 
                                              overall.mr$exposure=="height8.plain"),
                                      grep("mcv|corpuscular",overall.mr$outcome)),]

mr.height8.plain$hig=mr.height8.plain$b+mr.height8.plain$se*qnorm(1-0.025)
mr.height8.plain$low=mr.height8.plain$b+mr.height8.plain$se*qnorm(0.025)
library(ggplot2)
mr.height8.plain$outcome=factor(mr.height8.plain$outcome,levels=rev(c("Mean corpuscular volume || id:1250","mcv.or","mcv4.plain","mcv4.dsex",
                                                                    "mcv8.plain","mcv8.dsex")))

ggplot(data=mr.height8.plain,aes(x=outcome,y=b,ymin=low,ymax=hig))+geom_pointrange()+
  coord_flip()+ theme_minimal()+ggtitle("Height OR8 ")




mr.height8.dsex=overall.mr[intersect(which(overall.mr$method=="Inverse variance weighted" & 
                                             overall.mr$exposure=="height8.dsex"),
                                     grep("mcv|corpuscular",overall.mr$outcome)),]


mr.height8.dsex$hig=mr.height8.dsex$b+mr.height8.dsex$se*qnorm(1-0.025)
mr.height8.dsex$low=mr.height8.dsex$b+mr.height8.dsex$se*qnorm(0.025)
library(ggplot2)
mr.height8.dsex$outcome=factor(mr.height8.dsex$outcome,levels=rev(c("Mean corpuscular volume || id:1250","mcv.or","mcv4.plain","mcv4.dsex",
                                                                      "mcv8.plain","mcv8.dsex")))

ggplot(data=mr.height8.dsex,aes(x=outcome,y=b,ymin=low,ymax=hig))+geom_pointrange()+
  coord_flip()+ theme_minimal()+ggtitle("Height OR8 ~ dummy sex")



mr.height.ukb=rbind(mr.height,mr.height4.plain,mr.height4.dsex,mr.height8.plain,mr.height8.dsex)
mr.height.ukb$exposure=factor(mr.height.ukb$exposure,levels=c("height8.dsex","height8.plain","height4.dsex","height4.plain","Height || id:89"))
mr.height.ukb$exposure=recode(mr.height.ukb$exposure,"Height || id:89"="Height Woode",height4.plain="Height OR 4",
                              height4.dsex="Height OR 4 ~ sex",height8.plain="Height OR 8",height8.dsex="Height OR 8 ~  sex")
mr.height.ukb$outcome=as.character(mr.height.ukb$outcome)
mr.height.ukb$outcome[mr.height.ukb$outcome=="Mean corpuscular volume || id:1250"]="MCV"
mr.height.ukb$outcome=factor(mr.height.ukb$outcome,rev(c("MCV","mcv.or","mcv4.plain","mcv4.dsex",
      "mcv8.plain","mcv8.dsex")))
mr.height.ukb$outcome=recode(mr.height.ukb$outcome,MCV="MCV Astle 2016",
                             mcv.or="MCV Ukb",mcv4.plain="MCV OR 4",mcv4.dsex="MCV OR 4 ~  sex",mcv8.plain="MCV OR 8",mcv8.dsex="MCV OR 8 ~  sex")

pdf("simul/MR/Height_exposure_forestplot.pdf")

ggplot(data=mr.height.ukb,aes(x=outcome,y=b,ymin=low,ymax=hig,color=exposure))+geom_pointrange(position = position_dodge(width = 0.50))+
  coord_flip()+ theme_minimal()+ggtitle("Height exposure")+geom_hline(yintercept = 0,color="darkgrey") 
dev.off()

###### MR plots
library(TwoSampleMR)
load("simul/MR/harmonised_steig_data.Rdata")

p1=mr_scatter_plot(dat = st.filtered.ov[st.filtered.ov$exposure=="Mean corpuscular volume || id:1250" &st.filtered.ov$outcome=="height4.plain", ],mr_results=mr(st.filtered.ov[st.filtered.ov$exposure=="Mean corpuscular volume || id:1250" &st.filtered.ov$outcome=="height4.plain", ]))
p2=mr_scatter_plot(dat = st.filtered.ov[st.filtered.ov$exposure=="Mean corpuscular volume || id:1250" &st.filtered.ov$outcome=="height4.dsex", ],mr_results=mr(st.filtered.ov[st.filtered.ov$exposure=="Mean corpuscular volume || id:1250" &st.filtered.ov$outcome=="height4.dsex", ]))
p3=mr_scatter_plot(dat = st.filtered.ov[st.filtered.ov$exposure=="Mean corpuscular volume || id:1250" &st.filtered.ov$outcome=="height8.plain", ],mr_results=mr(st.filtered.ov[st.filtered.ov$exposure=="Mean corpuscular volume || id:1250" &st.filtered.ov$outcome=="height8.plain", ]))
p4=mr_scatter_plot(dat = st.filtered.ov[st.filtered.ov$exposure=="Mean corpuscular volume || id:1250" &st.filtered.ov$outcome=="height8.dsex", ],mr_results=mr(st.filtered.ov[st.filtered.ov$exposure=="Mean corpuscular volume || id:1250" &st.filtered.ov$outcome=="height8.dsex", ]))

library(ggpubr)
pdf("simul/MR/Scatter_MCV.pdf",height=14,width=14)
ggarrange(p1$`1250.Pznd10`+ylim(-0.055,0.07),p2$`1250.2eYU78`+ylim(-0.055,0.07),p3$`1250.cYnUyH`+ylim(-0.055,0.07),p4$`1250.lbuwme`+ylim(-0.055,0.07) )
dev.off()


p1=mr_scatter_plot(dat = st.filtered.ov[st.filtered.ov$exposure=="Height || id:89" &st.filtered.ov$outcome=="mcv4.plain", ],mr_results=mr(st.filtered.ov[st.filtered.ov$exposure=="Height || id:89" &st.filtered.ov$outcome=="mcv4.plain", ]))
p2=mr_scatter_plot(dat = st.filtered.ov[st.filtered.ov$exposure=="Height || id:89" &st.filtered.ov$outcome=="mcv4.dsex", ],mr_results=mr(st.filtered.ov[st.filtered.ov$exposure=="Height || id:89" &st.filtered.ov$outcome=="mcv4.dsex", ]))
p3=mr_scatter_plot(dat = st.filtered.ov[st.filtered.ov$exposure=="Height || id:89" &st.filtered.ov$outcome=="mcv8.plain", ],mr_results=mr(st.filtered.ov[st.filtered.ov$exposure=="Height || id:89" &st.filtered.ov$outcome=="mcv8.plain", ]))
p4=mr_scatter_plot(dat = st.filtered.ov[st.filtered.ov$exposure=="Height || id:89" &st.filtered.ov$outcome=="mcv8.dsex", ],mr_results=mr(st.filtered.ov[st.filtered.ov$exposure=="Height || id:89" &st.filtered.ov$outcome=="mcv8.dsex", ]))

pdf("simul/MR/Scatter_Height.pdf",height=14,width=14)
ggarrange(p1$`89.OGqI2u`+ylim(-0.04,0.05)+xlim(0,0.15),p2$`89.9ZYIcD`+ylim(-0.04,0.05)+xlim(0,0.15),p3$`89.WAKBZa`+ylim(-0.04,0.05)+xlim(0,0.15),p4$`89.AwDn54`+ylim(-0.04,0.05)+xlim(0,0.15) )
dev.off()




#### Edu rgs

file.list=system("ls simul/rgs/*",intern=T)
res.tot=c()
for(i in file.list){

  res=fread(i,data.table=F,skip=60,nrows = 1)
  res.tot=rbind(res.tot,res)

}

res.tot$p1=gsub("simul/files_ganna/","",gsub(".LDSC.sumstats.gz","",res.tot$p1))
res.tot$p2=recode(res.tot$p2,"simul/ldhub/edu_plain_mung.sumstats.gz"="Edu plain",
                  "simul/ldhub/edu_sex_mung.sumstats.gz"="Edu~sex")


plot(res.tot$rg[res.tot$p2=="Edu plain"],res.tot$rg[res.tot$p2=="Edu~sex"])


