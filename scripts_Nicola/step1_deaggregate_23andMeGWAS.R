library(data.table)
## Load GWAS all and GWAS <30
all=fread("data/23andMe_SexGWAS_AllSamples.assoc")
young=fread("data/23andMe_SexGWAS_Age30.assoc")
all$Z=all$Effect/all$se
young$Z=young$Effect/young$se
## Allign 2 GWAS
all=all[all$SNP%in%young$SNP,]
young=young[young$SNP%in%all$SNP,]

No=all$N-young$N # Numeber of older people
wo=sqrt(No)
wy=sqrt(young$N)
Zo=((all$Z*sqrt(wo^2+wy^2))-(young$Z*wy))/sqrt(No)







old=all
old$Z=Zo
old$P=pchisq(old$Z^2,df=1,lower=F)
old$N=No
#cor(old$Z,all$Z)
#cor(old$Z,young$Z)

-14.51245*(wo+wy) -6.528811*wy[1]


(-0.03555801     -   -0.01153482)/(1/sqrt(No))

-0.04709283*sqrt(No)

## prepare fo LD score regression
snp.list=fread("../../app_full_stuff/ldsc/eur_w_ld_chr/w_hm3.snplist")
old.ldsc=old[old$SNP%in%snp.list$SNP,]
young.ldsc=young[young$SNP%in%snp.list$SNP,]

old.ldsc=old.ldsc[,.(SNP,A1,A2,Z,P,N)]
young.ldsc=young.ldsc[,.(SNP,A1,A2,Z,P,N)]
names(old.ldsc)=c("snpid","A1","A2","Zscore","P-value","N")
names(young.ldsc)=c("snpid","A1","A2","Zscore","P-value","N")

write.table(young.ldsc,file="young_vs_old/Twenty_three_young.tsv",row.names=F,quote=F,sep="\t")
write.table(old.ldsc,file="young_vs_old/Twenty_three_old.tsv",row.names=F,quote=F,sep="\t")


library(TwoSampleMR)
token="../prj_070_CVD/.httr-oauth"
source("../prj_063_Taste_genes/fix_model/scripts/MR_functions_Metabolites_MR.R")
ao=available_outcomes(token)

ids=c(1001,835)
exposures=extract_instruments(ids,access_token = token)
outcome.y=format_data(young,type = "outcome",snps = unique(exposures$SNP),snp_col ="SNP",beta_col = "Effect",se_col = "se"
                      ,effect_allele_col = "A1",other_allele_col = "A2",eaf_col = "effect_freq",pval_col = "P",samplesize_col = "N" )
outcome.o=format_data(all,type = "outcome",snps = unique(exposures$SNP),snp_col ="SNP",beta_col = "Effect",se_col = "se"
                      ,effect_allele_col = "A1",other_allele_col = "A2",eaf_col = "effect_freq",pval_col = "P",samplesize_col = "N" )

harm.bmi.y=harmonise_data(exposures[exposures$exposure=="Body mass index || id:835",],outcome.y)
harm.bmi.o=harmonise_data(exposures[exposures$exposure=="Body mass index || id:835",],outcome.o)
harm.edu.y=harmonise_data(exposures[exposures$exposure=="Years of schooling || id:1001",],outcome.y)
harm.edu.o=harmonise_data(exposures[exposures$exposure=="Years of schooling || id:1001",],outcome.o)

res.bmi.y=mr.run(harm.bmi.y)
res.bmi.o=mr.run(harm.bmi.o)
res.edu.y=mr.run(harm.edu.y)
res.edu.o=mr.run(harm.edu.o)
RadialMR::plot_radial(mr(harm.edu.o,method_list = "ivw"))


diff_Z=(all$Effect-young$Effect)/sqrt(all$se^2+young$se^2)
all$diff_Z=diff_Z
all$p_diff=pchisq(all$diff_Z,df=1,lower=F)

poss.diff.y=young[which(  all$P<5e-8),]
poss.diff.o=all[which(  all$P<5e-8),]
to.clump=data.frame(SNP=poss.diff.o$SNP, chr_name=poss.diff.o$chrom,chrom_start=poss.diff.o$position_b37, pval.exposure=poss.diff.o$P)
pruned=clump_data(dat = to.clump)

poss.diff.y=poss.diff.y[poss.diff.y$SNP%in%pruned$SNP]
poss.diff.o=poss.diff.o[poss.diff.o$SNP%in%pruned$SNP]
res <- phenoscanner(snpquery=different$SNP[abs(different$Effect.x/different$Effect.y)>1])
plot()

snpid	A1	A2	Zscore	N	P-value


head(all)
head(young)

wa=1/sqrt(N)

Zscore total= Z(<30) * w(<30) + Z(>30)*w(>30) /(w(<30)+w(>30))





res[,test:=NA]
res[,.(qual)]=gsub("/","_",res[,.(qual)])

res$test[res$qual=='College or University degree']=20
res$test[res$qual=='A levels/AS levels or equivalent']=15

recoded= recode(res$qual,'College or University degree'=20 ,'A levels/AS levels or equivalent'=15 , 'O levels/GCSEs or equivalent'=13 , 'CSEs or equivalent'=12 ,
                                 'NVQ or HND or HNC or equivalent'=19 , 'Other professional qualifications eg: nursing, teaching'=17 , 'None of the above'=6)




res[,test:=recoded]

