library(bGWAS)
library(data.table)
source("../prj_063_Taste_genes/fix_model/scripts/Functions_MR_correction.R")

#R --vanilla --slave -q -f ../../scripts/step3_MR_correction.R --args 

tratto="Female"

trait.f=c("EDUyears_2016_sumstat.txt","All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq.gz",
          "tag.cpd.tbl.gz","tag.evrsmk.tbl.gz","tag.former.tbl.gz","bip1.scz1.ruderfer2014..ANY",
          "GPC_1.NEO_OPENNESS.full.txt","GPC_2.EXTRAVERSION.zip","GPC_2.NEUROTICISM.zip")

mystudies=select_priorGWASs(Z_matrices = "ZMatrices/",
                            include_files =trait.f)

r2.table=c()
all.mr.betas=c()
res.t=fread("data/23andMe_SexGWAS_AllSamples.assoc",data.table=F)
names(res.t)[c(1,3,2,5,6)]=c("rsid","a1","a0","beta1","se")
doppi=table(res.t$rsid)
doppi=names(doppi)[doppi>1]
res.t=res.t[which(!(res.t$rsid%in%doppi)),]

res.bg2=  bGWAS(name = "All_traits_23_andMe",
               GWAS = res.t,prior_studies = mystudies,
               Z_matrices = "ZMatrices/",
               prior_shrinkage = 1,
               MR_threshold = 5e-08,
               MR_pruning_dist = 1000,
               MR_pruning_LD = 0.1,
               MR_shrinkage = 1,
               stop_after_prior = T,
               save_files = FALSE,
               verbose=T)

pdf("Coefficient_plot_bonf_23_and_ME.pdf")
coefficients_plot_bGWAS(res.bg2)
dev.off()




#idx=grep("Out-of-sample R-squared across all chromosomes",res.bg$log_info)
#r2.table=c(tratto,res.bg$log_info[idx],res.bg$log_info[idx+1])
#write(r2.table,ncol=3,file="r2.log",append=T,sep="\t")

#all.mr.betas=cbind(tratto,res.bg$significant_studies)
#write.table(all.mr.betas,file="Betas.log",append=T,sep="\t",row.names=F,quote=F,col.names=F)


res.corr=gwma.data.frame(bGWAS.obj=res.bg2,chi.df=1
                         ,no.se=T,correct.se=T,rs.lab="rsid"
                         ,z_filter=1,results.data.frame=res.t)
gwama.plot(gwaa.obj = res.corr,bGWAS.obj =res.bg,output.file = "23andMe2.png")

write.table(res.corr,file=paste0(tratto,"_corrected_GWAS_23_andme.tsv"),row.names=F,quote=F,sep="\t")
system(paste0("gzip -f ",tratto,"_corrected_GWAS.tsv"))

snpid	A1	A2	Zscore	N	P-value
res.corr$N=452300
ref.snp=fread("/opt/working/wilson/app_full_stuff/ldsc/eur_w_ld_chr/w_hm3.snplist",data.table=F)

raw_gen=res.corr[which(res.corr$rs%in%ref.snp$SNP),c("rs","alt","ref","observed_Z","N","p")]
names(raw_gen)=c("snpid",	"A1",	"A2",	"Zscore",	"N",	"P-value")
corr_gen=res.corr[which(res.corr$rs%in%ref.snp$SNP),c("rs","alt","ref","z_diff","N","corr.p")]
names(corr_gen)=c("snpid",	"A1",	"A2",	"Zscore",	"N",	"P-value")


write.table(raw_gen,file="Raw_results_ldhub.tsv",row.names=F,quote=F,sep="\t")
write.table(corr_gen,file="Corr_results_ldhub.tsv",row.names=F,quote=F,sep="\t")
system("zip Raw_results_ldhub.zip Raw_results_ldhub.tsv")
system("zip Corr_results_ldhub.zip Corr_results_ldhub.tsv")





study    estimate  std_error         T            P
1   EDUyears_2016_sumstat.txt  0.10611060 0.02521357  4.208471 5.877646e-05
2 MAGIC_FastingGlucose.txt.gz  0.12305130 0.04159334  2.958437 3.911454e-03
3     DIAGRAMv3.2013MAY07.zip -0.06777658 0.02966009 -2.285111 2.455476e-02
4 cardiogram_gwas_results.txt -0.10270807 0.05040973 -2.037465 4.441434e-02
[230] "## Out-of-sample R-squared across all chromosomes is 0.1865"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
[231] "## Out-of-sample squared correlation across all chromosome is 0.2026"    

