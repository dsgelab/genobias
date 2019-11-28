library(data.table)
library(qqman)
library(gtools)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg19)


### Function to test the differences between <30 and all
betadiff <- function(beta1,beta2,se1,se2,N1,N2,cti)
{
	z1 <-  beta1/se1
	z2 <- beta2/se2

	w1 <- N1
	w2 <- N2

	sexdifz <- (1/sqrt(w1)*z1 - 1/sqrt(w2)*z2) / sqrt(1/w1 + 1/w2 - 2*sqrt(1/w1*1/w2)*cti)
	pvaldif <- 2*pnorm(-abs(sexdifz))
	return(list(sexdifz,pvaldif))
}


### Function for the manhattan plot
ggman_bw <- function(gwas, bp = NA, chrom = NA, pvalue = NA, intervals=NA, ymax=NA, minval=0.01) 
	{
		dfm <- as.data.frame(gwas)


		posmin <- tapply(dfm[,bp],dfm[,chrom], min)
    	posmax <- tapply(dfm[,bp],dfm[,chrom], max)
    	posshift <- head(c(0,cumsum(as.numeric(posmax))),-1)

	    names(posshift) <- names(posmin)
	    for (k in unique(dfm[,chrom]))
	    {
	        dfm$pos_new[dfm[,chrom]==k] <-  dfm[,bp][dfm[,chrom]==k] + posshift[names(posshift) == k]
	    }

	  	dfmsplit <- split(dfm, dfm[,chrom])
	    xbreaks <- sapply(dfmsplit,function(x) x$pos_new[length(x$pos_new)/2])

	    dfm$marker <- -log10(dfm[,pvalue])
	    df_manhattan <- dfm[dfm$marker > -log10(minval),]


	    ymax <- ifelse(max(df_manhattan$marker) < -log10(5e-8), -log10(5e-8), max(df_manhattan$marker)) + 0.2


		chrtable <- data.frame(table(df_manhattan[,chrom]))
		chrtable$Var1 <- as.character(chrtable$Var1)
		chrtable <- chrtable[mixedorder(chrtable$Var1),]
		oddchrom <- as.character(chrtable$Var1[seq(1,nrow(chrtable),2)])
		df_manhattan$chrom_alt <- replace(df_manhattan[,chrom], df_manhattan[,chrom] %in% oddchrom, 0)
		df_manhattan$chrom_alt <- replace(df_manhattan$chrom_alt, df_manhattan$chrom_alt != 0,1)
		df_manhattan$chrom_altA <- ifelse(df_manhattan$chrom_alt=="1",1,2)
		if (is.na(ymax)) {ymax <- max(-log10(df_manhattan[,pvalue])) + 1}

		df_manhattan$inint <- 1
		df_manhattan$shape <- 16
		df_manhattan$fill <- 0
		for (i in 1:nrow(intervals))
		{
			df_manhattan$inint[df_manhattan[,chrom]==intervals$chr[i] & df_manhattan[,bp]<=intervals$end[i] & df_manhattan[,bp]>=intervals$start[i]] <- 1
			df_manhattan$chrom_altA[df_manhattan[,chrom]==intervals$chr[i] & df_manhattan[,bp]<=intervals$end[i] & df_manhattan[,bp]>=intervals$start[i]] <- 1
			df_manhattan$fill[df_manhattan[,chrom]==intervals$chr[i] & df_manhattan[,bp]<=intervals$end[i] & df_manhattan[,bp]>=intervals$start[i]] <- 3
			df_manhattan$shape[df_manhattan[,chrom]==intervals$chr[i] & df_manhattan[,bp]<=intervals$end[i] & df_manhattan[,bp]>=intervals$start[i]] <- 16

			ind_min = intersect(which(df_manhattan[,pvalue] == min(df_manhattan[,pvalue][df_manhattan[,chrom]==intervals$chr[i] & df_manhattan[,bp]<=intervals$end[i] & df_manhattan[,bp]>=intervals$start[i]])), seq(df_manhattan[,pvalue])[df_manhattan[,chrom]==intervals$chr[i] & df_manhattan[,bp]<=intervals$end[i] & df_manhattan[,bp]>=intervals$start[i]])

			df_manhattan$shape[ind_min] <- 23

			df_manhattan$inint[ind_min] <- 3
			df_manhattan <- rbind(df_manhattan,df_manhattan[ind_min,])

		}
		p1 <- ggplot(df_manhattan, aes(x = pos_new,y = marker)) +
	                geom_point(aes(size=as.factor(inint),
	                	shape=shape,
	                	fill=as.factor(fill),
	                	alpha=as.factor(chrom_altA),
	                	colour=as.factor(chrom))) + 
	                scale_shape_identity() +
	                scale_x_continuous(breaks = xbreaks, labels = names(xbreaks), expand=c(0,0)) +
	                scale_y_continuous(expand = c(0,0), limits = c(-log10(minval),ymax), breaks=c(2,4,6,8,seq(10,ymax,4)),labels=c(2,4,6,8,seq(10,ymax,4)))+
	                expand_limits(x = 22.3, y=(ymax+3)) +
	                guides(colour = FALSE,alpha=FALSE, size=FALSE, fill=FALSE) +
	                labs(x = "chromosome", y = "-log10(P value)") + 
	                theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_rect(fill = "transparent"),  plot.background = element_rect(fill = "transparent", color = NA),axis.line = element_line(colour = "red")) +
	            	geom_hline(aes(yintercept= -log10(0.00000005)),colour = "red", lwd=0.6, linetype = 5)  + 
	            	scale_colour_manual(values = rep(c("grey67","grey22"),11)) + 
	            	scale_alpha_manual(values=c(1,1)) +
	            	scale_size_manual(values=c(0.4,3)) + 
	            	scale_fill_manual(values=c("0"="white","3"="red"))
	}



###############
### 23andMe ###
###############

######### STEP 1. Import data and check QC #######

andme <- fread("/stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/23andMe_SexGWAS_AllSamples.assoc")
andme30 <- fread("/stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/23andMe_SexGWAS_Age30.assoc")
andme_direct <- fread("/stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/23andMe_Sex_AllSamples_Geno_direct_geno.assoc")

## Keep common variants and imputation rsqr of 0.8
andmeQC <- andme[andme$effect_freq > 0.01 & andme$effect_freq < 0.99 & andme$imp_rsqr > 0.8,]

andme30QC <- andme30[andme30$effect_freq > 0.01 & andme30$effect_freq < 0.99 & andme30$imp_rsqr > 0.8,]


dim(andmeQC[andmeQC$P < 5e-8,])


### Write files to pass to FUMA
write.table(andmeQC[,c("SNP","A1","A2","N","Effect","se","P")], file="/stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/23andMe_SexGWAS_AllSamples_edited.tsv", row.names=F, quote=F, sep="\t")
write.table(andme30QC[,c("SNP","A1","A2","N","Effect","se","P")], file="/stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/23andMe_SexGWAS_30younger_edited.tsv", row.names=F, quote=F, sep="\t")



######### STEP 2. ADDITIONAL QC ON DIRECTLY GENOTYPED SNPS #######
 intervals <- fread("/stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/FUMA/GenomicRiskLoci.txt")


# Extract directly genotyped SNP with lowest P-value for each loci and add sequence 50 and 250 bp upstream and downstream 
direct_geno_per_locus <- NULL
for (i in 1:nrow(intervals))
{
  andme_directS <- andme_direct[andme_direct$position_b37 < as.numeric(intervals[i,"end"]) & andme_direct$position_b37 > as.numeric(intervals[i,"start"]) & andme_direct$chrom == as.numeric(intervals[i,"chr"]),]
  smallest_P <- andme_directS[which.min(andme_directS$pvalue),]
  if(nrow(smallest_P)>0)
  {
    #seqout250 <- toString(getSeq(Hsapiens, paste0("chr",as.character(smallest_P$chrom)), start = smallest_P$position_b37 - 250, end = smallest_P$position_b37 + 250))
    seqout50 <- toString(getSeq(Hsapiens, paste0("chr",as.character(smallest_P$chrom)), start = smallest_P$position_b37 - 50, end = smallest_P$position_b37 + 50))
    direct_geno_per_locus <- rbind(direct_geno_per_locus,data.frame(locus_start=as.numeric(intervals[i,"start"]),locus_end=as.numeric(intervals[i,"end"]),seqout50=seqout50,smallest_P))
  }
}

# Keep only if P-value < 10x-8
direct_geno_per_locus <- direct_geno_per_locus[direct_geno_per_locus$pvalue < 0.00000005,]

# Write sequences for BLAST
TO_EXP <- NULL
for(k in 1:nrow(direct_geno_per_locus))
{TO_EXP <- rbind(TO_EXP,rbind(paste0(">",direct_geno_per_locus$varID[k]),as.character(direct_geno_per_locus$seqout50[k])))}
write.table(TO_EXP,file="/stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/blat_50_seq.txt",col.names=F,row.names=F,quote=F)


## Read in results and add flags to the final results
seq50_res <- read.csv("/stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/blat_50_seq_res.csv", sep=";")
colnames(seq50_res)[1] <- "QUERY"
seq50_res$similar <- as.numeric(gsub(" %","",seq50_res$IDENTITY))


snp_exclude <- unique(seq50_res$QUERY[seq50_res$CHROM %in% c("chrX","chrY") & seq50_res$similar > 95])

direct_geno_per_locus$flag_homology <- ifelse(direct_geno_per_locus$varID %in% snp_exclude,1,0)
direct_geno_per_locus$flag_maf <- ifelse(direct_geno_per_locus$geno_freq_a < 0.05 | direct_geno_per_locus$geno_freq_a > 0.95,1,0)
direct_geno_per_locus$flag_hwe <- ifelse(direct_geno_per_locus$geno_hwe_p < 0.000001 ,1,0)
direct_geno_per_locus$flag_callrate <- ifelse(direct_geno_per_locus$gt.rate < 0.98 ,1,0)

dim(direct_geno_per_locus[direct_geno_per_locus$flag_maf==1 | direct_geno_per_locus$flag_hwe==1 | direct_geno_per_locus$flag_callrate==1,])

dim(direct_geno_per_locus[direct_geno_per_locus$flag_maf==0 & direct_geno_per_locus$flag_hwe==0 & direct_geno_per_locus$flag_callrate==0 & direct_geno_per_locus$flag_homology==0,])


######### STEP 3. PLOT MANHATTAN PLOT #######

p_fin <- ggman_bw(gwas=andmeQC,pvalue="P",chrom="chrom",bp="position_b37",intervals=data.frame(chr=intervals$chr,start=intervals$start,end=intervals$end),ymax=-log10(0.0000000005))


pdf("/stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/23andme_mm_grey.pdf", width=11, height=5)
p_fin
dev.off()


######### STEP 4. CHECK FTO LOCUS #######

andme[andme$SNP=="rs10468280",]
wget https://www.dropbox.com/s/gulqmkaighh8w3b/21001_irnt.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 21001_irnt.gwas.imputed_v3.both_sexes.tsv.gz
grep 16:53827479 21001_irnt.gwas.imputed_v3.both_sexes.tsv



######### STEP 5. COMPARE SNPs IN INDIVIDUALS < 30 vs ALL #######

## Compare SNPs all ages and younger than 30
andmeall <- merge(andmeQC,andme30QC,by="SNP")
andmeall <- andmeall[andmeall$SNP %in% intervals$rsID,]

andmeall$pval_all<- -log10(andmeall$P.x)
andmeall$pval_30 <- -log10(andmeall$P.y)

andmeall$betaallmin <- andmeall$Effect.x - 1.96*andmeall$se.x
andmeall$betaallmax <- andmeall$Effect.x + 1.96*andmeall$se.x
andmeall$beta30min <- andmeall$Effect.y - 1.96*andmeall$se.y
andmeall$beta30max <- andmeall$Effect.y + 1.96*andmeall$se.y

andmeall$beta30min[is.nan(andmeall$beta30min)] <- 0
andmeall$beta30max[is.nan(andmeall$beta30max)] <- 0

## Statistics

table(sign(andmeall$Effect.x),sign(andmeall$Effect.y))
cor(andmeall$Effect.x,andmeall$Effect.y)

testeffectdiff <- betadiff(andmeall$Effect.x,andmeall$Effect.y,andmeall$se.x,andmeall$se.y,cti=0.3648,2461431,320366)
min(testeffectdiff[[2]])


pdf("/stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/snp_replication_andme_andmelt30.pdf", width=5, height=5)
ggplot(aes(y=Effect.x,x=Effect.y,ymin=betaallmin,ymax=betaallmax,xmin=beta30min,xmax=beta30max),data=andmeall) + geom_point(color="red",size=3)  +  geom_errorbar(aes(ymin = betaallmin,ymax = betaallmax),size=0.001) + geom_errorbarh(aes(xmin = beta30min,xmax = beta30max),size=0.001) + theme_bw() + ylab("Coefficient for all individuals") + xlab("Coefficient for individuals < 30 years old") + geom_abline(intercept=0, slope=1)  
dev.off()


######### STEP 5. RUN LDSCORE #######


## Heritability for all
/stanley/genetics/analysis/software/aganna/ldsc-master/munge_sumstats.py --sumstats /stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/23andMe_SexGWAS_AllSamples_edited.tsv --N-col N --out /stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/23andMe_SexGWAS_AllSamples_edited --snp SNP --a1 A1 --a2 A2 --p P --signed-sumstats Effect,0 --merge-alleles /stanley/genetics/analysis/software/aganna/ldsc-master/w_hm3.snplist

/stanley/genetics/analysis/software/aganna/ldsc-master/ldsc.py --h2 /stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/23andMe_SexGWAS_AllSamples_edited.sumstats.gz \
--ref-ld-chr /stanley/genetics/analysis/software/aganna/ldsc-master/1000G_Phase3_baselineLD_ldscores/baselineLD. \
--w-ld-chr /stanley/genetics/analysis/software/aganna/ldsc-master/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--overlap-annot \
--frqfile-chr /stanley/genetics/analysis/software/aganna/ldsc-master/1000G_Phase3_frq/1000G.EUR.QC. \
--print-coefficients \
--out /stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/23andMe_SexGWAS_AllSamples_edited


## Heritability for younger than 30
/stanley/genetics/analysis/software/aganna/ldsc-master/munge_sumstats.py --sumstats /stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/23andMe_SexGWAS_30younger_edited.tsv --N-col N --out /stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/23andMe_SexGWAS_30younger_edited --snp SNP --a1 A1 --a2 A2 --p P --signed-sumstats Effect,0 --merge-alleles /stanley/genetics/analysis/software/aganna/ldsc-master/w_hm3.snplist

/stanley/genetics/analysis/software/aganna/ldsc-master/ldsc.py --h2 /stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/23andMe_SexGWAS_AllSamples_edited.sumstats.gz \
--ref-ld-chr /stanley/genetics/analysis/software/aganna/ldsc-master/1000G_Phase3_baselineLD_ldscores/baselineLD. \
--w-ld-chr /stanley/genetics/analysis/software/aganna/ldsc-master/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--overlap-annot \
--frqfile-chr /stanley/genetics/analysis/software/aganna/ldsc-master/1000G_Phase3_frq/1000G.EUR.QC. \
--print-coefficients \
--out /stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/23andMe_SexGWAS_AllSamples_edited


## Genetic correlation between all and younger than 30
/stanley/genetics/analysis/software/aganna/ldsc-master/ldsc.py --rg /stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/23andMe_SexGWAS_AllSamples_edited.sumstats.gz,/stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/23andMe_SexGWAS_30younger_edited.sumstats.gz --ref-ld-chr /stanley/genetics/analysis/software/aganna/ldsc-master/eur_w_ld_chr/ --w-ld-chr /stanley/genetics/analysis/software/aganna/ldsc-master/eur_w_ld_chr/ --out /stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/23andMe_SexGWAS_AllSamples_edited_vs_23andMe_SexGWAS_30younger_edited






#################
### UK BIOBANK ##
#################


ukbb <- fread("/stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/UKBB_AgeAdj_sex_23andMe_coded_Imputed")
dim(ukbb[ukbb$P_BOLT_LMM_INF < 5e-8,])


## save for FUM
write.table(ukbb[,c("SNP","ALLELE1","ALLELE0","BETA","P_BOLT_LMM_INF")], file="/stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/UKBB_AgeAdj_sex_23andMe_coded_Imputed_edit.tsv", row.names=F, quote=F, sep="\t")


## Munge
/stanley/genetics/analysis/software/aganna/ldsc-master/munge_sumstats.py --sumstats /stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/UKBB_AgeAdj_sex_23andMe_coded_Imputed_edit.tsv --N 452302 --out /stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/UKBB_AgeAdj_sex_23andMe_coded_Imputed_edit --snp rsid --a1 ALLELE1 --a2 ALLELE0 --p P_BOLT_LMM_INF --signed-sumstats BETA,0 --merge-alleles /stanley/genetics/analysis/software/aganna/ldsc-master/w_hm3.snplist


#h2
/stanley/genetics/analysis/software/aganna/ldsc-master/ldsc.py --h2 /stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/UKBB_AgeAdj_sex_23andMe_coded_Imputed_edit.sumstats.gz \
--ref-ld-chr /stanley/genetics/analysis/software/aganna/ldsc-master/1000G_Phase3_baselineLD_ldscores/baselineLD. \
--w-ld-chr /stanley/genetics/analysis/software/aganna/ldsc-master/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--overlap-annot \
--frqfile-chr /stanley/genetics/analysis/software/aganna/ldsc-master/1000G_Phase3_frq/1000G.EUR.QC. \
--print-coefficients \
--out /stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/UKBB_AgeAdj_sex_23andMe_coded_Imputed_edit




###########
### BBJ ###
###########


bbj1_cov <- fread("/stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/BBJ_covD_an1_230619.tsv")


bbj1_cov$P <- 2*pnorm(-abs(bbj1_cov$beta/bbj1_cov$se))
dim(bbj1_cov[bbj1_cov$P < 5e-8,])

write.table(bbj1_cov[,c("rsid","reference_allele","effect_allele","n","beta","se","P")], file="/stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/BBJ_covD_an1_230619_edited.tsv", row.names=F, quote=F, sep="\t")

/stanley/genetics/analysis/software/aganna/ldsc-master/munge_sumstats.py --sumstats /stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/BBJ_covD_an1_230619_edited.tsv --N-col n --out /stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/BBJ_covD_an1_230619_edited --snp rsid --a1 effect_allele --a2 reference_allele --p P --signed-sumstats beta,0 --merge-alleles /stanley/genetics/analysis/software/aganna/ldsc-master/w_hm3.snplist


/stanley/genetics/analysis/software/aganna/ldsc-master/ldsc.py --h2 /stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/BBJ_covD_an1_230619_edited.sumstats.gz \
--ref-ld-chr /stanley/genetics/analysis/software/aganna/ldsc-master/1000G_Phase3_baselineLD_ldscores/baselineLD. \
--w-ld-chr /stanley/genetics/analysis/software/aganna/ldsc-master/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--overlap-annot \
--frqfile-chr /stanley/genetics/analysis/software/aganna/ldsc-master/1000G_Phase3_frq/1000G.EUR.QC. \
--print-coefficients \
--out /stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/BBJ_covD_an1_230619_edited


##############
### iPSYCH ###
##############


pyshc_cov <- fread("/stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/iPSYCHwPsychCov_an1_09072019.tsv")


pyshc_cov$P <- 2*pnorm(-abs(pyshc_cov$beta/pyshc_cov$se))
dim(pyshc_cov[pyshc_cov$P < 5e-8,])


write.table(pyshc_cov[,c("rsid","reference_allele","effect_allele","n","beta","se","P")], file="/stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/iPSYCHwPsychCov_an1_09072019_edited.tsv", row.names=F, quote=F, sep="\t")

/stanley/genetics/analysis/software/aganna/ldsc-master/munge_sumstats.py --sumstats /stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/iPSYCHwPsychCov_an1_09072019_edited.tsv --N-col n --out /stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/iPSYCHwPsychCov_an1_09072019_edited --snp rsid --a1 effect_allele --a2 reference_allele --p P --signed-sumstats beta,0 --merge-alleles /stanley/genetics/analysis/software/aganna/ldsc-master/w_hm3.snplist


/stanley/genetics/analysis/software/aganna/ldsc-master/ldsc.py --h2 /stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/iPSYCHwPsychCov_an1_09072019_edited.sumstats.gz \
--ref-ld-chr /stanley/genetics/analysis/software/aganna/ldsc-master/1000G_Phase3_baselineLD_ldscores/baselineLD. \
--w-ld-chr /stanley/genetics/analysis/software/aganna/ldsc-master/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--overlap-annot \
--frqfile-chr /stanley/genetics/analysis/software/aganna/ldsc-master/1000G_Phase3_frq/1000G.EUR.QC. \
--print-coefficients \
--out /stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/iPSYCHwPsychCov_an1_09072019_edited



##############
### Finngen ##
##############

snpid <- fread("/stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/r3_rsids.tsv")

finngen_an1 <- fread("/stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/Finngen_r3_SEX_IMPUTED.pheweb")
finngen_an1$V2 <- paste0(finngen_an1$"#chrom","_",finngen_an1$pos)


finngen_an1_m <- merge(finngen_an1,snpid,by="V2")
finngen_an1_m <- finngen_an1_m[!duplicated(finngen_an1_m$V1),]
finngen_an1_m <- finngen_an1_m[finngen_an1_m$"#chrom" != "X",]


dim(finngen_an1_m[finngen_an1_m$pval < 5e-8,])

write.table(finngen_an1_m[,c("V1","ref","alt","beta","sebeta","pval")], file="/stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/Finngen_r3_SEX_IMPUTED.pheweb_edited.tsv", row.names=F, quote=F, sep="\t")


## Munge 
/stanley/genetics/analysis/software/aganna/ldsc-master/munge_sumstats.py --sumstats /stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/Finngen_r3_SEX_IMPUTED.pheweb_edited.tsv --N 135638 --out /stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/Finngen_r3_SEX_IMPUTED.pheweb_edited --snp V1 --a1 alt --a2 ref --p pval --signed-sumstats beta,0 --merge-alleles /stanley/genetics/analysis/software/aganna/ldsc-master/w_hm3.snplist


#h2
/stanley/genetics/analysis/software/aganna/ldsc-master/ldsc.py --h2 /stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/Finngen_r3_SEX_IMPUTED.pheweb_edited.sumstats.gz \
--ref-ld-chr /stanley/genetics/analysis/software/aganna/ldsc-master/1000G_Phase3_baselineLD_ldscores/baselineLD. \
--w-ld-chr /stanley/genetics/analysis/software/aganna/ldsc-master/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--overlap-annot \
--frqfile-chr /stanley/genetics/analysis/software/aganna/ldsc-master/1000G_Phase3_frq/1000G.EUR.QC. \
--print-coefficients \
--out /stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/Finngen_r3_SEX_IMPUTED.pheweb_edited



#########################################################################################
#### PLOT OF HERITABILITIES - STILL UNCLEAR WHICH SCALE TO USE, UP FOR DISCUSSION #######
#########################################################################################

liabscale <- function(coef,se,prev,samp)
{coefe <- coef*prev*(1-prev)/(dnorm(qnorm(prev))^2)
see <- se*prev*(1-prev)/(dnorm(qnorm(prev))^2)
pval <- 2*pnorm(-abs(coefe/see))
coefe2 <- coef*(prev*(1-prev)^2/(samp*(1-samp)*(dnorm(qnorm(prev))^2)))
see2 <- se*(prev*(1-prev)^2/(samp*(1-samp)*(dnorm(qnorm(prev))^2)))
pval2 <- 2*pnorm(-abs(coefe2/see2))
return(list(coefe,see,pval,coefe2,see2,pval2))}


rg <- c(0.0092,-0.0074,0.0087,0.0145,0.0192)
se <- c(0.005,0.0109,0.0054,0.0019,0.0008)
prev <- c(0.46,0.47,0.56,0.54,0.53)
samp <- c(rep(0.5,5))
pval <- 2*pnorm(-abs(rg/se))

rg_liab <- liabscale(rg,se,prev,samp)[[1]]
se_liab <- liabscale(rg,se,prev,samp)[[2]]
pval_liab <- liabscale(rg,se,prev,samp)[[3]]


rg_liab2 <- liabscale(rg,se,prev,samp)[[4]]
se_liab2 <- liabscale(rg,se,prev,samp)[[5]]
pval_liab2 <- liabscale(rg,se,prev,samp)[[6]]


label <- c("Biobank Japan","iPSYCH","Finngen","UK Biobank","23andMe")
labelsig <- c("","","","***","***")
ymin <- rg_liab2 - 1.96*se_liab2
ymax <- rg_liab2 + 1.96*se_liab2
df <- data.frame(rg_liab2,se_liab2,label,ymin,ymax,pval_liab2,labelsig)
df$rg_liab[df$label=="iPSYCH"] <- 0.0001

pdf("/stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/heritability_plot.pdf", width=5,height=3)
ggplot(aes(y=rg_liab,x=label), data=df) + geom_bar(stat="identity", aes(fill=labelsig)) + theme_bw() + ylab("SNP-heritability for sex") + coord_cartesian(ylim = c(0,0.06)) + xlab("") +  geom_text(aes(label=paste0("P==",gsub('e-0*', ' %*% 10^-', prettyNum(df$pval_liab, digits=2)))), parse=TRUE,vjust=-0.4, size=3, hjust=+0.4) + scale_fill_manual(values=c("blue","red")) + theme(legend.position = "none") 
dev.off()



