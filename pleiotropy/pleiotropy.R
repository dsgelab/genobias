##########
# Check if SNPs significantly associated with SEX in 23andMe are more likely to be pleiotropic than other SNPs and also check results in the GWAS catalog
#########
library(data.table)
library(plyr)
library(ggplot2)
library(LDlinkR) # devtools::install_github("CBIIT/LDlinkR")
library(gwascat)
library(ontologyIndex)


# Extract genobias/pleiotropy/uniqID.zip first
snps_all <- fread("genobias/pleiotropy/uniqID.txt", header=F)
snps_all$chr <- sapply(strsplit(snps_all$V1,":"),"[[",1)
snps_all$pos <- sapply(strsplit(snps_all$V1,":"),"[[",2)
snps_all$chr_pos <- paste0(snps_all$chr,"_",snps_all$pos)

rsid_translator <- fread("path/ukb.gwas_variants.tsv")
rsid_translator$chr_pos <- paste0(rsid_translator$chr,"_",rsid_translator$pos)

snps_allM <- merge(snps_all,rsid_translator,by="chr_pos")

snps_allMS <- snps_allM[,c("rsid")]
snps_allMS$X.traits <- 0
snps_allMS$X.domains <- 0

####### this?-
pleiotropy <- read.csv("genobias/pleiotropy/41588_2019_481_MOESM3_ESM.csv", header=TRUE, sep=";")
colnames(pleiotropy)[2] <- "rsid"
snps_allMS <- snps_allMS[!snps_allMS$rsid %in% pleiotropy$rsid,]
pleiotropyF <- rbind.fill(snps_allMS,pleiotropy)
colnames(pleiotropyF)[1:3] <- c("rsid","n_of_traits","n_of_domains")

write.table(pleiotropyF,file="pleiotropy_main_file.tsv", col.names=T, row.names=F, quote=F, sep="\t")



### read main pleiotropy file
pleiotropyF <- fread("genobias/pleiotropy/pleiotropy_main_file.tsv")
andme <- fread("23andMe_SexGWAS_AllSamples.assoc")

## Test if SNPs that are genome-wide significant are more enriched than expected
pleiotropyF$andme <- ifelse(pleiotropyF$rsid %in% andme$SNP[andme$P<0.00000005],1,0)
pleiotropyF$n_of_traits_cut <- ifelse(pleiotropyF$n_of_traits>=5,"5+",as.character(pleiotropyF$n_of_traits))

obs <- table(pleiotropyF$n_of_traits_cut[pleiotropyF$andme==1])
exp <- table(pleiotropyF$n_of_traits_cut)/length(pleiotropyF$n_of_traits_cut)

# Main test
chisq.test(obs, p=exp,correct=TRUE)$p.value



#### CREATE SUPPLEMENTARI TABLE 1 WITH THE RESULTS FROM GWAS of 23ANDME 
intervals <- fread("genobias/gwas_sex/GenomicRiskLoci.txt")
intervals$combine_locus <- paste0(intervals$chr,"_",intervals$start,"_",intervals$end)


andme <- fread("23andMe_SexGWAS_AllSamples.assoc")
andme30 <- fread("23andMe_SexGWAS_Age30.assoc")

## Keep common variants and imputation rsqr of 0.8
andmeQC <- andme[andme$effect_freq > 0.01 & andme$effect_freq < 0.99 & andme$imp_rsqr > 0.8,]
andme30QC <- andme30[andme30$effect_freq > 0.01 & andme30$effect_freq < 0.99 & andme30$imp_rsqr > 0.8,]

## Read additional files
directly_genotypes <- fread("genobias/gwas_sex/direct_genotyped_SNPs_results.tsv")

intervalsmm <- merge(intervals,directly_genotypes,by="combine_locus",all.x=T)
intervalsmm2 <- merge(intervalsmm,andme,by.x="rsID",by.y="SNP")

write.table(intervalsmm2[,c("combine_locus","rsID","A1","A2","Effect","se","P","varID")],file="imputed_SNPs_results.tsv",row.names=F, col.names=T, quote=F, sep="\t")


#################################################################
## Check loci associated with other traits in GWAS catalog ######
#################################################################


intervals <- fread("genobias/gwas_sex/GenomicRiskLoci.txt")

in_pleio_ana <- NULL
for (i in 1:nrow(intervals))
{

	andmeS <- andme[andme$position_b37 < as.numeric(intervals[i,"end"]) & andme$position_b37 > as.numeric(intervals[i,"start"]) & andme$chrom == as.numeric(intervals[i,"chr"]) & andme$SNP %in%  pleiotropyF$rsid,]

	if (nrow(andmeS)>0)
	{
	smallest_P <- andmeS[which.min(andmeS$P),]
	r2out <- tryCatch(expr={LDpair(var1 = smallest_P$SNP, var2 = as.character(intervals[i,"rsID"]), pop="CEU", token="dfedfc272aef")$r2},
    error = function(c) {NA})

	in_pleio_ana <- rbind(in_pleio_ana,cbind(smallest_P,r2out))
	}

	print(i)
}

## Keep only if r2 > 0.2
in_pleio_anaS <- in_pleio_ana[in_pleio_ana$r2out > 0.2,]
pleiotropyF$gwloci <- ifelse(pleiotropyF$rsid %in% in_pleio_anaS$SNP,1,0)

one_plus <- ifelse(pleiotropyF$X.traits >0,1,0)
five_plus <- ifelse(pleiotropyF$X.traits >4,1,0)
table(pleiotropyF$gwloci,one_plus)
table(pleiotropyF$gwloci,five_plus)


## Check most associated traits using the GWAS catalog
setwd("path")
#LDproxy_batch(intervals$rsID, pop = "CEU", r2d = "r2", token = "dfedfc272aef", append = TRUE)

################### query_snp_list
all_proxies <- fread("genobias/pleiotropy/combined_query_snp_list.txt")
all_proxiesS <- all_proxies[all_proxies$R2 > 0.2,]


ebicat37_latest <- makeCurrentGwascat(table.url =
"http://www.ebi.ac.uk/gwas/api/search/downloads/alternative",
fixNonASCII = TRUE, genome="GRCh37",
withOnt = TRUE)


RESPROXY <- NULL
for (snp in unique(all_proxiesS$query_snp))
{
	intes <- intersect(getRsids(ebicat37_latest) , all_proxiesS$RS_Number[all_proxiesS$query_snp==snp])
	if (length(intes) > 0)
	{ttout <-  ebicat37_latest[ intes]
	traits <- ttout@elementMetadata@listData$MAPPED_TRAIT
	pvals <- ttout@elementMetadata@listData['P-VALUE'][[1]]
	snpproxy <- getRsids(ttout)[!duplicated(traits) & pvals < 0.00000005]
	if(length(snpproxy)>0)
	{traits <- traits[!duplicated(traits) & pvals < 0.00000005]
	RESPROXY <- rbind(RESPROXY,cbind(snp,snpproxy,traits))}}
}


to_exp <- data.frame(aggregate(list(RESPROXY[,2],RESPROXY[,3]), list(RESPROXY[,1]), function(x){paste(x,collapse = ",")}))
colnames(to_exp) <- c("rsID","rsid_proxy","traits")
intervals <- fread("genobias/gwas_sex/GenomicRiskLoci.txt")
intervals$combine_locus <- paste0(intervals$chr,"_",intervals$start,"_",intervals$end)

to_exp <- merge(to_exp,intervals[,c("rsID","combine_locus")],by.x="rsID")

write.table(to_exp,file="gwas_catalog_look_up.tsv", row.names=F, quote=F, sep="\t")


