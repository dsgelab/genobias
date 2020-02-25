linbrary(ggplot2)
library(ggrepel)

traitmap <- data.frame(rbind(
c("ADHD_2017","ADHD"),
c("AgeFirstBirth_Female","Age at first birth (Females)"),
c("AgeFirstBirth_Male","Age at first birth (Males)"),
c("BMI","BMI"),
c("CAD","Coronary Artery Disease"),
c("EA","Education attainment"),
c("IBD","Inflammatory Bowel Disease"),
c("LDL","LDL Cholesterol"),
c("MDD2018_ex23andMe.sumstats","Major depressive disorder"),
c("Neuroticism_Full","Neuroticism"),
c("NumberChildrenEverBorn_Female","Number of children (Females)"),
c("NumberChildrenEverBorn_Male","Number of children (Males)"),
c("SCZ2","Schizophrenia"),
c("SWB_Full","Subjective well-being"),
c("T2D_2018","Type 2 diabetes"),
c("TG","Triglycerides"),
c("UKB.self_rated_health","Self-rated health"),
c("age_at_menarche","Age at menarche"),
c("age_at_menopauze","Age at menopause"),
c("alcohol_clarke","Alcohol use"),
c("anorexia_2019","Anorexia"),
c("anxiety","Anxiety"),
c("autism_2017.ipsych.pgc","Autism"),
c("bipolar","Bipolar disorder"),
c("birth_weight","Birth weight"),
c("breast_cancer","Breast Cancer"),
c("cannabis_ever_2018.no23andMe","Cannabis use"),
c("finger2d4d.average.sumstats","2D:4D digit ratio"),
c("fasting_insuline","Fasting Insuline"),
c("height_combined","Height"),
c("loneliness.rapid_UKB.sumstats","Loneliness"),
c("openness","Openness to experience"),
c("prostate_cancer","Prostate Cancer"),
c("risk_taking_PC1","Risk behaviour"),
c("smoking_ever_vs_never","Smoking: ever smoking"),
c("smoking_cigs_per_day","Smoking: cigs/day"),
c("whr_females","Waist-to-hip ratio (Females)"),
c("whr_males","Waist-to-hip ratio (Males)")))



### Function to extract genetic correlaition
get_h2_se <- function(x)
{
	tt <- system( paste0('grep "/stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/" ', x), intern = T)[5]


	tt2 <- strsplit(tt," ")
	
	if (tt2[[1]][4] == "")
	{rg <- as.numeric(tt2[[1]][5])
	se <- as.numeric(tt2[[1]][7])
	p <- as.numeric(tt2[[1]][11])}
	else
	{rg <- as.numeric(tt2[[1]][4])
	se <- as.numeric(tt2[[1]][6])
	p <- as.numeric(tt2[[1]][9])}

	return(c(rg,se,p))
}




#### Run genetic correlations ####

files <- list.files("/stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/genetic_correlation/",pattern="*.sumstats.gz")

namesv <- c("UKBB_AgeAdj_sex_23andMe_coded_Imputed_edit.sumstats.gz","23andMe_SexGWAS_30younger_edited.sumstats.gz" ,"23andMe_SexGWAS_AllSamples_edited.sumstats.gz")

rgs <- NULL

for (names in namesv)
{
	for (f in files)
	{
	system(paste0("/stanley/genetics/analysis/software/aganna/ldsc-master/ldsc.py --rg /stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/",names,",/stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/genetic_correlation/",f," --ref-ld-chr /stanley/genetics/analysis/software/aganna/ldsc-master/eur_w_ld_chr/ --w-ld-chr /stanley/genetics/analysis/software/aganna/ldsc-master/eur_w_ld_chr/ --out /stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/genetic_correlation/",gsub(".sumstats.gz","",names),"?",gsub(".sumstats.gz","",f),"_rg"))

	temp <- get_h2_se(paste0("/stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/genetic_correlation/",gsub(".sumstats.gz","",names),"?",gsub(".sumstats.gz","",f),"_rg.log"))

	rgs <- rbind(rgs,c(names,f,temp))

	}
	print(names)
}


df <- data.frame(rgs)
df$matchname <- gsub(".LDSC.sumstats.gz","",df$X2)
dfm <- merge(df,traitmap,by.x="matchname", by.y="X1")
dfm$group <- NA
dfm$group[dfm$X1=="UKBB_AgeAdj_sex_23andMe_coded_Imputed_edit.sumstats.gz"] <- "UKBB"
dfm$group[dfm$X1=="23andMe_SexGWAS_AllSamples_edited.sumstats.gz"] <- "23andMe"
dfm$group[dfm$X1=="23andMe_SexGWAS_30younger_edited.sumstats.gz"] <- "23andMe_lt30"

dfm_toexp <- dfm[,colnames(dfm) %in% c("matchname","X3","X4","X5","X2.y","group")]
colnames(dfm_toexp) <- c("matchname","rg","se","pval","trait","group")

write.csv(dfm_toexp,"/stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/genetic_correlation/SexGWAS2019_GeneticCorrelations_final.csv", row.names=F, quote=F)



### Compare 23andme and UK Biobank 
d <- read.csv("/stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/genetic_correlation/SexGWAS2019_GeneticCorrelations_final.csv")

mm <- merge(d,traitmap,by.x="matchname", by.y="X1")
mm <- mm[!mm$group %in% c("23andMe_lt30"),]
mm$ymin <- mm$rg - 1.96*mm$se
mm$ymax <- mm$rg + 1.96*mm$se

# How many significant in Uk Biobank, how many in 23andMe and in Uk Biobank
dim(mm[mm$group=="UKBB" & mm$pval < (0.05/38),])
dim(mm[mm$group=="23andMe" & mm$pval < (0.05/38),])

mmr <- data.frame(cbind(mm[mm$group=="UKBB",],mm[mm$group=="23andMe",]),stringsAsFactors=F)
mmr$pvalind <- factor(ifelse(mmr$pval < (0.05/nrow(mm)) | mmr$pval.1 < (0.05/nrow(mm)),1,0))
mmr$X2[mmr$pvalind==0] <- NA

##
pdf("/stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/genetic_correlation/rg_ukb_23andMe.pdf",width=6,height=5)
ggplot(aes(x=rg,y=rg.1),data=mmr) + geom_point(aes(col=pvalind,alpha=pvalind)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + ylab("Genetic correlation with sex in 23andMe") + xlab("Genetic correlation with sex in UK Biobank")  + geom_abline(intercept =0 , slope = 1, size=0.1) + geom_errorbarh(aes(xmin = ymin,xmax = ymax, col=pvalind), alpha=0.5) + geom_errorbar(aes(ymin = ymin.1,ymax = ymax.1, col=pvalind), alpha=0.5) + scale_colour_manual(values=c("grey","steelblue4"), guide=FALSE) + scale_alpha_manual(values=c(0.3,1), guide=FALSE) + xlim(-0.5,0.5) + ylim(-0.5,0.5) + geom_text_repel(aes(label=X2),size=3) + geom_vline(xintercept =0, size=0.1) + geom_hline(yintercept =0, size=0.1)
dev.off()


### genetic correlation between UK Biobank and 23andMe 

/stanley/genetics/analysis/software/aganna/ldsc-master/ldsc.py --rg /stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/23andMe_SexGWAS_AllSamples_edited.sumstats.gz,/stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/UKBB_AgeAdj_sex_23andMe_coded_Imputed_edit.sumstats.gz --ref-ld-chr /stanley/genetics/analysis/software/aganna/ldsc-master/eur_w_ld_chr/ --w-ld-chr /stanley/genetics/analysis/software/aganna/ldsc-master/eur_w_ld_chr/ --out /stanley/genetics/analysis/ukbb/aganna/uk_bio/bias/gwas_of_sex/23andMe_vs_ukBiobank

