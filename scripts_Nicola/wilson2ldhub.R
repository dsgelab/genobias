wilson2ldhub=function(file="simul/gwas/edu_sex.tsv.gz",out="simul/ldhub/edu_sex.sumstat"){
 require(data.table)
 data=fread(file,data.table=F)
 snp.list=fread("/opt/working/wilson/app_full_stuff/ldsc/eur_w_ld_chr/w_hm3.snplist")
 data=data[which(data$rsid%in%snp.list$SNP),]
 data$Zscore=data$beta1/data$se
 data=data[,c("rsid","a1","a0","Zscore","n","p")]
 names(data)=c("snpid","A1","A2","Zscore","N","P-value")
 write.table(data,file=out,sep="\t",row.names=F,quote=F)
 
 system(paste0("munge_sumstats.py --sumstats ",out,"  --merge-alleles /opt/working/wilson/app_full_stuff/ldsc/eur_w_ld_chr/w_hm3.snplist --out ",out,"_mung"))

}
file.list=system("ls simul/gwas/*",intern=T)
for(i in file.list){
  
  out=gsub("simul/gwas/","simul/ldhub/",gsub(".gz","",i))
  
  wilson2ldhub(file=i,out = out)
  
  
}

### Heritability
file.list=system("ls simul/ldhub/*mung*.gz",intern=T)
for(i in file.list){
system(paste0("ldsc.py --h2 ",i," --ref-ld-chr /opt/working/wilson/app_full_stuff/ldsc/eur_w_ld_chr/ --w-ld-chr /opt/working/wilson/app_full_stuff/ldsc/eur_w_ld_chr/  --out ",gsub(".sumstats.gz","_h2",i)))
}
file.list=system("ls simul/ldhub/*h2*.log",intern=T)

h2.estimates=c()
for(i in file.list){
  
  res=fread(i,data.table=F,skip=25,nrow=1)
  res=cbind(i,res)
  h2.estimates=rbind(h2.estimates,res)
}


##rgs

file.list=system("ls simul/ldhub/*_mung.sumstats.gz",intern=T)
file.list=file.list[grep("height|mcv",file.list)]

file.height=file.list[grep("height",file.list)]
file.mcv=file.list[grep("mcv",file.list)]

library(magrittr)
for(i in file.height){
  for(j in file.mcv){
    tratto1=gsub("_mung.sumstats.gz","",i)%>%gsub("simul/ldhub/","",.)
    tratto2=gsub("_mung.sumstats.gz","",j)%>%gsub("simul/ldhub/","",.)
    system(paste0("ldsc.py --rg ",i,",",j," --ref-ld-chr /opt/working/wilson/app_full_stuff/ldsc/eur_w_ld_chr/ --w-ld-chr /opt/working/wilson/app_full_stuff/ldsc/eur_w_ld_chr/  --out simul/ldhub/rg_",tratto1,"_",tratto2,".out"))
  }
}


rgs.list=system("ls simul/ldhub/rg_*.out.log",intern=T)
rg.all=c()
for(i in rgs.list){
  
  rg.all=rbind(rg.all,fread(i,skip=60,nrows = 1,data.table=F))
  
  
}
library(dplyr)
rg.all$p1=recode(rg.all$p1,"simul/ldhub/height_4_dummy_mung.sumstats.gz"="Height OR 4 ~ dummy sex",
                 "simul/ldhub/height_4_plain_mung.sumstats.gz"="Height OR 4",
                 "simul/ldhub/height_8_dummy_mung.sumstats.gz"="Height OR 8 ~ dummy sex",
                 "simul/ldhub/height_8_plain_mung.sumstats.gz"="Height OR 8",
                 "simul/ldhub/height_original_mung.sumstats.gz"="Height all")


rg.all$p2=recode(rg.all$p2,"simul/ldhub/mcv_4_dummy_mung.sumstats.gz"="MCV OR 4 ~ dummy sex",
                 "simul/ldhub/mcv_4_plain_mung.sumstats.gz"="MCV OR 4",
                 "simul/ldhub/mcv_8_dummy_mung.sumstats.gz"="MCV OR 8 ~ dummy sex",
                 "simul/ldhub/mcv_8_plain_mung.sumstats.gz"="MCV OR 8",
                 "simul/ldhub/mcv_original_mung.sumstats.gz"="MCV all")
rg.all$lab=paste0(round(rg.all$rg,2),"\np=",formatC(rg.all$p, format = "e", digits = 1))

pdf("simul/ldhub/RG_hetamap.pdf",height=9,width=9)
ggplot(rg.all, aes(x=p1, y=p2,fill=rg)) + 
  geom_tile(colour = "black") + 
  scale_fill_gradient(low = "steelblue",  high = "white")+
  geom_text(aes(label = lab))+theme_minimal()+theme(axis.text.x = element_text(angle = 90))
dev.off()


pheno=fread("simul/phenotype_simul.txt")

png("simul/all_samples.png")
plot(pheno$st.height,pheno$st.mcv)
dev.off()
png("simul/OR4_samples.png")
plot(pheno$st.height_4 ,pheno$st.mcv_4)
dev.off()

png("simul/OR8_samples.png")
plot(pheno$st.height_8 ,pheno$st.mcv_8)
dev.off()

pdf("st.height_density.pdf")
ggarrange(ggplot(pheno, aes(x=st.height)) +
  geom_density()+theme_minimal()+xlim(-10,5),
ggplot(pheno, aes(x=st.height_4)) +
  geom_density()+theme_minimal()+xlim(-10,5),
ggplot(pheno, aes(x=st.height_8)) +
  geom_density()+theme_minimal()+xlim(-10,5),nrow=3)
dev.off()

pdf("st.mcv_density.pdf")
ggarrange(ggplot(pheno, aes(x=st.mcv)) +
  geom_density()+theme_minimal()+xlim(-5,15),
ggplot(pheno, aes(x=st.mcv_4)) +
  geom_density()+theme_minimal()+xlim(-5,15),
ggplot(pheno, aes(x=st.mcv_8)) +
  geom_density()+theme_minimal()+xlim(-5,15),nrow=3)
dev.off()

pdf("st.height_density_sex.pdf")
ggarrange(ggplot(pheno, aes(x=st.height,color=as.factor(dummy.sex))) +
            geom_density()+theme_minimal()+xlim(-10,5),
          ggplot(pheno, aes(x=st.height_4,color=as.factor(dummy.sex))) +
            geom_density()+theme_minimal()+xlim(-10,5),
          ggplot(pheno, aes(x=st.height_8,color=as.factor(dummy.sex))) +
            geom_density()+theme_minimal()+xlim(-10,5),nrow=3)
dev.off()
pdf("st.mcv_density_sex.pdf")
ggarrange(ggplot(pheno, aes(x=st.mcv,color=as.factor(dummy.sex))) +
            geom_density()+theme_minimal()+xlim(-5,15),
          ggplot(pheno, aes(x=st.mcv_4,color=as.factor(dummy.sex))) +
            geom_density()+theme_minimal()+xlim(-5,15),
          ggplot(pheno, aes(x=st.mcv_8,color=as.factor(dummy.sex))) +
            geom_density()+theme_minimal()+xlim(-5,15),nrow=3)
dev.off()


png("scatter_plot.png",height=1024,width=1024)
ggplot(pheno, aes(x=st.mcv,y=st.height,color=as.factor(dummy.sex)))+geom_point(alpha=0.5)+theme_minimal()+ylim(-6,3)
dev.off()
png("scatter_plot_4.png",height=1024,width=1024)
ggplot(pheno, aes(x=st.mcv_4,y=st.height_4,color=as.factor(dummy.sex)))+geom_point(alpha=0.5)+theme_minimal()+ylim(-6,3)
dev.off()
png("scatter_plot_8.png",height=1024,width=1024)
ggplot(pheno, aes(x=st.mcv_8,y=st.height_8,color=as.factor(dummy.sex)))+geom_point(alpha=0.5)+theme_minimal()+ylim(-6,3)
dev.off()






