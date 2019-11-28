wilson2ldhub=function(file="simul/gwas/edu_sex.tsv.gz",out="simul/ldhub/edu_sex.sumstat"){
  require(data.table)
  data=fread(file,data.table=F)
  snp.list=fread("/exports/igmm/eddie/wilson-lab/apps_by_us_full_stuff/ldsc/data/w_hm3.snplist")
  data=data[which(data$rsid%in%snp.list$SNP),]
  data$Zscore=as.numeric(data$beta1)/as.numeric(data$se)
  data=data[,c("rsid","a1","a0","Zscore","n","p")]
  names(data)=c("snpid","A1","A2","Zscore","N","P-value")
  data.table:::fwrite(data,file=out,sep="\t",row.names=F,quote=F)
  system(paste0("munge_sumstats.py --sumstats ",out,"  --merge-alleles /exports/igmm/eddie/wilson-lab/apps_by_us_full_stuff/ldsc/data/w_hm3.snplist --out ",out,"_mung"))
  #system(paste("rm",out))
  
}

file.list=system("ls ../../GWAS/p05_comb_chr/*",intern=T)
#done.files=system("ls sumstats/*.sumstats.gz",intern=T)

missing=which(!(gsub("../../GWAS/p05_comb_chr/","",
               gsub(".gz","",file.list))%in%
          gsub("sumstats/","",
               gsub("_mung.sumstats.gz","",done.files))))
file.list=file.list[missing]

for(i in file.list){
  
  out=gsub("../../GWAS/p05_comb_chr/","sumstats/",gsub(".gz","",i))
  
  wilson2ldhub(file=i,out = out)
  
}

### Heritability
file.list=system("ls sumstats/*.gz",intern=T)
file.list=file.list[grep("bmi|edu",file.list)]

for(i in file.list){
  system(paste0("ldsc.py --h2 ",i," --ref-ld-chr /exports/igmm/eddie/wilson-lab/apps_by_us_full_stuff/ldsc/data/eur_w_ld_chr/ --w-ld-chr /exports/igmm/eddie/wilson-lab/apps_by_us_full_stuff/ldsc/data/eur_w_ld_chr/  --out ",gsub(".sumstats.gz","_h2",gsub("sumstats/","h2/",i))))
}
file.list=system("ls sumstats/*.gz",intern=T)
file.list=file.list[grep("bmi|edu",file.list)]
file.bmi=system("ls ../sumstats/bmi*_mung.sumstats.gz",intern=T)
library(magrittr)
for(i in file.list){
  for(j in file.bmi){
    tratto1=gsub(".tsv_mung.sumstats.gz","",i)%>%gsub("sumstats/","",.)
    tratto2=gsub(".tsv_mung.sumstats.gz","",j)%>%gsub("../sumstats/","",.)
      system(paste0("ldsc.py --rg ",i,",",j," --ref-ld-chr /exports/igmm/eddie/wilson-lab/apps_by_us_full_stuff/ldsc/data/eur_w_ld_chr/ --w-ld-chr /exports/igmm/eddie/wilson-lab/apps_by_us_full_stuff/ldsc/data/eur_w_ld_chr/  --out rg/rg_",tratto1,"_",tratto2,".out"))
  
  }
}




file.list=system("ls sumstats/*.gz",intern=T)
file.list=file.list[grep("height|mcv",file.list)]
file.height=file.list[grep("height",file.list)]
file.mcv=file.list[grep("mcv",file.list)]

library(magrittr)
for(i in file.height){
  for(j in file.mcv){
    tratto1=gsub(".tsv_mung.sumstats.gz","",i)%>%gsub("sumstats/","",.)
    tratto2=gsub(".tsv_mung.sumstats.gz","",j)%>%gsub("sumstats/","",.)
    if(gsub("height_","",i)==gsub("mcv_","",j)){
      system(paste0("ldsc.py --rg ",i,",",j," --ref-ld-chr /exports/igmm/eddie/wilson-lab/apps_by_us_full_stuff/ldsc/data/eur_w_ld_chr/ --w-ld-chr /exports/igmm/eddie/wilson-lab/apps_by_us_full_stuff/ldsc/data/eur_w_ld_chr/  --out rg/rg_",tratto1,"_",tratto2,".out"))
    }
  }
}


#### 

### rg matrices
library(data.table)
rg.files=system("ls simul/new_simul/rg/*",intern=T)

h2.tot=data.frame(sc=rep(NA,length(rg.files)),height=rep(NA,length(rg.files))
                  ,mcv=rep(NA,length(rg.files)),sex=rep(NA,length(rg.files)),rg=rep(NA,length(rg.files))
                  ,se=rep(NA,length(rg.files)),p=rep(NA,length(rg.files)))
library(magrittr)
k=1
for(i in rg.files){
  
  
  tmp=fread(i,data.table=F,nrows=1,skip=60)
  gsub("sumstats/","",tmp$p1)%>% strsplit(.,split="sc") %>% unlist -> a
  par=as.numeric(strsplit(a[2],split="_h_|_m_|.tsv")[[1]][-4])
  if(length(grep("sex",i))==1){
    param=c(par,1,tmp$rg,tmp$se,tmp$p)
  }else{
    param=c(par,0,tmp$rg,tmp$se,tmp$p)
    
  }
  
  h2.tot[k,]=param
  k=k+1
  
}

h2.tot$lab=paste0(round(h2.tot$rg,2),"\np=",formatC(h2.tot$p, format = "e", digits = 1))

library(ggplot2)
pdf("simul/new_simul/Simulation_results_rgs.pdf",height=14,width=14)
ggplot(h2.tot, aes(x=as.factor(mcv), y=as.factor(height),fill=rg)) + 
  geom_tile(colour = "black") + 
  scale_fill_gradient(low = "steelblue",  high = "white")+
  geom_text(aes(label = lab))+
  theme_minimal()+theme(axis.text.x = element_text(angle = 90))+facet_grid(sc~sex)
dev.off()




file.list=system("ls simul/new_simul/h2/*mung*.log",intern=T)
h2.estimates=c()
for(i in file.list){
  
  res=fread(i,data.table=F,skip=25,nrow=1)
  res=cbind(i,res)
  h2.estimates=rbind(h2.estimates,res)
}

h2.estimates$sc=NA
h2.estimates$h_or=NA
h2.estimates$m_or=NA
h2.estimates$sex=0
h2.estimates$trait=NA


h2.estimates$sc[grep("sc1",h2.estimates$i)]=1
h2.estimates$sc[grep("sc2",h2.estimates$i)]=2
h2.estimates$sc[grep("sc3",h2.estimates$i)]=3

h2.estimates$h_or[grep("h_1.2",h2.estimates$i)]=1.2
h2.estimates$h_or[grep("h_1.5",h2.estimates$i)]=1.5
h2.estimates$h_or[grep("h_1.8",h2.estimates$i)]=1.8
h2.estimates$h_or[grep("h_2",h2.estimates$i)]=2
h2.estimates$h_or[grep("h_4",h2.estimates$i)]=4

h2.estimates$m_or[grep("m_1.2",h2.estimates$i)]=1.2
h2.estimates$m_or[grep("m_1.5",h2.estimates$i)]=1.5
h2.estimates$m_or[grep("m_1.8",h2.estimates$i)]=1.8
h2.estimates$m_or[grep("m_2",h2.estimates$i)]=2
h2.estimates$m_or[grep("m_4",h2.estimates$i)]=4

h2.estimates$sex[grep("_sex",h2.estimates$i)]=1

h2.estimates$trait[grep("dummy.sex",h2.estimates$i)]="dummy.sex"
h2.estimates$trait[grep("height",h2.estimates$i)]="height"
h2.estimates$trait[grep("mcv",h2.estimates$i)]="mcv"

h2.estimates$V6=gsub("(","",h2.estimates$V6,fixed=T)%>%gsub(")","",.,fixed=T)%>%as.numeric
h2.estimates$high=h2.estimates$V5+(qnorm(p = 0.975)*h2.estimates$V6)
h2.estimates$low=h2.estimates$V5+(qnorm(p = 1-0.975)*h2.estimates$V6)



ggplot(h2.estimates[h2.estimates$trait=="height",], aes(x=as.factor(m_or), y=as.factor(h_or),fill=V5)) + 
  geom_tile(colour = "black") + 
  scale_fill_gradient(low = "white",  high = "steelblue")+
  #geom_text(aes(label = lab))+
  theme_minimal()+theme(axis.text.x = element_text(angle = 90))+facet_grid(sc~sex)
h2.estimates=h2.estimates[,-c(1:5)]
names(h2.estimates)[1:2]=c("h2","se")

library(ggpubr)
pdf("Heritability_plots.pdf",height=14,width=14)

  tmp=h2.estimates[h2.estimates$trait=="dummy.sex",]
  p1=ggplot(tmp[which(tmp$sc==1),],aes(y=h2,ymin=low,ymax=high,x=1))+geom_pointrange()+
    facet_grid(h_or~m_or)+geom_hline(yintercept=0,color="darkgrey")+xlim(0,2)+ylim(-0.01,0.09)+theme_minimal()+theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
  p2=ggplot(tmp[which(tmp$sc==2),],aes(y=h2,ymin=low,ymax=high,x=1))+geom_pointrange()+
    facet_grid(h_or~m_or)+geom_hline(yintercept=0,color="darkgrey")+xlim(0,2)+ylim(-0.01,0.09)+theme_minimal()+theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
  p3=ggplot(tmp[which(tmp$sc==3),],aes(y=h2,ymin=low,ymax=high,x=1))+geom_pointrange()+
    facet_grid(h_or~m_or)+geom_hline(yintercept=0,color="darkgrey")+xlim(0,2)+ylim(-0.01,0.09)+theme_minimal()+theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
  p4=ggplot(tmp,aes(x=1,y=1))+geom_text(label="Dummy sex",size=20)+theme_void()
  ggarrange(p1,p2,p3,p4,labels = c("Scenario 1","Scenario 2","Scenario 3",""),common.legend = TRUE)

  
  tmp=h2.estimates[h2.estimates$trait=="height",]
  p1=ggplot(tmp[which(tmp$sc==1 ),],aes(y=h2,ymin=low,ymax=high,x=1,shape=as.factor(sex),color=as.factor(sex)))+geom_pointrange(position=position_dodge(width = 0.5))+
  facet_grid(h_or~m_or)+geom_hline(yintercept=0,color="darkgrey")+xlim(0,2)+ylim(-0.1,0.3)+theme_minimal()+theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

  p2=ggplot(tmp[which(tmp$sc==2 ),],aes(y=h2,ymin=low,ymax=high,x=1,shape=as.factor(sex),color=as.factor(sex)))+geom_pointrange(position=position_dodge(width = 0.5))+
    facet_grid(h_or~m_or)+geom_hline(yintercept=0,color="darkgrey")+xlim(0,2)+ylim(-0.1,0.3)+theme_minimal()+theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

  p3=ggplot(tmp[which(tmp$sc==3 ),],aes(y=h2,ymin=low,ymax=high,x=1,shape=as.factor(sex),color=as.factor(sex)))+geom_pointrange(position=position_dodge(width = 0.5))+
    facet_grid(h_or~m_or)+geom_hline(yintercept=0,color="darkgrey")+xlim(0,2)+ylim(-0.1,0.3)+theme_minimal()+theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
  p4=ggplot(tmp,aes(x=1,y=1))+geom_text(label="Height",size=20)+theme_void()
  
  ggarrange(p1,p2,p3,p4,labels = c("Scenario 1","Scenario 2","Scenario 3",""),common.legend = TRUE)
  
  tmp=h2.estimates[h2.estimates$trait=="mcv",]
  p1=ggplot(tmp[which(tmp$sc==1 ),],aes(y=h2,ymin=low,ymax=high,x=1,shape=as.factor(sex),color=as.factor(sex)))+geom_pointrange(position=position_dodge(width = 0.5))+
    facet_grid(h_or~m_or)+geom_hline(yintercept=0,color="darkgrey")+xlim(0,2)+ylim(-0.1,0.3)+theme_minimal()+theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
  
  p2=ggplot(tmp[which(tmp$sc==2 ),],aes(y=h2,ymin=low,ymax=high,x=1,shape=as.factor(sex),color=as.factor(sex)))+geom_pointrange(position=position_dodge(width = 0.5))+
    facet_grid(h_or~m_or)+geom_hline(yintercept=0,color="darkgrey")+xlim(0,2)+ylim(-0.1,0.3)+theme_minimal()+theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
  
  p3=ggplot(tmp[which(tmp$sc==3 ),],aes(y=h2,ymin=low,ymax=high,x=1,shape=as.factor(sex),color=as.factor(sex)))+geom_pointrange(position=position_dodge(width = 0.5))+
    facet_grid(h_or~m_or)+geom_hline(yintercept=0,color="darkgrey")+xlim(0,2)+ylim(-0.1,0.3)+theme_minimal()+theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
  p4=ggplot(tmp,aes(x=1,y=1))+geom_text(label="MCV",size=20)+theme_void()
  
  ggarrange(p1,p2,p3,p4,labels = c("Scenario 1","Scenario 2","Scenario 3",""),common.legend = TRUE)
  
dev.off()

plot(tmp$V5,x=rep(1:25,3),col=tmp$sc,pch=19,type="n")
text(labels = paste(tmp$h_or,tmp$m_or),tmp$V5,x=rep(1:25,3),col=tmp$sc)



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




system(paste0("ldsc.py --h2  --ref-ld-chr /opt/working/wilson/app_full_stuff/ldsc/eur_w_ld_chr/ --w-ld-chr /opt/working/wilson/app_full_stuff/ldsc/eur_w_ld_chr/  --out ",gsub(".sumstats.gz","_h2",i)))


