library(data.table)
res=fread("simul/real_sex/phenotype_simul.txt")
library(ggplot2)

summary(lm(st.mcv_sc1_h_4_m_4~dummy.sex,data=res))


library(ggplot2)
library(ggpubr)
cols <- c("1" = "blue", "0" = "red", "NA" = NA)


for(i in c(1,2,3)){
  
  idx=names(res)[grep(paste0("sc",i),names(res))]
  res.tmp=res[,..idx]
  res.scen=data.table()
  for(j in c(1.2,1.5,1.8,2,4)){
     for(k in c(1.2,1.5,1.8,2,4)){
       
       h.name=paste0("st.height_sc",i,"_h_",j,"_m_",k)
       h.val=res.tmp[,..h.name]
       m.name=paste0("st.mcv_sc",i,"_h_",j,"_m_",k)
       m.val=res.tmp[,..m.name]
       d.name=paste0("dummy.sex_sc",i,"_h_",j,"_m_",k)
       d.val=res.tmp[,..d.name]
       tot=data.table(scenario=i,height=j,mcv=k,h.val,m.val,d.val)
       names(tot)=c("scenario","height","mcv","height.phen","mcv.phen","sex.phen")
       tot=na.omit(tot)
       res.scen=rbind(res.scen,tot)
       
    }
  }
  p1=ggplot(res.scen,aes(x=height.phen,y=mcv.phen))+stat_density_2d(aes(fill = as.factor(sex.phen)), geom = "polygon", bins=6, alpha=.3)+
    geom_hline(yintercept=0, linetype="dashed",color="black")+geom_vline(xintercept=0, linetype="dashed",color="black")
    theme_minimal()+scale_fill_manual(values=c("red","blue"),labels=c("Women","Men"))+facet_grid(as.factor(mcv)~as.factor(height))
    
  +theme(legend.title = element_text("Gender"),axis.title.x = element_text("Height"), axis.title.y =element_text("MCV") )
  
  pdf(paste0("pheno_sc",i,".pdf"),width=18,height=14)
  print(p1)
  dev.off()
  
  
}


for(i in c(1,2,3)){
  
  idx=names(res)[grep(paste0("sc",i),names(res))]
  res.tmp=res[,..idx]
  res.scen=data.table()
  for(j in c(1.2,1.5,1.8,2,4)){
    for(k in c(1.2,1.5,1.8,2,4)){
      
      h.name=paste0("st.height_sc",i,"_h_",j,"_m_",k)
      d.name=paste0("dummy.sex_sc",i,"_h_",j,"_m_",k)
      m.name=paste0("st.mcv_sc",i,"_h_",j,"_m_",k)
      h.val=resid(lm(as.formula(paste(h.name,"~",d.name)),data=res.tmp,na.action="na.exclude"))
      m.val=resid(lm(as.formula(paste(m.name,"~",d.name)),data=res.tmp,na.action="na.exclude"))
      d.val=res.tmp[,..d.name]
      tot=data.table(scenario=i,height=j,mcv=k,h.val,m.val,d.val)
      names(tot)=c("scenario","height","mcv","height.phen","mcv.phen","sex.phen")
      tot=na.omit(tot)
      res.scen=rbind(res.scen,tot)
      
    }
  }
  res.scen$mcv=as.factor(res.scen$mcv)
  res.scen$height=as.factor(res.scen$height)
  idx=sample(1:nrow(res.scen),size = 10000)
  p1=ggplot(res.scen[idx,],aes(x=height.phen,y=mcv.phen))+stat_density_2d(aes(fill = as.factor(sex.phen)), geom = "polygon", bins=6, alpha=.3)+
    theme_minimal()+theme(legend.title = element_text("Gender"),axis.title.x = element_text("Height"), axis.title.y =element_text("MCV") , )+
    scale_fill_manual(values=c("red","blue"),labels=c("Women","Men"))+facet_grid(mcv~height)
    
  pdf(paste0("pheno_sc",i,"_sex.pdf"),width=18,height=14)
  print(p1)
  dev.off()
  
  
}



p1=ggplot(iris,aes(x=Sepal.Length,y=Sepal.Width))+stat_density_2d(aes(fill = as.factor(Species)), geom = "polygon", bins=6, alpha=.3)+geom_hline(yintercept = 3)+
  theme_minimal()+theme(legend.title = element_text("Gender"),axis.title.x = element_text("Height"), axis.title.y =element_text("MCV")  )+facet_grid(iris$Species)




