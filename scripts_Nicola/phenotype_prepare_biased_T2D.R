source("/exports/igmm/eddie/wilson-lab/apps_by_us_full_stuff/gwas_pipeline/scripts/gwas_helper_functions.R") # Load the functions
library(data.table) 
traslation.file=fread("/exports/igmm/eddie/wilson-lab/data/base_data/ukbb_19655/phenotypes/ukbb_19655_translation.csv") # load translation file
ph=read_in_pheno( file_name="/exports/igmm/eddie/wilson-lab/data/base_data/ukbb_19655/phenotypes/ukbb_19655_base_phenotypes.tsv" # BAse file
                  #,nrows=10000
                  # number of rows
                  ,required_variables=c("body_mass_index_bmi","qual","sex","yob","ethnicity") # wanted variables
                  ,translation_df=traslation.file
) # translation file

ph.g=fread("zcat /exports/igmm/eddie/wilson-lab/data/base_data/ukbb_19655/phenotypes/ukbb_19655_genetic_phenotypes.tsv")
ph=merge(ph,ph.g,by="iid")
ph.t2d=fread("/exports/igmm/eddie/wilson-lab/projects/prj_109_sex_gwas/MR_simul/T2D_phenotype.tsv")
ph=merge(ph,ph.t2d,by="iid")


#xform_bespoke <- function(ph){
  

library("reshape2")

yn_maybe_to_tf <- function(answer_full) {
  
  ans<-tolower(substr(answer_full,1,1))
  
  res<-ifelse(ans=="y",T,NA)
  res<-ifelse(ans=="n",F,res)
  
  print(table(answer_full,res,exclude=NULL))
  return(res)
}



make_nice_var <- function(x) {
  x <- make.names(x)
  x <- tolower(x)
  x<- gsub("\\.","_",x)
  return(x)
}


make_array_of_values_binary <- function(df,array_name,prefix="",false_to=FALSE) {
  
  to_false_to <- function(x,false_to){
    x[which(x==FALSE)]<-false_to
    return(x)
  }
  
  array_cols <- grepl(array_name,names(df))
  wanted<-df[,(names(df)=="iid")|array_cols]
  
  #df<-df[,!array_cols]
  
  wanted_mel <- melt(wanted,id="iid")
  res<- dcast(wanted_mel, iid~value,fun.aggregate=length)
  res<-res[,names(res) != "NA"]
  val_cols<- 2:ncol(res)
  res[,val_cols] <- lapply(res[,val_cols],as.logical)
  res[,val_cols] <- lapply(res[,val_cols],to_false_to,false_to=false_to)
  if(prefix != "") names(res)[val_cols]<- paste(prefix,names(res)[val_cols])
  names(res)[val_cols] <- make_nice_var(names(res)[val_cols])
  
  heading("made the following and about to add to df")
  print(head(res))
  df<-merge_pj_many(df,res,id_var="iid")
  
  
  return(df)
}


string_to_bool <- function(string_vec,trues,falses,char_len=1) {
  
  string_vec <- tolower(string_vec)
  trues <- tolower(trues)
  falses <- tolower(falses)
  
  string_vec <- substr(string_vec,1,char_len)
  
  res<-ifelse(string_vec %in% trues, T,NA)
  
  res<-ifelse(string_vec %in% falses, F,res)
  
  return(res)
}

proc_yn_missing <- function(x) {
  
  y<- string_to_bool(x,"y","n")
  
  return(y)
}


proc_num_missing <- function(x) {
  
  y<-ifelse(x==-1 | x== -3,NA,x)
  return(y)
}



year_from_date <- function(dates) {
  if(class(dates)=="character") dates <- as.Date(dates)
  #print(dates)
  #print(class(dates))
  return( as.numeric(format(dates, "%Y")))
}

heading("make eth =TRUE/NA vars")

ph$ethnicity[which(ph$gen_eth==1)] <- "g_bri"
ph$ethnicity<-substr(ph$ethnicity,1,5)
ph<- make_array_of_values_binary(ph,"ethnicity",false_to=NA)

ph$yob <- ph$yob-1900



### Start here simulations
## Define traits


ph$bmi=ph$body_mass_index_bmi

set.seed(123456)

ph$sex=ifelse(ph$sex=="Male",1,0)

#ph$dummy.sex=rbinom(n=nrow(ph),size=1,prob=0.5)
#print(table(ph$dummy.sex))

library(dplyr)
recoded= recode(ph$qual,'College or University degree'=20 ,'A levels/AS levels or equivalent'=15 , 'O levels/GCSEs or equivalent'=13 , 'CSEs or equivalent'=12 ,
                                 'NVQ or HND or HNC or equivalent'=19 , 'Other professional qualifications eg: nursing, teaching'=17 , 'None of the above'=6)


ph$eduYears=recoded


ph$eduYears=(ph$eduYears-mean(ph$eduYears,na.rm=T))/sd(ph$eduYears,na.rm=T)

ph$bmi=(ph$bmi-mean(ph$bmi,na.rm=T))/sd(ph$bmi,na.rm=T)
ph$bmi[ph$bmi>4]=NA



############## No bias
analysis.plan=data.frame(anal_id="standard_covariates", 
                         formula="batch + array + g_bri + unrelated + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9+ pc10 + pc11 + pc12 + pc13 + pc14 + pc15 + pc16 + pc17 + pc18 + pc19 + pc20 + pc21 + pc22 + pc23 + pc24 + pc25 + pc26 + pc27 + pc28 + pc29 + pc30 + pc31 + pc32 + pc33 + pc34 + pc35 + pc36 + pc37 + pc38 + pc39 + pc40", 
                         analysis_type=NA,
                         residual_transform=NA,
                         fit_mixed_model=NA,stringsAsFactors = F)

### T2D
# All
ph$t2d_or=ph$fc1_poss_t2dm
an.pl=c("t2d_or","t2d_or~age+age2+standard_covariates","quantitative_trait","none","FALSE")
analysis.plan=rbind(analysis.plan,an.pl)
an.pl=c("t2d_or_sex","t2d_or~sex+age+age2+standard_covariates","quantitative_trait","none","FALSE")
analysis.plan=rbind(analysis.plan,an.pl)

# Female
ph$t2d_or_f=ph$fc1_poss_t2dm
ph$t2d_or_f[ph$sex==1]=NA
an.pl=c("t2d_or_f","t2d_or_f~age+age2+standard_covariates","quantitative_trait","none","FALSE")
analysis.plan=rbind(analysis.plan,an.pl)

#Malew
ph$t2d_or_m=ph$fc1_poss_t2dm
ph$t2d_or_m[ph$sex==0]=NA
an.pl=c("t2d_or_m","t2d_or_m~age+age2+standard_covariates","quantitative_trait","none","FALSE")
analysis.plan=rbind(analysis.plan,an.pl)

### BMI
# All
an.pl=c("bmi_or","bmi~age+age2+standard_covariates","quantitative_trait","none","FALSE")
analysis.plan=rbind(analysis.plan,an.pl)
an.pl=c("bmi_or_sex","bmi~sex+age+age2+standard_covariates","quantitative_trait","none","FALSE")
analysis.plan=rbind(analysis.plan,an.pl)

# Female
ph$bmi_or_f=ph$bmi
ph$bmi_or_f[ph$sex==1]=NA
an.pl=c("bmi_or_f","bmi_or_f~age+age2+standard_covariates","quantitative_trait","none","FALSE")
analysis.plan=rbind(analysis.plan,an.pl)

#Malew
ph$bmi_or_m=ph$bmi
ph$bmi_or_m[ph$sex==0]=NA
an.pl=c("bmi_or_m","bmi_or_m~age+age2+standard_covariates","quantitative_trait","none","FALSE")
analysis.plan=rbind(analysis.plan,an.pl)


### Edu
# All
an.pl=c("edu_or","eduYears~age+age2+standard_covariates","quantitative_trait","none","FALSE")
analysis.plan=rbind(analysis.plan,an.pl)
an.pl=c("edu_or_sex","eduYears~sex+age+age2+standard_covariates","quantitative_trait","none","FALSE")
analysis.plan=rbind(analysis.plan,an.pl)

# Female
ph$edu_or_f=ph$eduYears
ph$edu_or_f[ph$sex==1]=NA
an.pl=c("edu_or_f","edu_or_f~age+age2+standard_covariates","quantitative_trait","none","FALSE")
analysis.plan=rbind(analysis.plan,an.pl)

#Malew
ph$edu_or_m=ph$eduYears
ph$edu_or_m[ph$sex==0]=NA
an.pl=c("edu_or_m","edu_or_m~age+age2+standard_covariates","quantitative_trait","none","FALSE")
analysis.plan=rbind(analysis.plan,an.pl)


#####################
#
#  BIAS
######################

parameter.values=c(exp(-log(rev(c(1,1.5,1.8,2,4)))),exp(log(c(1.5,1.8,2,4))))

### Scenario 1
for( OR.edu in parameter.values){
  for( OR.bmi in parameter.values){
    
    ext=paste0("sc1_bmi_",OR.bmi,"_edu_",OR.edu) 
    
    set.seed(123456)
    z.pmen=(ph$eduYears*log(OR.edu))+(ph$bmi*log(OR.bmi))
    z.pmen[ph$dummy.sex==0]=z.pmen[ph$dummy.sex==0]*0
    prob.men=1/(1+exp(-z.pmen))
    ph$select=rbinom(n=nrow(ph),size = 1,prob = prob.men)
   
    ### New Edu
#    ph$new=ph$eduYears
#    ph$new[which(ph$select==0)]=NA
#    names(ph)[names(ph)=="new"]=paste(ext,"edu_bias",sep="_")
#    an.pl=c(paste(ext,"edu_bias",sep="_"),paste(paste(ext,"edu_bias",sep="_"),"~age+age2+standard_covariates"),"quantitative_trait","none","FALSE")
#    analysis.plan=rbind(analysis.plan,an.pl)
#    an.pl=c(paste(ext,"edu_bias_sex",sep="_"),paste(paste(ext,"edu_bias",sep="_"),"~sex+age+age2+standard_covariates"),"quantitative_trait","none","FALSE")
#    analysis.plan=rbind(analysis.plan,an.pl)
            
#    # New edu F
#    
#    ph$new=ph$eduYears
#    ph$new[which(ph$select==0 & ph$sex==1)]=NA
#    names(ph)[names(ph)=="new"]=paste(ext,"edu_bias_f",sep="_")
#    an.pl=c(paste(ext,"edu_bias_f",sep="_"),paste(paste(ext,"edu_bias_f",sep="_"),"~age+age2+standard_covariates"),"quantitative_trait","none","FALSE")
#    analysis.plan=rbind(analysis.plan,an.pl)
#    
    # New Edu M
    
#    ph$new=ph$eduYears
#    ph$new[which(ph$select==0 & ph$sex==0)]=NA
#    names(ph)[names(ph)=="new"]=paste(ext,"edu_bias_m",sep="_")
#    an.pl=c(paste(ext,"edu_bias_m",sep="_"),paste(paste(ext,"edu_bias_m",sep="_"),"~age+age2+standard_covariates"),"quantitative_trait","none","FALSE")
#    analysis.plan=rbind(analysis.plan,an.pl)
#    
    ### New BMI
    
#    ph$new=ph$bmi
#    ph$new[which(ph$select==0)]=NA
#    names(ph)[names(ph)=="new"]=paste(ext,"bmi_bias",sep="_")
#    an.pl=c(paste(ext,"bmi_bias",sep="_"),paste(paste(ext,"bmi_bias",sep="_"),"~age+age2+standard_covariates"),"quantitative_trait","none","FALSE")
#    analysis.plan=rbind(analysis.plan,an.pl)
#    an.pl=c(paste(ext,"bmi_bias_sex",sep="_"),paste(paste(ext,"bmi_bias",sep="_"),"~sex+age+age2+standard_covariates"),"quantitative_trait","none","FALSE")
#    analysis.plan=rbind(analysis.plan,an.pl)
    
    # New edu F
    
#    ph$new=ph$bmi
#    ph$new[which(ph$select==0 & ph$sex==1)]=NA
#    names(ph)[names(ph)=="new"]=paste(ext,"bmi_bias_f",sep="_")
#    an.pl=c(paste(ext,"bmi_bias_f",sep="_"),paste(paste(ext,"bmi_bias_f",sep="_"),"~age+age2+standard_covariates"),"quantitative_trait","none","FALSE")
#    analysis.plan=rbind(analysis.plan,an.pl)
    
    # New Edu M
    
#    ph$new=ph$bmi
#    ph$new[which(ph$select==0 & ph$sex==0)]=NA
#    names(ph)[names(ph)=="new"]=paste(ext,"bmi_bias_m",sep="_")
#    an.pl=c(paste(ext,"bmi_bias_m",sep="_"),paste(paste(ext,"bmi_bias_m",sep="_"),"~age+age2+standard_covariates"),"quantitative_trait","none","FALSE")
#    analysis.plan=rbind(analysis.plan,an.pl)
    
    ### New T2D
    
    ph$new=ph$fc1_poss_t2dm
    ph$new[which(ph$select==0)]=NA
    names(ph)[names(ph)=="new"]=paste(ext,"t2d_bias",sep="_")
    an.pl=c(paste(ext,"t2d_bias",sep="_"),paste(paste(ext,"t2d_bias",sep="_"),"~age+age2+standard_covariates"),"quantitative_trait","none","FALSE")
    analysis.plan=rbind(analysis.plan,an.pl)
    an.pl=c(paste(ext,"t2d_bias_sex",sep="_"),paste(paste(ext,"t2d_bias",sep="_"),"~sex+age+age2+standard_covariates"),"quantitative_trait","none","FALSE")
    analysis.plan=rbind(analysis.plan,an.pl)
    
    # New T2D F
    
    ph$new=ph$fc1_poss_t2dm
    ph$new[which(ph$select==0 & ph$sex==1)]=NA
    names(ph)[names(ph)=="new"]=paste(ext,"t2d_bias_f",sep="_")
    an.pl=c(paste(ext,"t2d_bias_f",sep="_"),paste(paste(ext,"t2d_bias_f",sep="_"),"~age+age2+standard_covariates"),"quantitative_trait","none","FALSE")
    analysis.plan=rbind(analysis.plan,an.pl)
    
    # New T2D M
    
    ph$new=ph$fc1_poss_t2dm
    ph$new[which(ph$select==0 & ph$sex==0)]=NA
    names(ph)[names(ph)=="new"]=paste(ext,"t2d_bias_m",sep="_")
    an.pl=c(paste(ext,"t2d_bias_m",sep="_"),paste(paste(ext,"t2d_bias_m",sep="_"),"~age+age2+standard_covariates"),"quantitative_trait","none","FALSE")
    analysis.plan=rbind(analysis.plan,an.pl)
  
    names(ph)[which(names(ph)=="select")]=paste(ext,"select",sep="_")
    #an.pl=c(paste(ext,"select",sep="_"),paste(paste(ext,"select",sep="_"),"~age+age2+standard_covariates"),"quantitative_trait","none","FALSE")
    #analysis.plan=rbind(analysis.plan,an.pl)
    
    
  }
}



### Scenario 2
for( OR.edu in parameter.values){
  for( OR.bmi in parameter.values){
    
    ext=paste0("sc2_bmi_",OR.bmi,"_edu_",OR.edu) 
    
    set.seed(123456)
    z.pmen=(ph$eduYears*log(OR.edu))+(ph$bmi*log(OR.bmi))
    z.pmen[ph$dummy.sex==0]=z.pmen[ph$dummy.sex==0]*-1
    prob.men=1/(1+exp(-z.pmen))
    ph$select=rbinom(n=nrow(ph),size = 1,prob = prob.men)
    
    ### New Edu
#    ph$new=ph$eduYears
#    ph$new[which(ph$select==0)]=NA
#    names(ph)[names(ph)=="new"]=paste(ext,"edu_bias",sep="_")
#    an.pl=c(paste(ext,"edu_bias",sep="_"),paste(paste(ext,"edu_bias",sep="_"),"~age+age2+standard_covariates"),"quantitative_trait","none","FALSE")
#    analysis.plan=rbind(analysis.plan,an.pl)
#    an.pl=c(paste(ext,"edu_bias_sex",sep="_"),paste(paste(ext,"edu_bias",sep="_"),"~sex+age+age2+standard_covariates"),"quantitative_trait","none","FALSE")
#    analysis.plan=rbind(analysis.plan,an.pl)
#    
#    # New edu F
#    
#    ph$new=ph$eduYears
#    ph$new[which(ph$select==0 & ph$sex==1)]=NA
#    names(ph)[names(ph)=="new"]=paste(ext,"edu_bias_f",sep="_")
#    an.pl=c(paste(ext,"edu_bias_f",sep="_"),paste(paste(ext,"edu_bias_f",sep="_"),"~age+age2+standard_covariates"),"quantitative_trait","none","FALSE")
#    analysis.plan=rbind(analysis.plan,an.pl)
#    
#    # New Edu M
#    
#    ph$new=ph$eduYears
#    ph$new[which(ph$select==0 & ph$sex==0)]=NA
#    names(ph)[names(ph)=="new"]=paste(ext,"edu_bias_m",sep="_")
#    an.pl=c(paste(ext,"edu_bias_m",sep="_"),paste(paste(ext,"edu_bias_m",sep="_"),"~age+age2+standard_covariates"),"quantitative_trait","none","FALSE")
#    analysis.plan=rbind(analysis.plan,an.pl)
#    
#    ### New BMI
#    
#    ph$new=ph$bmi
#    ph$new[which(ph$select==0)]=NA
#    names(ph)[names(ph)=="new"]=paste(ext,"bmi_bias",sep="_")
#    an.pl=c(paste(ext,"bmi_bias",sep="_"),paste(paste(ext,"bmi_bias",sep="_"),"~age+age2+standard_covariates"),"quantitative_trait","none","FALSE")
#    analysis.plan=rbind(analysis.plan,an.pl)
#    an.pl=c(paste(ext,"bmi_bias_sex",sep="_"),paste(paste(ext,"bmi_bias",sep="_"),"~sex+age+age2+standard_covariates"),"quantitative_trait","none","FALSE")
#    analysis.plan=rbind(analysis.plan,an.pl)
#    
#    # New edu F
#    
#    ph$new=ph$bmi
#    ph$new[which(ph$select==0 & ph$sex==1)]=NA
#    names(ph)[names(ph)=="new"]=paste(ext,"bmi_bias_f",sep="_")
#    an.pl=c(paste(ext,"bmi_bias_f",sep="_"),paste(paste(ext,"bmi_bias_f",sep="_"),"~age+age2+standard_covariates"),"quantitative_trait","none","FALSE")
#    analysis.plan=rbind(analysis.plan,an.pl)
#    
#    # New Edu M
#    
#    ph$new=ph$bmi
#    ph$new[which(ph$select==0 & ph$sex==0)]=NA
#    names(ph)[names(ph)=="new"]=paste(ext,"bmi_bias_m",sep="_")
#    an.pl=c(paste(ext,"bmi_bias_m",sep="_"),paste(paste(ext,"bmi_bias_m",sep="_"),"~age+age2+standard_covariates"),"quantitative_trait","none","FALSE")
#    analysis.plan=rbind(analysis.plan,an.pl)
#    
#    ### New T2D
    
    ph$new=ph$fc1_poss_t2dm
    ph$new[which(ph$select==0)]=NA
    names(ph)[names(ph)=="new"]=paste(ext,"t2d_bias",sep="_")
    an.pl=c(paste(ext,"t2d_bias",sep="_"),paste(paste(ext,"t2d_bias",sep="_"),"~age+age2+standard_covariates"),"quantitative_trait","none","FALSE")
    analysis.plan=rbind(analysis.plan,an.pl)
    an.pl=c(paste(ext,"t2d_bias_sex",sep="_"),paste(paste(ext,"t2d_bias",sep="_"),"~sex+age+age2+standard_covariates"),"quantitative_trait","none","FALSE")
    analysis.plan=rbind(analysis.plan,an.pl)
    
    # New T2D F
    
    ph$new=ph$fc1_poss_t2dm
    ph$new[which(ph$select==0 & ph$sex==1)]=NA
    names(ph)[names(ph)=="new"]=paste(ext,"t2d_bias_f",sep="_")
    an.pl=c(paste(ext,"t2d_bias_f",sep="_"),paste(paste(ext,"t2d_bias_f",sep="_"),"~age+age2+standard_covariates"),"quantitative_trait","none","FALSE")
    analysis.plan=rbind(analysis.plan,an.pl)
    
    # New T2D M
    
    ph$new=ph$fc1_poss_t2dm
    ph$new[which(ph$select==0 & ph$sex==0)]=NA
    names(ph)[names(ph)=="new"]=paste(ext,"t2d_bias_m",sep="_")
    an.pl=c(paste(ext,"t2d_bias_m",sep="_"),paste(paste(ext,"t2d_bias_m",sep="_"),"~age+age2+standard_covariates"),"quantitative_trait","none","FALSE")
    analysis.plan=rbind(analysis.plan,an.pl)
    
    names(ph)[which(names(ph)=="select")]=paste(ext,"select",sep="_")
    #an.pl=c(paste(ext,"select",sep="_"),paste(paste(ext,"select",sep="_"),"~age+age2+standard_covariates"),"quantitative_trait","none","FALSE")
    #analysis.plan=rbind(analysis.plan,an.pl)
    
    
  }
}





## Scenario  males -- females ++
#OR.edu=0.25 ## Bias OR for edu
#OR.bmi=0.25 ## Bias OR for bmi

#set.seed(123456)
#z.pmen=(ph$eduYears*log(OR.edu))*ph$dummy.sex+(ph$bmi*log(OR.bmi))*ph$dummy.sex
#prob.men=1/(1+exp(-z.pmen))
#ph$select=rbinom(n=nrow(ph),size = 1,prob = prob.men)

#ph$st.edu_new2=ph$eduYears
#ph$st.bmi_new2=ph$bmi
#ph$sex_new2=ph$dummy.sex

#ph$st.edu_new2[ph$select==0]=NA
#ph$st.bmi_new2[ph$select==0]=NA
#ph$sex_new2[ph$select==0]=NA
#### Create interaction term

#ph$st.edu_new2_male=ph$st.edu_new2
#ph$st.edu_new2_male[which(!is.na(ph$st.edu_new2_male) & ph$dummy.sex==0)]=0
#ph$st.edu_new2_female=ph$st.edu_new2
#ph$st.edu_new2_female[which(!is.na(ph$st.edu_new2_female) & ph$dummy.sex==1)]=0

#ph$st.bmi_new2_male=ph$st.bmi_new2
#ph$st.bmi_new2_male[which(!is.na(ph$st.bmi_new2_male) & ph$dummy.sex==0)]=0
#ph$st.bmi_new2_female=ph$st.bmi_new2
#ph$st.bmi_new2_female[which(!is.na(ph$st.bmi_new2_female) & ph$dummy.sex==1)]=0

write.table(analysis.plan,file="analysis_plan.txt",row.names=F,quote=F,sep="\t")

write.table(ph,file="phenotype_simul.txt",row.names=F,quote=F,sep="\t")


return(ph)

#}


