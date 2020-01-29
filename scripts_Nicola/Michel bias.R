#library(devtools)
#install_github("MichelNivard/GenomicSEM")
require(GenomicSEM)
#require(semPlot)

sumstats=system("ls ../gencor/sumstats/*.gz",intern=T)
[1] "../gencor/sumstats/bmi_or.tsv_mung.sumstats.gz"            
[2] "../gencor/sumstats/bmi_sex.tsv_mung.sumstats.gz"           
[3] "../gencor/sumstats/dummy.sex.tsv_mung.sumstats.gz"         
[4] "../gencor/sumstats/edu_or.tsv_mung.sumstats.gz"            
[5] "../gencor/sumstats/edu_sex.tsv_mung.sumstats.gz"           
[6] "../gencor/sumstats/sex_new1.tsv_mung.sumstats.gz"          
[7] "../gencor/sumstats/st.bmi_female.tsv_mung.sumstats.gz"     
[8] "../gencor/sumstats/st.bmi_male.tsv_mung.sumstats.gz"       
[9] "../gencor/sumstats/st.bmi_new1_female.tsv_mung.sumstats.gz"
[10] "../gencor/sumstats/st.bmi_new1_male.tsv_mung.sumstats.gz"  
[11] "../gencor/sumstats/st.bmi_new1_sex.tsv_mung.sumstats.gz"   
[12] "../gencor/sumstats/st.bmi_new1.tsv_mung.sumstats.gz"       
[13] "../gencor/sumstats/st.edu_female.tsv_mung.sumstats.gz"     
[14] "../gencor/sumstats/st.edu_male.tsv_mung.sumstats.gz"       
[15] "../gencor/sumstats/st.edu_new1_female.tsv_mung.sumstats.gz"
[16] "../gencor/sumstats/st.edu_new1_male.tsv_mung.sumstats.gz"  
[17] "../gencor/sumstats/st.edu_new1_sex.tsv_mung.sumstats.gz"   
[18] "../gencor/sumstats/st.edu_new1.tsv_mung.sumstats.gz"       
> 
  
labels=c( "BMI_origin"            ,
          "BMI_origin_sex"           ,
          "dummy_sex"         ,
          "Edu_origin"            ,
          "Edu_origin_sex"           ,
          "dummy.sex_bias"          ,
          "BMI_origin_Female"     ,
          "BMI_origin_Male"       ,
          "BMI_bias_Female",
          "BMI_bias_Male"  ,
          "BMI_bias_sex"   ,
          "BMI_bias"       ,
          "Edu_origin_Female"     ,
          "Edu_origin_Male"       ,
          "Edu_bias_Female",
          "Edu_bias_Male"  ,
          "Edu_bias_sex"   ,
          "Edu_bias"       )

ld_out <- ldsc( traits = sumstats[c(2,3,5,6,9,10,15,16)] 
               ,trait.names = labels[c(2,3,5,6,9,10,15,16)]
               ,sample.prev=c(NA,NA,0.4999443,NA,NA,0.4965508,rep(NA,times=12))
               ,population.prev=c(NA,NA,0.4999443,NA,NA,0.4999443,rep(NA,times=12))
               , ld = "/exports/igmm/eddie/wilson-lab/apps_by_us_full_stuff/ldsc/data/eur_w_ld_chr/"
               ,wld = "/exports/igmm/eddie/wilson-lab/apps_by_us_full_stuff/ldsc/data/eur_w_ld_chr/")

save(ld.out,file="ld_out.Rdata")



model_true <- ' dummy_sex ~ BMI_origin_sex + Edu_origin_sex'
model_bias <- 'dummy.sex_bias ~ BMI_origin_sex + Edu_origin_sex '
model_Heckman <- 'dummy.sex_bias ~   BMI_origin_sex + Edu_origin_sex +BMI_bias_Female + BMI_bias_Male + Edu_bias_Female + Edu_bias_Male'
model_bias <- 'BMI_bias_Female ~  Edu_bias_Female '

truth <- usermodel(ld_out,model=model_true,estimation = "DWLS")
biased <- usermodel(ld_out,model=model_bias,estimation = "DWLS")

Heckman <- usermodel(ld_out,model=model_Heckman,estimation = "DWLS")

covstruc=ld_out
model=model_Heckman
estimation = "DWLS"
CFIcalc = TRUE
std.lv = FALSE 
imp_cov = FALSE


time <- proc.time()
test <- c(str_detect(model, "~"), str_detect(model, "="), 
          str_detect(model, "\\+"))
if (any(test) != TRUE) {
  warning("Your model name may be listed in quotes; please remove the quotes and try re-running if the function has returned an error about not locating the ReorderModel.")
}
rearrange <- function(k, fit, names) {
  order1 <- names
  order2 <- rownames(inspect(fit)[[1]])
  kst <- k * (k + 1)/2
  covA <- matrix(NA, k, k)
  covA[lower.tri(covA, diag = TRUE)] <- 1:kst
  covA <- t(covA)
  covA[lower.tri(covA, diag = TRUE)] <- 1:kst
  colnames(covA) <- rownames(covA) <- order1
  covA <- covA[order2, order2]
  vec2 <- lav_matrix_vech(covA)
  return(vec2)
}
V_LD <- as.matrix(covstruc[[1]])
S_LD <- as.matrix(covstruc[[2]])
k <- ncol(S_LD)
z <- (k * (k + 1))/2
write.names <- function(k, label = "V") {
  varnames <- vector(mode = "character", length = k)
  for (i in 1:k) {
    varnames[i] <- paste(label, i, sep = "")
  }
  return(varnames)
}







