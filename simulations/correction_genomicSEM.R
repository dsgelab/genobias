# # GenomicSEM: https://github.com/MichelNivard/GenomicSEM
require("GenomicSEM")

rgs <- c('-0.3','-0.1','0','0.1','0.3')
ORs <- c('1.2','1.5','1.8','2','3')

for (r in rgs) {
  for (o in ORs) {
    
    # Load munged sumstats for X*, Y*, S, pred(S)*, for each subsample (rg_OR) 
    setwd('/../../simulations_heckman/sumstats/')
    Ystar <- paste0("Y_star_",r,"_",o,".sumstats.gz")
    Xstar <- paste0("X_star_",r,"_",o,".sumstats.gz")
    predstar <- paste0("pred_star_",r,"_",o,".sumstats.gz")
    sel <- paste0("sel_",r,"_",o,".sumstats.gz")
    
    # Run multivariable LD-Score regression to obtain the genetic covariance (S) matrix 
    # and corresponding sampling covariance matrix (V)
    ldsc <- ldsc(traits = c(Ystar,Xstar,predstar,sel),
                 sample.prev = c(NA,NA,NA,.5),population.prev = c(NA,NA,NA,.5),
                 ld = "/ldsc/eur_w_ld_chr/",
                 wld= "/ldsc/eur_w_ld_chr/",
                 trait.names = c("Y","X","pred", "S"))
    
    model_bias <-  "Y ~~ X"
    biased <- usermodel(covstruc = ldsc,estimation = "DWLS",model = model_bias)
    
    # Biased rg
    rgB <- biased$results[2,]
    
    
    # # # # # # # # # # # # # # # # # 
    # # # Path model correction # # #
    # # # # # # # # # # # # # # # # # 
    
    model_condition_on_S_1 <- "
    X ~ Y
    X ~~ 0*S
    Y ~~ S
    "
    model_condition_on_S_2 <- "
    Y ~ X
    X ~~ S
    Y ~~ 0*S
    "
    a <- usermodel(covstruc = ldsc,estimation = "DWLS",model = model_condition_on_S_1)
    b <- usermodel(covstruc = ldsc,estimation = "DWLS",model = model_condition_on_S_2)
    chisq <- c(a$modelfit[1,1],b$modelfit[1,1])
    
    # ! Parameter order might change between different installation of genomicSEM:
    # Inspect objects a and b and be sure to take the line X ~ Y (or the reverse)
    rg <- rbind(a$results[4,],b$results[2,])
    
    # Corrected rg is the value from the model with best fit out of a and b:
    rgM <- rg[which.min(chisq),]
    
    
    # # # # # # # # # # # # # # # # # # # #
    # # # Heckman derived correction  # # #
    # # # # # # # # # # # # # # # # # # # #
    
    mod_correct <- " 
    Y~~a*X
    Y~~b*pred
    X~~c*pred
    fin:=(a-(b*c))/(1-c^2)
    "
    
    c <- usermodel(covstruc = ldsc,estimation = "DWLS",model = mod_correct)
    
    # Corrected rg
    rgH <- c$results[1,]
    
    RES <- NULL
    
    res <- unlist(c(r,o,rgB,rgM,rgH))
    
    RES <- rbind(RES, res)
    RES <- data.frame(RES)
    names(RES)[1:2] <- c("rg", "or")
    
    write.table(RES, '/../../simulations_heckman/results_correction.tsv', sep = '\t', row.names = F, append = T)
    
  }
}