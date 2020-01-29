library(RadialMR)
library(mr.raps)
create.instrument=function(data=res,ratio.col="corr2raw_ratio",p.col,a1.col="alt",
                           a0.col="ref",p.threshold=5e-8,
                           ratio.threshold=0.05,
                           phen.name="beef",beta.col="corr.beta",se.col="se"){
  
  idx=which(data[,p.col]<=p.threshold & 
              data[,ratio.col]>(1-ratio.threshold) &
              data[,ratio.col]<(1+ratio.threshold))
  
  if(length(idx)==0){
    instr=c()
  }else{
    data=data[idx,]
    
    prune=TwoSampleMR::clump_data(data.frame(SNP=data$rs,
                                             chr_name=data$chrm,
                                             chrom_start=data$pos,
                                             pval.exposure =data[,p.col]))
    data=data[data$rs%in%prune$SNP,]
    data$Phenotype=phen.name
    
    instr=TwoSampleMR::format_data(data,
                                   type ="exposure",
                                   snp_col="rs",
                                   beta_col=beta.col,
                                   se_col=se.col,eaf_col = "freq1"
                                   ,effect_allele_col = a1.col,other_allele_col = a0.col,
                                   pval_col = p.col,samplesize_col="n")
  }
  instr
}


estimate.standard.beta=function(z,n){
  
  r2 = (z^2)/(n-2-(z^2))
  beta=sqrt(r2)*sign(z)
  se=sqrt(1-r2)/sqrt(n)
  res=data.frame(beta=beta,se=se)
  res
}


het.pruner=function(harmon=harmon,exposure.id="fnNokw",outcome.id="9F0GHV",method="mr_ivw",target.p=0.05){
  current.p=NA
  tmp=harmon
  start=mr_heterogeneity(tmp,method_list = method,)$Q_pval
  bad=c()
  
  if(!is.null(start)){
    current.p=start
    while(current.p<target.p & nrow(tmp)>2){
      pval.string=c()
      for(i in 1:nrow(tmp)){
        pval.string=c(pval.string,mr_heterogeneity(tmp[-i,],method_list = method)$Q_pval)
      }
      bad=c(bad,as.character(tmp$SNP[which.max(pval.string)]))
      tmp=tmp[-which.max(pval.string),]
      current.p=mr_heterogeneity(tmp,method_list = method)$Q_pval
    }
  }
  output=list()
  output[["excluded"]]=bad
  output[["nsnps"]]=nrow(tmp)
  output[["perc.retained"]]=nrow(tmp)/(nrow(tmp)+length(bad))
  output[["Het_p"]]=current.p
  output[["harmonised_data"]]=tmp
  output
}



mr.run=function(harm.data,label="raw",mr.raps=TRUE){
  
  harm.data=harm.data[which(!is.na(harm.data$beta.outcome)),]
  res.raw=data.frame(exposure=first(harm.data$exposure),outcome=first(harm.data$outcome),method=NA,nsnp=NA,b=NA,se=NA,pval=NA,type=paste("No_instr",label,sep="_"),F=NA,outliers=NA,stringsAsFactors=F)
  
  n.init=nrow(harm.data)
  #harm.data=harm.data[harm.data$pval.outcome>1e-5,]
  n.excl=n.init-nrow(harm.data)
  het.p=NA
  ### Pruning 
  bad=NULL
  if(nrow(harm.data)>2){
    prun=ivw_radial2(format_radial(harm.data$beta.exposure,harm.data$beta.outcome,
                                   seBXG = harm.data$se.exposure,seBYG = harm.data$se.outcome,RSID = harm.data$SNP)
                     ,alpha = 0.05/nrow(harm.data))
    print(1)
    if(prun$outliers!="No significant outliers"){
      bad=prun$outliers$SNP[which.min(prun$outliers$p.value)]
      print(2)
      
      ### Remove the worse case and rerun
      tmp=harm.data[-which(harm.data$SNP%in%bad),]
      print(3)
      
      if(nrow(tmp)>2){
        prun=ivw_radial2(format_radial(tmp$beta.exposure,tmp$beta.outcome,seBXG = tmp$se.exposure,seBYG = tmp$se.outcome,RSID = tmp$SNP),alpha = 0.05/nrow(harm.data))
        print(4)
        
        if(prun$outliers!="No significant outliers"){
          prun=ivw_radial2(format_radial(harm.data$beta.exposure,harm.data$beta.outcome,
                                         seBXG = harm.data$se.exposure,seBYG = harm.data$se.outcome,RSID = harm.data$SNP)
                           ,alpha = 0.05/nrow(harm.data))
          print(5)
          
          bad=as.character(prun$outliers$SNP)
          harm.data=harm.data[which(!(harm.data$SNP%in%bad)),]
          print(6)
          
        }else{
          harm.data=tmp
          print(7)
          
        }
      }else{
        harm.data=tmp
        print(8)
        
      }
    }
    if(nrow(harm.data)>1){
      anal=ivw_radial2(format_radial(harm.data$beta.exposure,harm.data$beta.outcome,
                                   seBXG = harm.data$se.exposure,seBYG = harm.data$se.outcome,RSID = harm.data$SNP)
                     ,alpha = 0.05/nrow(harm.data))
      het.p=anal$Total_Q_chi
      if(anal$Total_Q_chi<(0.05/78)){
        anal=anal$re.coef
      }else{
        anal=anal$fe.coef
      
      }
    }
  }
   
  print(bad)
  n.out=length(bad)
  
  if(n.excl!=n.init){
    res.raw=data.frame(exposure=first(harm.data$exposure),outcome=first(harm.data$outcome),method=NA,nsnp=NA,b=NA,se=NA,pval=NA,type=paste("No_instr",label,sep="_"),F=NA,outliers=n.out,stringsAsFactors=F)
  }
  if(nrow(harm.data)>0){
    #F.R2=sum(((harm.data$beta.exposure/harm.data$se.exposure)^2)/harm.data$samplesize.exposure)
    #F.raw=(F.R2*(harm.data$samplesize.exposure[1]-1-nrow(harm.data)))/((1-F.R2)*nrow(harm.data))
    F.R2=NA
    F.raw=NA
    
    
    raw.pleio=mr_pleiotropy_test(harm.data)
    if(nrow(raw.pleio)==0) raw.pleio=data.frame(pval=NA)
  
    if(nrow(harm.data)>2){
      sens.anal=mr(harm.data,method_list = c("mr_egger_regression","mr_weighted_median"))
      meth="IVW (FE)"
      if(het.p<(0.05/129)){
        meth="IVW (RE)"
      }
      res.raw=data.frame(outcome=unique(harm.data$outcome),exposure=unique(harm.data$exposure[1]),method=meth
                         ,nsnp=nrow(harm.data),b=anal$Estimate,se=anal$Std.Error,pval=anal$`Pr(>|t|)`)
      if(nrow(sens.anal)>0){
        sens.anal=data.frame(beta_egger=sens.anal$b[which(sens.anal$method=="MR Egger")]
                             ,se_egger=sens.anal$se[which(sens.anal$method=="MR Egger")]
                             ,beta_median=sens.anal$b[which(sens.anal$method=="Weighted median")]
                             ,se_median=sens.anal$se[which(sens.anal$method=="Weighted median")])
      }else{
        sens.anal=data.frame(beta_egger=NA
                             ,se_egger=NA
                             ,beta_median=NA
                             ,se_median=NA)
      }
      if(mr.raps==T){
        
        
        mr.raps.res=mr.raps(harm.data)
        mr.raps.res=c(mr.raps.res$beta.hat,mr.raps.res$beta.se)
        
      }else{
        
        mr.raps.res=c(NA,NA)
        
      }
      
    }else if(nrow(harm.data)==2){
      res.raw=mr(harm.data,method_list = "mr_ivw")
      sens.anal=data.frame(beta_egger=NA
                           ,se_egger=NA
                           ,beta_median=NA
                           ,se_median=NA)
      mr.raps.res=c(NA,NA)
      
    }else{
      res.raw=mr(harm.data,method_list = "mr_wald_ratio")
      sens.anal=data.frame(beta_egger=NA
                           ,se_egger=NA
                           ,beta_median=NA
                           ,se_median=NA)
      mr.raps.res=c(NA,NA)
      
    }
    names(mr.raps.res)=c("mr.raps.beta","mr.raps.se")
    res.raw$type=label
    res.raw$F=F.raw
    res.raw$outliers=n.out
    res.raw$egger.pleio=raw.pleio$pval
    res.raw=data.frame(res.raw[,c("exposure","outcome","method","nsnp","b","se","pval","type","F","outliers","egger.pleio")],t(data.frame(mr.raps.res)),sens.anal,het.p)
  }
  res.raw$excl.p=n.excl
  res.raw
  
}


ivw_radial2<-function (r_input, alpha, weights, tol) 
{
  Ratios <- r_input[, 3]/r_input[, 2]
  F <- r_input[, 2]^2/r_input[, 4]^2
  mf <- mean(F)
  cat()
  if (missing(alpha)) {
    alpha <- 0.05
    warning("Significance threshold for outlier detection not specified: Adopting a 95% threshold")
  }
  if (missing(weights)) {
    weights <- 3
    warning("Weights not specified: Adopting modified second-order weights")
  }
  if (missing(tol)) {
    tol <- 1e-05
  }
  summary <- TRUE
  if (weights == 1) {
    W <- ((r_input[, 2]^2)/(r_input[, 5]^2))
  }
  if (weights == 2) {
    W <- ((r_input[, 5]^2/r_input[, 2]^2) + ((r_input[, 3]^2 * 
                                                r_input[, 4]^2)/r_input[, 2]^4))^-1
  }
  if (weights == 3) {
    W <- ((r_input[, 2]^2)/(r_input[, 5]^2))
  }
  Wj <- sqrt(W)
  BetaWj <- Ratios * Wj
  IVW.Model <- lm(BetaWj ~ -1 + Wj)
  EstimatesIVW <- summary(lm(IVW.Model))
  IVW.Slope <- EstimatesIVW$coefficients[1]
  IVW.SE <- EstimatesIVW$coefficients[2]
  IVW_CI <- confint(IVW.Model)
  DF <- length(r_input[, 1]) - 1
  Qj <- W * (Ratios - IVW.Slope)^2
  Total_Q <- sum(Qj)
  Total_Q_chi <- pchisq(Total_Q, length(r_input[, 2]) - 1, 
                        lower.tail = FALSE)
  if (weights == 3) {
    W <- ((r_input[, 5]^2 + (IVW.Slope^2 * r_input[, 4]^2))/r_input[, 
                                                                    2]^2)^-1
    Wj <- sqrt(W)
    BetaWj <- Ratios * Wj
    IVW.Model <- lm(BetaWj ~ -1 + Wj)
    EstimatesIVW <- summary(lm(BetaWj ~ -1 + Wj))
    IVW.Slope <- EstimatesIVW$coefficients[1]
    IVW.SE <- EstimatesIVW$coefficients[2]
    IVW_CI <- confint(IVW.Model)
    Qj <- W * (Ratios - IVW.Slope)^2
    Total_Q <- sum(Qj)
    Total_Q_chi <- pchisq(Total_Q, length(r_input[, 2]) - 
                            1, lower.tail = FALSE)
  }
  Iterative_ivw <- function(int.tol) {
    Diff <- 1
    Bhat1.Iterative <- 0
    count <- 0
    while (Diff >= tol) {
      W <- 1/(r_input[, 5]^2/r_input[, 2]^2 + (Bhat1.Iterative^2) * 
                r_input[, 4]^2/r_input[, 2]^2)
      Wj <- sqrt(W)
      BetaWj <- Ratios * Wj
      new.IVW.Model <- lm(BetaWj ~ -1 + Wj)
      new.EstimatesIVW <- summary(lm(BetaWj ~ -1 + Wj))
      new.IVW.Slope <- new.EstimatesIVW$coefficients[1]
      new.IVW.SE <- new.EstimatesIVW$coefficients[2]
      new.IVW_CI <- confint(new.IVW.Model)
      new.Qj <- W * (Ratios - new.IVW.Slope)^2
      new.Total_Q <- sum(new.Qj)
      new.Total_Q_chi <- pchisq(new.Total_Q, length(r_input[, 
                                                            2]) - 1, lower.tail = FALSE)
      Diff <- abs(Bhat1.Iterative - new.IVW.Slope)
      Bhat1.Iterative <- new.IVW.Slope
      Bhat1.SE <- new.IVW.SE
      Bhat1.t <- summary(new.IVW.Model)$coefficients[1, 
                                                     3]
      Bhat1.p <- summary(new.IVW.Model)$coefficients[1, 
                                                     4]
      count <- count + 1
    }
    It.Dat <- data.frame(Bhat1.Iterative, Bhat1.SE, Bhat1.t, 
                         Bhat1.p)
    multi_return2 <- function() {
      Out_list2 <- list(It.Res = It.Dat, count = count, 
                        It.CI = new.IVW_CI)
      return(Out_list2)
    }
    OUT2 <- multi_return2()
  }
  Bhat1.Iterative <- Iterative_ivw(tol)
  PL2 = function(a) {
    b = a[1]
    w = 1/((phi) * r_input[, 5]^2/r_input[, 2]^2 + (b^2) * 
             r_input[, 4]^2/r_input[, 2]^2)
    q = sum(w * (Ratios - b)^2)
  }
  PLfunc = function(a) {
    phi = a[1]
    PL2 = function(a) {
      beta = a[1]
      w = 1/(phi * r_input[, 5]^2/r_input[, 2]^2 + (beta^2) * 
               r_input[, 4]^2/r_input[, 2]^2)
      q = (sum(w * (Ratios - beta)^2))
    }
    b = optimize(PL2, interval = c(lb, ub))$minimum
    w = 1/(phi * r_input[, 5]^2/r_input[, 2]^2 + (b^2) * 
             r_input[, 4]^2/r_input[, 2]^2)
    q = (sum(w * (Ratios - b)^2) - DF)^2
  }
  BootVar = function(sims = 1000) {
    B = NULL
    pp = NULL
    for (hh in 1:sims) {
      L = length(r_input[, 2])
      choice = sample(seq(1, L), L, replace = TRUE)
      while(length(unique(choice))==1){
        choice = sample(seq(1, L), L, replace = TRUE)
      }
      bxg = r_input[, 2][choice]
      seX = r_input[, 4][choice]
      byg = r_input[, 3][choice]
      seY = r_input[, 5][choice]
      Ratios = byg/bxg
      W1 = 1/(seY^2/bxg^2)
      BIVw1 = Ratios * sqrt(W1)
      sW1 = sqrt(W1)
      IVWfitR1 = summary(lm(BIVw1 ~ -1 + sW1))
      phi_IVW1 = IVWfitR1$sigma^2
      W2 = 1/(seY^2/bxg^2 + (byg^2) * seX^2/bxg^4)
      BIVw2 = Ratios * sqrt(W2)
      sW2 = sqrt(W2)
      IVWfitR2 = summary(lm(BIVw2 ~ -1 + sW2))
      phi_IVW2 = IVWfitR2$sigma^2
      phi_IVW2 = max(1, phi_IVW2)
      phi_IVW1 = max(1, phi_IVW1)
      lb = IVWfitR1$coef[1] - 10 * IVWfitR1$coef[2]
      ub = IVWfitR1$coef[1] + 10 * IVWfitR1$coef[2]
      PL2 = function(a) {
        b = a[1]
        w = 1/((phi) * seY^2/bxg^2 + (b^2) * seX^2/bxg^2)
        q = sum(w * (Ratios - b)^2)
      }
      PLfunc = function(a) {
        phi = a[1]
        PL2 = function(a) {
          beta = a[1]
          w = 1/(phi * seY^2/bxg^2 + (beta^2) * seX^2/bxg^2)
          q = (sum(w * (Ratios - beta)^2))
        }
        b = optimize(PL2, interval = c(-lb, ub))$minimum
        w = 1/(phi * seY^2/bxg^2 + (b^2) * seX^2/bxg^2)
        q = (sum(w * (Ratios - b)^2) - DF)^2
      }
      phi = optimize(PLfunc, interval = c(phi_IVW2, phi_IVW1 + 
                                            0.001))$minimum
      B[hh] = optimize(PL2, interval = c(lb, ub))$minimum
      
    }
    se = sd(B)
    mB = mean(B)
    return(list(mB = mB, se = se))
  }
  CIfunc = function() {
    z = qt(df = DF, 0.975)
    z2 = 2 * (1 - pnorm(z))
    PL3 = function(a) {
      b = a[1]
      w = 1/(r_input[, 5]^2/r_input[, 2]^2 + (b^2) * r_input[, 
                                                             4]^2/r_input[, 2]^2)
      q = (sum(w * (Ratios - b)^2) - qchisq(1 - z2, DF))^2
    }
    lb = Bhat - 10 * SE
    ub = Bhat + 10 * SE
    low = optimize(PL3, interval = c(lb, Bhat))$minimum
    high = optimize(PL3, interval = c(Bhat, ub))$minimum
    CI = c(low, high)
    return(list(CI = CI))
  }
  phi = 1
  Bhat = optimize(PL2, interval = c(-2, 2))$minimum
  W = 1/(r_input[, 5]^2/r_input[, 2]^2 + (Bhat^2) * r_input[, 
                                                            4]^2/r_input[, 2]^2)
  SE = sqrt(1/sum(W))
  FCI = CIfunc()
  QIVW = sum(W * (Ratios - Bhat)^2)
  Qp = 1 - pchisq(QIVW, DF)
  Qind = W * (Ratios - Bhat)^2
  ExactQ = c(QIVW, Qp)
  ExactQind = Qind
  FE_EXACT = t(c(Bhat, SE, Bhat/SE, 2 * (1 - pt(abs(Bhat/SE), 
                                                DF))))
  FE_EXACT <- data.frame(FE_EXACT)
  names(FE_EXACT) <- c("Estimate", "Std.Error", "t value", 
                       "Pr(>|t|)")
  BIVW1 = Ratios * sqrt(1/(r_input[, 5]^2/r_input[, 2]^2))
  IVWfit1 = summary(lm(BIVW1 ~ -1 + sqrt(1/(r_input[, 5]^2/r_input[, 
                                                                   2]^2))))
  phi_IVW1 = IVWfit1$sigma^2
  BIVW2 <- Ratios * sqrt(1/(r_input[, 5]^2/r_input[, 2]^2 + 
                              (r_input[, 3]^2) * r_input[, 4]^2/r_input[, 2]^4))
  IVWfit2 = summary(lm(BIVW2 ~ -1 + sqrt(1/(r_input[, 5]^2/r_input[, 
                                                                   2]^2 + (r_input[, 3]^2) * r_input[, 4]^2/r_input[, 2]^4))))
  phi_IVW2 = IVWfit2$sigma^2
  phi_IVW2 = max(1, phi_IVW2)
  phi_IVW1 = max(1, phi_IVW1) + 0.001
  lb = Bhat - 10 * SE
  ub = Bhat + 10 * SE
  
  phi = optimize(PLfunc, interval = c(phi_IVW2, phi_IVW1))$minimum
  Bhat = optimize(PL2, interval = c(lb, ub))$minimum

  Boot = BootVar()
  SE = Boot$se
  RCI = Bhat + c(-1, 1) * qt(df = DF, 0.975) * SE
  RE_EXACT = t(c(Bhat, SE, Bhat/SE, 2 * (1 - pt(abs(Bhat/SE), 
                                                DF))))
  RE_EXACT <- data.frame(RE_EXACT)
  names(RE_EXACT) <- c("Estimate", "Std.Error", "t value", 
                       "Pr(>|t|)")
  Qj_Chi <- 0
  for (i in 1:length(Qj)) {
    Qj_Chi[i] <- pchisq(Qj[i], 1, lower.tail = FALSE)
  }
  r_input$Qj <- Qj
  r_input$Qj_Chi <- Qj_Chi
  Out_Indicator <- rep(0, length(r_input[, 2]))
  for (i in 1:length(r_input[, 2])) {
    if (Qj_Chi[i] < alpha) {
      Out_Indicator[i] <- 1
    }
  }
  r_input$Outliers <- factor(Out_Indicator)
  levels(r_input$Outliers)[levels(r_input$Outliers) == "0"] <- "Variant"
  levels(r_input$Outliers)[levels(r_input$Outliers) == "1"] <- "Outlier"
  if (sum(Out_Indicator == 0)) {
    outlier_status <- "No significant outliers"
    outtab <- "No significant outliers"
  }
  if (sum(Out_Indicator > 0)) {
    outlier_status <- "Outliers detected"
    Out_Dat <- subset(r_input, Outliers == "Outlier")
    outtab <- data.frame(Out_Dat[, 1], Out_Dat$Qj, Out_Dat$Qj_Chi)
    colnames(outtab) = c("SNP", "Q_statistic", "p.value")
  }
  if (summary == TRUE) {
    cat("\n")
    cat("Radial IVW\n")
    cat("\n")
    Sum.Dat <- data.frame(coef(EstimatesIVW))
    names(Sum.Dat) <- c("Estimate", "Std.Error", "t value", 
                        "Pr(>|t|)")
    names(Bhat1.Iterative$It.Res) <- names(Sum.Dat)
    combined.dat <- (rbind(Sum.Dat, Bhat1.Iterative$It.Res))
    combined.dat <- rbind(combined.dat, FE_EXACT)
    combined.dat <- rbind(combined.dat, RE_EXACT)
    row.names(combined.dat) <- c("Effect", "Iterative", "Exact (FE)", 
                                 "Exact (RE)")
    if (weights == 1) {
      row.names(combined.dat)[1] <- "Effect (1st)"
    }
    if (weights == 2) {
      row.names(combined.dat)[1] <- "Effect (2nd)"
    }
    if (weights == 3) {
      row.names(combined.dat)[1] <- "Effect (Mod.2nd)"
    }
    print(combined.dat)
    cat("\n")
    cat("\nResidual standard error:", round(EstimatesIVW$sigma, 
                                            3), "on", EstimatesIVW$df[2], "degrees of freedom")
    cat("\n")
    cat(paste(c("\nF-statistic:", " on", " and"), round(EstimatesIVW$fstatistic, 
                                                        2), collapse = ""), "DF, p-value:", format.pval(pf(EstimatesIVW$fstatistic[1L], 
                                                                                                           EstimatesIVW$fstatistic[2L], EstimatesIVW$fstatistic[3L], 
                                                                                                           lower.tail = FALSE), digits = 3))
    cat("\n")
    cat("Q-Statistic for heterogeneity:", Total_Q, "on", 
        length(r_input[, 2]) - 1, "DF", ",", "p-value:", 
        Total_Q_chi)
    cat("\n")
    cat("\n", outlier_status, "\n")
    cat("Number of iterations =", Bhat1.Iterative$count)
    cat("\n")
  }
  out_data <- data.frame(r_input[, 1], r_input[, 6], r_input[, 
                                                             7], r_input[, 8])
  out_data$Wj <- Wj
  out_data$BetaWj <- BetaWj
  out_data <- out_data[c(1, 5, 6, 2, 3, 4)]
  names(out_data) <- c("SNP", "Wj", "BetaWj", "Qj", "Qj_Chi", 
                       "Outliers")
  multi_return <- function() {
    Out_list <- list(coef = EstimatesIVW$coef, qstatistic = Total_Q, 
                     df = length(r_input[, 2]) - 1, outliers = outtab, 
                     data = out_data, confint = confint(IVW.Model), it.coef = combined.dat[2,]
                     , fe.coef = combined.dat[3, ], re.coef = combined.dat[4,], it.confint = Bhat1.Iterative$It.CI, fe.confint = FCI$CI, 
                     re.confint = RCI, meanF = mf,Total_Q_chi=Total_Q_chi)
    class(Out_list) <- "IVW"
    return(Out_list)
  }
  OUT <- multi_return()
}


PC.MR=function(trait.names=groups$food[1:6]          ## Trait names in the group to analise
               ,trait.file.name="../step3_MR_GWAS_correction/%%_corrected_GWAS.tsv.gz" ## Template for the file name %% is substituted by the trait name
               ,label="test"         # label for output
               ,p.limit=5e-8          # pvalue threshold
               ,cor2raw.tol=0.05      # 1+/- tollerance accepted to cor2raw ratio
               ,rg.matrix=rg.matrix   #genetic correlaiton matrix
               ,outcomes.list=outcomes ## OUTCOME list
               ,token="../../../f001_ciara_supervisormeeting/MR/Two_way_MR/Part_B/Z_scores_based/CAD/mrbase.oauth" # GOOGLE token 
               ,std=std   #vector with SD for each trait with trait names as names of vector
){
  
  
  g.idx=trait.names
  
  selection=c()
  single.traits=g.idx
  #selection=selection[,c("rs","chr","pos","corr.p")]
  for(j in single.traits){
    
    mv.res=fread(gsub("%%",j,trait.file.name),data.table=F)
    mv.res.sel=mv.res[which(mv.res$corr.p<p.limit 
                            & mv.res$corr2raw_ratio>1- cor2raw.tol
                            & mv.res$corr2raw_ratio<1+cor2raw.tol),]
    print(j)
    mv.res.sel=mv.res.sel[,c("rs","chrm","pos","corr.p")]
    names(mv.res.sel)=c("rs","chr","pos","corr.p")
    selection=rbind(selection,mv.res.sel)
    
  }
  selection=selection[order(selection$rs),]
  mappa=unique(selection[,c("rs","chr","pos")])
  mappa$p=NA
  for(j in 1:nrow(mappa)){
    
    mappa$p[j]=min(selection$corr.p[which(selection$rs==mappa$rs[j])])
    
    
  }
  
  prune=TwoSampleMR::clump_data(data.frame(SNP=mappa$rs,
                                           chr_name=mappa$chr,
                                           chrom_start=mappa$pos,
                                           pval.exposure =mappa$p))
  
  MHC=which(prune$chr_name==6 & prune$chrom_start>28000000 & prune$chrom_start<35000000)
  if(length(MHC)>0)prune=prune[-MHC,]
  # create intrument objects
  for(j in g.idx){
    
    tmp.res=fread(gsub("%%",j,trait.file.name),data.table=F)
    tmp.res=tmp.res[match(prune$SNP,tmp.res$rs),]
    tmp.res$beta1=tmp.res$beta1/std[j]
    tmp.res$corr.beta=tmp.res$corr.beta/std[j]
    tmp.res$se=tmp.res$se/std[j]
    if(j == g.idx[1]){
      
      instrument.unc=tmp.res[,c("rs","ref","alt","beta1")]
      names(instrument.unc)[ncol(instrument.unc)]=paste0("beta.",j)
      instrument.se=tmp.res[,c("rs","ref","alt","se")]
      names(instrument.se)[ncol(instrument.se)]=paste0("se.",j)
      instrument.corr=tmp.res[,c("rs","ref","alt","corr.beta")]
      names(instrument.corr)[ncol(instrument.corr)]=paste0("beta.",j)
      
      
    }else{
      
      ## Check alleles
      if(all(instrument.unc$ref==tmp.res$ref & instrument.unc$alt==tmp.res$alt)){
        instrument.unc$beta.new=tmp.res$beta1
        names(instrument.unc)[ncol(instrument.unc)]=paste0("beta.",j)
        instrument.se$se.new=tmp.res$se
        names(instrument.se)[ncol(instrument.se)]=paste0("se.",j)
        
        instrument.corr$beta.new=tmp.res$corr.beta
        names(instrument.corr)[ncol(instrument.corr)]=paste0("beta.",j)
      }
    }
  }
  
  ## Rotation matrix and loadings
  rotation.matrix=eigen(rg.matrix[g.idx,g.idx])
  explained.var=(rotation.matrix$values)/sum(rotation.matrix$values)
  colnames(instrument.unc)=gsub("beta.","",colnames(instrument.unc))
  colnames(instrument.corr)=gsub("beta.","",colnames(instrument.corr))
  
  library(psych)
  load=eigen.loadings(rotation.matrix)
  max.PC=min(which(cumsum(explained.var)>0.95))
  row.names(load)=g.idx
  load=load[,1:max.PC]
  colnames(load)=paste0("PC.",1:ncol(load))
  load=load[g.idx,]
  rotation.matrix=rotation.matrix$vectors
  row.names(rotation.matrix)=g.idx
  ## PCs.calc
  raw.val.unc=instrument.unc[,-c(1:3)]
  raw.val.unc=raw.val.unc[,row.names(rotation.matrix)]
  PCs=as.matrix(raw.val.unc)%*%(rotation.matrix)
  PCs=PCs[,1:max.PC]
  colnames(PCs)=paste0("PC.",1:ncol(PCs))
  instr.pcs.unc=data.frame(instrument.unc[,c(1:3)],PCs,stringsAsFactors = F)
  
  raw.val.corr=instrument.corr[,-c(1:3)]
  raw.val.corr=raw.val.corr[,row.names(rotation.matrix)]
  PCs=as.matrix(raw.val.corr)%*%(rotation.matrix)
  PCs=PCs[,1:max.PC]
  colnames(PCs)=paste0("PC.",1:ncol(PCs))
  instr.pcs.corr=data.frame(instrument.corr[,c(1:3)],PCs,stringsAsFactors = F)
  
  
  
  
  
  
  pdf(paste0(label,"_results_PCs.pdf"),width=14)
  res.tot=c()
  
  for(k in outcomes.list[,1]){
    t.id=strsplit(k,split="id:")[[1]][2]
    outcome.data=try(TwoSampleMR::extract_outcome_data(snps=instrument.unc$rs,proxies = F,access_token = "../../../f001_ciara_supervisormeeting/MR/Two_way_MR/Part_B/Z_scores_based/CAD/mrbase.oauth",outcomes = t.id))
    while(is(object =outcome.data,"try-error")){
      outcome.data=try(TwoSampleMR::extract_outcome_data(snps=instrument.unc$rs,proxies = F
                                                         ,access_token = token
                                                         ,outcomes = t.id))
    }
    
    if(!is.null(outcome.data)){
      outcome.data$effect_allele.outcome=toupper( outcome.data$effect_allele.outcome)
      outcome.data$other_allele.outcome=toupper(outcome.data$other_allele.outcome)
      outcome.data=outcome.data[which(outcome.data$pval.outcome>1e-5),]  
      ord.instr.unc=instr.pcs.unc[match(outcome.data$SNP,instr.pcs.unc$rs),]
      ord.instr.corr=instr.pcs.corr[match(outcome.data$SNP,instr.pcs.corr$rs),]
      
      if(nrow(ord.instr.unc)>(max.PC)*1.5){
        alleles=cbind(ord.instr.unc[,c("ref","alt")],outcome.data[,c("effect_allele.outcome","other_allele.outcome")])
        bad=apply(alleles,1,function(x)
          ifelse( length( which(x[1:2]%in%x[3:4]))==2,"ok","bad"))
        if(any(bad=="bad")){
          outcome.data=outcome.data[-which(bad=="bad"),]
          ord.instr.unc=ord.instr.unc[-which(bad=="bad"),]
          ord.instr.corr=ord.instr.corr[-which(bad=="bad"),]
          
        }
        to.flip=which(outcome.data$effect_allele.outcome==ord.instr.unc$ref & outcome.data$other_allele.outcome==ord.instr.unc$alt)
        if(length(to.flip)>0){
          
          outcome.data$beta.outcome[to.flip]=outcome.data$beta.outcome[to.flip]*-1
          
        }
        
        ## uncorrected
        ord.instr.unc$beta.outcome=outcome.data$beta.outcome
        ord.instr.unc$se.outcome=outcome.data$se.outcome
        row.names(ord.instr.unc)=ord.instr.unc$rs
        regress.unc=lm(ord.instr.unc$beta.outcome~as.matrix(ord.instr.unc[,grep("PC",names(ord.instr.unc))])-1,weights=ord.instr.unc$se.outcome^-2)
        
        beta.unc = summary(regress.unc)$coef[,1]
        se.unc=summary(regress.unc)$coef[,2]/
          min(1,summary(regress.unc)$sigma)
        pval.unc=pchisq((beta.unc/se.unc)^2,df=1,lower=F)
        names(beta.unc)=gsub('as.matrix(ord.instr.unc[, grep("PC", names(ord.instr.unc))])','',names(beta.unc),fixed=T)
        results.unc=data.frame(trait=names(beta.unc),beta.unc,se.unc,p=pval.unc,stringsAsFactors = F)
        
        ## corrected
        
        ord.instr.corr$beta.outcome=outcome.data$beta.outcome
        ord.instr.corr$se.outcome=outcome.data$se.outcome
        row.names(ord.instr.corr)=ord.instr.corr$rs
        regress.corr=lm(ord.instr.corr$beta.outcome~as.matrix(ord.instr.corr[,grep("PC",names(ord.instr.corr))])-1,weights=ord.instr.corr$se.outcome^-2)
        beta.corr = summary(regress.corr)$coef[,1]
        se.corr=summary(regress.corr)$coef[,2]/
          min(1,summary(regress.corr)$sigma)
        pval.corr=pchisq((beta.corr/se.corr)^2,df=1,lower=F)
        names(beta.corr)=gsub('as.matrix(ord.instr.corr[, grep("PC", names(ord.instr.corr))])','',names(beta.corr),fixed=T)
        results.corr=data.frame(trait=names(beta.corr),beta.corr,se.corr,p=pval.corr,stringsAsFactors = F)
        
        
        #directed.load=t((t(unclass(load))*()))
        directed.load=load[g.idx,1:max.PC]
        
        row.names(directed.load)=gsub("beta.","",row.names(directed.load))
        row.names(directed.load)=nice.labels(x = row.names(directed.load),lab.table = lab.table)
        a=ggcorrplot(directed.load,method="circle",lab = TRUE,lab_size = 2.5)
        
        forestplot.data.unc=data.frame(label=colnames(directed.load),beta=beta.unc,se=se.unc,stringsAsFactors = F)
        #forestplot.data$beta=forestplot.data$beta*sign(forestplot.data$beta)
        forestplot.data.unc$lower=forestplot.data.unc$beta+qnorm(0.025)*forestplot.data.unc$se
        forestplot.data.unc$upper=forestplot.data.unc$beta+qnorm(1-0.025)*forestplot.data.unc$se
        forestplot.data.unc$label=factor(forestplot.data.unc$label,levels = paste0("PC.",1:nrow(forestplot.data.unc)))
        forestplot.data.unc$type="Uncorrected"
        
        forestplot.data.corr=data.frame(label=colnames(directed.load),beta=beta.corr,se=se.corr,stringsAsFactors = F)
        #forestplot.data$beta=forestplot.data$beta*sign(forestplot.data$beta)
        forestplot.data.corr$lower=forestplot.data.corr$beta+qnorm(0.025)*forestplot.data.corr$se
        forestplot.data.corr$upper=forestplot.data.corr$beta+qnorm(1-0.025)*forestplot.data.corr$se
        forestplot.data.corr$label=factor(forestplot.data.corr$label,levels = paste0("PC.",1:nrow(forestplot.data.corr)))
        forestplot.data.corr$type="Corrected"
        
        forestplot.data=rbind(forestplot.data.corr,forestplot.data.unc)
        plab=paste(formatC(as.numeric(pval.unc),format="e",digits=1),formatC(as.numeric(pval.corr),format="e",digits=1),sep="\n")
        #plab=paste(formatC(as.numeric(pval),format="e",digits=1),sep="\n")
        
        f.plot=ggplot(data=forestplot.data,aes(x=label,y=beta,ymin=lower,ymax=upper,color=type))+
          geom_pointrange(position=position_dodge(width = 0.5))+coord_flip()+theme_classic()+
          scale_x_discrete(labels= plab)+
          geom_hline(yintercept=0, lty=2)+labs( x = "", y = "Effect size")+ 
          scale_colour_manual(values = c("Corrected" = "goldenrod", "Uncorrected" = "deepskyblue"))
        
        f.plot=f.plot+theme(plot.margin=unit(c(5,1,55,1),"pt"))
        forestplot.data$outcome=k
        a=a+coord_fixed(0.8)
        library(ggpubr)
        gigi=ggarrange(a,f.plot,ncol=2,align ="h")
        plot(annotate_figure(gigi,top=paste(i,outcome.data$outcome[1],sep="->")))
        forestplot.data$n.snp=nrow(ord.instr.unc)
        res.tot=rbind(res.tot,forestplot.data)
        
      }
    }
  }
  write.table(res.tot,paste0(label,"_PC_MR_results.tsv"),row.names=F,quote=F,sep="\t")
  dev.off()
  
}  


PC.mr.uni=function( trait.names="dr.fruit"          ## Trait names in the group to analise
                   ,trait.file.name="../step3_MR_GWAS_correction/%%_corrected_GWAS.tsv" ## Template for the file name %% is substituted by the trait name
                   ,label="test"         # label for output
                   ,p.limit=5e-8          # pvalue threshold
                   ,cor2raw.tol=0.05      # 1+/- tollerance accepted to cor2raw ratio
                   ,rg.matrix=rg.matrix
                   ,ph.matrix=ph.matrix#genetic correlaiton matrix
                   ,std=std
                   ,outcomes=outcomes
                   ,run.PC=2
                   ,mr.raps=TRUE
){
  
  ## Select IVs
  g.idx=trait.names
  if(length(g.idx)==1){
    run.PC=1
    #label=g.idx
  }
  
  selection=c()
  single.traits=g.idx
  #selection=selection[,c("rs","chr","pos","corr.p")]
  for(j in single.traits){
    
    mv.res=fread(gsub("%%",j,trait.file.name),data.table=F)
    mv.res.sel=mv.res[which(mv.res$corr.p<p.limit 
                            & mv.res$corr2raw_ratio>1- cor2raw.tol
                            & mv.res$corr2raw_ratio<1+cor2raw.tol),]
    print(j)
    mv.res.sel=mv.res.sel[,c("rs","chrm","pos","corr.p")]
    names(mv.res.sel)=c("rs","chr","pos","corr.p")
    selection=rbind(selection,mv.res.sel)
    
  }
  selection=selection[order(selection$rs),]
  mappa=unique(selection[,c("rs","chr","pos")])
  mappa$p=NA
  for(j in 1:nrow(mappa)){
    
    mappa$p[j]=min(selection$corr.p[which(selection$rs==mappa$rs[j])])
    
    
  }
  
  MHC=which(mappa$chr==6 & mappa$pos>28000000 & mappa$pos<35000000)
  if(length(MHC)>0)mappa=mappa[-MHC,]
  
  for(j in g.idx){
    
    tmp.res=fread(gsub("%%",j,trait.file.name),data.table=F)
    tmp.res=tmp.res[match(mappa$rs,tmp.res$rs),]
    tmp.res$beta1=tmp.res$beta1/std[j]
    tmp.res$corr.beta=tmp.res$corr.beta/std[j]
    tmp.res$se=tmp.res$se/std[j]
    if(j == g.idx[1]){
      
      instrument.unc=tmp.res[,c("rs","ref","alt","freq1","beta1")]
      names(instrument.unc)[ncol(instrument.unc)]=paste0("beta.",j)
      instrument.se=tmp.res[,c("rs","ref","alt","se")]
      names(instrument.se)[ncol(instrument.se)]=paste0("se.",j)
      
      instrument.corr=tmp.res[,c("rs","ref","alt","freq1","corr.beta")]
      names(instrument.corr)[ncol(instrument.corr)]=paste0("beta.",j)
      
      
    }else{
      
      ## Check alleles
      if(all(instrument.unc$ref==tmp.res$ref & instrument.unc$alt==tmp.res$alt)){
        instrument.unc$beta.new=tmp.res$beta1
        names(instrument.unc)[ncol(instrument.unc)]=paste0("beta.",j)
        instrument.se$se.new=tmp.res$se
        names(instrument.se)[ncol(instrument.se)]=paste0("se.",j)
        
        instrument.corr$beta.new=tmp.res$corr.beta
        names(instrument.corr)[ncol(instrument.corr)]=paste0("beta.",j)
      }
    }
  }
  
  
  
  
  
  
  
  # create intrument objects
  
  if(length(g.idx)>1){
    rotation.matrix=eigen(rg.matrix[g.idx,g.idx])
    explained.var=(rotation.matrix$values)/sum(rotation.matrix$values)
    colnames(instrument.unc)=gsub("beta.","",colnames(instrument.unc))
    colnames(instrument.corr)=gsub("beta.","",colnames(instrument.corr))
    colnames(instrument.se)=gsub("se.","",colnames(instrument.se))
    
    library(psych)
    load=eigen.loadings(rotation.matrix)
    max.PC=min(which(cumsum(explained.var)>0.95))
    row.names(load)=g.idx
    load=load[,1:max.PC]
    colnames(load)=paste0("PC.",1:ncol(load))
    load=load[g.idx,]
    write.table(load,file=paste0(label,"_loadings.tsv"),sep="\t")
    rotation.matrix=rotation.matrix$vectors
    row.names(rotation.matrix)=g.idx
    ## PCs.calc
    raw.val.unc=instrument.unc[,-c(1:4)]
    raw.val.unc=raw.val.unc[,row.names(rotation.matrix)]
    PCs=as.matrix(raw.val.unc)%*%(rotation.matrix)
    PCs=PCs[,1:max.PC]
    colnames(PCs)=paste0("PC.",1:ncol(PCs))
    instr.pcs.unc=data.frame(instrument.unc[,c(1:4)],PCs,stringsAsFactors = F)
    
    raw.val.corr=instrument.corr[,-c(1:4)]
    raw.val.corr=raw.val.corr[,row.names(rotation.matrix)]
    PCs=as.matrix(raw.val.corr)%*%(rotation.matrix)
    PCs=PCs[,1:max.PC]
    colnames(PCs)=paste0("PC.",1:ncol(PCs))
    instr.pcs.corr=data.frame(instrument.corr[,c(1:4)],PCs,stringsAsFactors = F)
    
    
    ### se calculation
    red.ph.matr=ph.matrix[g.idx,g.idx]
    se.matrix=matrix(ncol=ncol(rotation.matrix),nrow=nrow(instrument.se))
    for(i in 1:ncol(rotation.matrix)){
      
      vec=t(t(instrument.se[,g.idx])*rotation.matrix[,i])
      se.vec=rep(NA,length=nrow(instrument.se))
      for(j in 1:nrow(instrument.se)) se.vec[j]=sqrt(vec[j,]  %*%    red.ph.matr %*% (matrix(vec[j,])))
      se.matrix[,i]=se.vec
    }
    instr.pcs.se=data.frame(instr.pcs.corr$rs,se.matrix[,1:max.PC],stringsAsFactors=F)
    row.names(instr.pcs.se)=instr.pcs.corr$rs
    colnames(instr.pcs.se)=c("rs",paste0("PC.",1:max.PC))
    
    
    library(reshape2)
    
    instr.pcs.unc.melted=melt(instr.pcs.unc,id.vars = c("rs", "ref", "alt" ,   "freq1"))
    instr.pcs.corr.melted=melt(instr.pcs.corr,id.vars = c("rs", "ref", "alt" ,   "freq1"))
    instr.pcs.se=melt(instr.pcs.se)
    names(instr.pcs.se)[3]="se"
    names(instr.pcs.unc.melted)[6]="beta"
    names(instr.pcs.corr.melted)[6]="beta"
    
    instr.pcs.unc.melted=merge(instr.pcs.unc.melted,instr.pcs.se,by=c("rs", "variable"))
    instr.pcs.corr.melted=merge(instr.pcs.corr.melted,instr.pcs.se,by=c("rs", "variable"))
    
    instr.pcs.unc.melted$p=pchisq((instr.pcs.unc.melted$beta/instr.pcs.unc.melted$se)^2,df=1,lower=F)
    instr.pcs.corr.melted$p=pchisq((instr.pcs.corr.melted$beta/instr.pcs.corr.melted$se)^2,df=1,lower=F)
    
    unc.instr=format_data(instr.pcs.unc.melted,type = "exposure",phenotype_col = "variable",snp_col = "rs"
                          ,beta_col = "beta",se_col = "se",effect_allele_col = "alt",other_allele_col = "ref",eaf_col = "freq1"
                          ,pval_col = "p")
    
    corr.instr=format_data(instr.pcs.corr.melted,type = "exposure",phenotype_col = "variable",snp_col = "rs"
                           ,beta_col = "beta",se_col = "se",effect_allele_col = "alt",other_allele_col = "ref",eaf_col = "freq1"
                           ,pval_col = "p")
    
    
    
  }else{
    
    instrument.corr=merge(instrument.corr,instrument.se,by=c("rs", "ref", "alt"))
    instrument.corr$variable=g.idx
    instrument.unc=merge(instrument.unc,instrument.se,by=c("rs", "ref", "alt"))
    instrument.unc$variable=g.idx
    names(instrument.corr)=gsub(paste0(".",g.idx),"",names(instrument.corr))
    names(instrument.unc)=gsub(paste0(".",g.idx),"",names(instrument.unc))
    instrument.corr$p=pchisq((instrument.corr$beta/instrument.corr$se)^2,df=1,lower=F)
    instrument.unc$p=pchisq((instrument.unc$beta/instrument.unc$se)^2,df=1,lower=F)
    
    unc.instr=format_data(instrument.unc,type = "exposure",phenotype_col = "variable",snp_col = "rs"
                          ,beta_col = "beta",se_col = "se",effect_allele_col = "alt",other_allele_col = "ref",eaf_col = "freq1"
                          ,pval_col = "p")
    
    corr.instr=format_data(instrument.corr,type = "exposure",phenotype_col = "variable",snp_col = "rs"
                           ,beta_col = "beta",se_col = "se",effect_allele_col = "alt",other_allele_col = "ref",eaf_col = "freq1"
                           ,pval_col = "p")
    
    
  }
  
  
  res.pc.unc=c()
  res.pc.corr=c()
  token = "../../../f001_ciara_supervisormeeting/MR/Two_way_MR/Part_B/Z_scores_based/CAD/mrbase.oauth"
  #ao=available_outcomes(token)
  
  for(k in outcomes$id){
    if(length(g.idx)==1){
      run.trait=g.idx
    }else{
      run.trait=paste0("PC.",1:run.PC)
    }
    instrument.unc.tmp=unc.instr[unc.instr$exposure==run.trait[1],]
    outcome.data=c()
    splits=ceiling(nrow(instrument.unc.tmp)/400)
    index=rep_len(1:splits, length.out=nrow(instrument.unc.tmp))
    for(sp in 1:splits){
      outcome.data.tmp=try(extract_outcome_data(snps=instrument.unc.tmp$SNP[index==sp],outcomes = k,access_token = token,proxies = FALSE))
      while(is(outcome.data.tmp,class2 = "try-error")){
        outcome.data.tmp=try(extract_outcome_data(snps=instrument.unc.tmp$SNP[index==sp],outcomes = k,access_token = token,proxies = FALSE))
      }
      outcome.data=rbind(outcome.data,outcome.data.tmp)
    }
    instrument.unc=instrument.unc[instrument.unc$SNP%in%outcome.data$SNP,]
    if(!is.null(outcome.data) & nrow(instrument.unc.tmp)>0){
      mappa.tmp=unique(mappa[,1:3])[match(instrument.unc.tmp$SNP,mappa$rs),]
      
      prune=try(TwoSampleMR::clump_data(data.frame(SNP=mappa.tmp$rs,    #### Continua qui 
                                                   chr_name=mappa.tmp$chr,
                                                   chrom_start=mappa.tmp$pos,
                                                   pval.exposure =instrument.unc.tmp$pval.exposure)))
      while(is(prune,class2 = "try-error")){
        prune=try(TwoSampleMR::clump_data(data.frame(SNP=mappa.tmp$rs,    #### Continua qui 
                                                     chr_name=mappa.tmp$chr,
                                                     chrom_start=mappa.tmp$pos,
                                                     pval.exposure =instrument.unc.tmp$pval.exposure)))
        
      }
      
    for(p in run.trait){
         unc.instr.pruned=unc.instr[which(unc.instr$SNP%in%prune$SNP & unc.instr$exposure==p),] 
         corr.instr.pruned=corr.instr[which(corr.instr$SNP%in%prune$SNP & corr.instr$exposure==p),] 
         
       ## MR uncorrected betas

          harm=harmonise_data(unc.instr.pruned,outcome.data,action = 1)
          harm=harm[which(harm$beta.outcome!=0 & harm$pval.outcome>1e-5),]
  
          if(nrow(harm)>0){
            res.tmp=mr.run(harm,mr.raps=mr.raps)
            res.pc.unc=rbind(res.pc.unc,res.tmp)
          }
          
        ## MR corrected betas
          
          
          harm=harmonise_data(corr.instr.pruned,outcome.data,action = 1)
          harm=harm[which(harm$beta.outcome!=0 & harm$pval.outcome>1e-5),]
          
          if(nrow(harm)>0){
            res.tmp=mr.run(harm,mr.raps=mr.raps)
            res.pc.corr=rbind(res.pc.corr,res.tmp)
          }
          
        }
      }
  }
  
  
  res.pc.unc$type="Uncorrected"
  res.pc.unc$exposure=paste0(label,"_",res.pc.unc$exposure)
  res.pc.corr$type="Corrected"
  res.pc.corr$exposure=paste0(label,"_",res.pc.corr$exposure)
  
  res.tot=rbind(res.pc.unc,res.pc.corr)
  write.table(res.tot,file=paste0(label,"_PC_MR_results.tsv"),row.names=F,quote=F,sep="\t")
  
  
}






