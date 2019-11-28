



res.tot=c()

for(j in 5:10){
n=j
for(i in 1:1000){
 prop.fat.wt=rnorm(n=n,0.2,0.1)
 prop.prot.wt=rnorm(n=n,0.5,0.1)
 prop.norm.wt=1-(prop.fat.wt+prop.prot.wt)

 prop.fat.mc4r=rnorm(n=n,0.3,0.1)
 prop.prot.mc4r=rnorm(n=n,0.3,0.1)
 prop.norm.mc4r=1-(prop.fat.mc4r+prop.prot.mc4r)


 wildtype=rnorm(n =n ,mean = 3,sd = 0.1)
 mc4r=rnorm(n=n,mean = 4,sd = 0.1)

 data.wt=data.frame(type="wt",overall.cons=wildtype,fat=wildtype*prop.fat.wt,prot=wildtype*prop.prot.wt,norm=wildtype*prop.norm.wt)
 data.mc4r=data.frame(type="mc4r",overall.cons=mc4r,fat=mc4r*prop.fat.mc4r,prot=mc4r*prop.prot.mc4r,norm=mc4r*prop.norm.mc4r)

 data.all=rbind(data.wt,data.mc4r)
 res=summary(aov(manova(cbind(fat/overall.cons,prot/overall.cons,norm/overall.cons)~type,data=data.all)))
 res.tmp=c(n,res$` Response 1`$`Pr(>F)`[1],res$` Response 2`$`Pr(>F)`[1],res$` Response 3`$`Pr(>F)`[1])
 res.tot=rbind(res.tot,res.tmp)
}
}



alpha=0.05

by(res.tot[,-1]<alpha,INDICES = res.tot[,1],function(x)colSums(x)/1000)
  











