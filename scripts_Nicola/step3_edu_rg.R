
file.list.edu=system("ls simul/ldhub/edu*.gz",intern=T)
file.list.ganna=system("ls simul/files_ganna/*.gz",intern=T)


for(i in file.list.ganna){
    tratto2=gsub(".LDSC.sumstats.gz","",gsub("simul/files_ganna/","",i))
    for(j in file.list.edu){
      tratto1=gsub("_mung.sumstats.gz","",gsub("simul/ldhub/","",j))
      
      
      system(paste0("ldsc.py --rg ",i,",",j," --ref-ld-chr /opt/working/wilson/app_full_stuff/ldsc/eur_w_ld_chr/ --w-ld-chr /opt/working/wilson/app_full_stuff/ldsc/eur_w_ld_chr/  --out simul/rgs/rg_",tratto1,"_",tratto2,".out"))

      
    }
}

