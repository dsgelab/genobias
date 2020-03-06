## Set up a conda enviroment
# use Anaconda3
# conda create --name myenv python=3.7
# source activate myenv
# conda install -c conda-forge r-rgdal

install.packages(c("data.table","dplyr","survey","SDMTools","ggplot2","rgdal","sp"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("preprocessCore")

## READ RAW DATA AND CREATE PHENOTYPE ##
library(data.table)
library(dplyr)
library(survey)
library(SDMTools)
library(ggplot2)
library(broom)

# ! Access to UK Biobank individual-level data is required, we can only provide the variable names here.
d <- fread('path/ukbb_raw_data', header=T)
colnames(d) <- paste0("f.",gsub("-",".",colnames(d)))

t1 <- c("f.eid","f.31.0.0","f.53.0.0","f.6138.0.0","f.22704.0.0","f.22702.0.0","f.21022.0.0","f.54.0.0","f.22700.0.0","f.21001.0.0","f.2040.0.0",colnames(d)[grepl("f.22704.0.",colnames(d))],colnames(d)[grepl("f.22702.0.",colnames(d))],colnames(d)[grepl("f.22700.0.",colnames(d))])

bdE4 <- d[,t1, with=FALSE]

bdE4$edu <- ifelse(bdE4$f.6138.0.0==1,6,
ifelse(bdE4$f.6138.0.0==2,5,
ifelse(bdE4$f.6138.0.0==3,2.5,
ifelse(bdE4$f.6138.0.0==4,2.5,
ifelse(bdE4$f.6138.0.0==5,6,
ifelse(bdE4$f.6138.0.0==6,6,
ifelse(bdE4$f.6138.0.0==-7,1, 
ifelse(bdE4$f.6138.0.0==-3,NA,NA))))))))

bdE4$age_2011 <- as.numeric(bdE4$f.21022.0.0 + ((as.Date("2011-03-27") - as.Date(bdE4$f.53.0.0,  format="%Y-%m-%d"))/365.25))

bdE4$age_2011_cat <- ifelse(bdE4$age_2011 <= 49,"35-49",
                     ifelse(bdE4$age_2011 <= 64,"50-64","65+"))

tt <- bdE4[,grepl("f.22700.0.",colnames(bdE4)), with=FALSE]
tt2 <- tt[ , lapply( .SD, as.Date , format="%Y-%m-%d") ]
tt2d <- as.data.frame(tt2)

tt_east <- bdE4[,grepl("f.22702.0.",colnames(bdE4)), with=FALSE]
tt_north <- bdE4[,grepl("f.22704.0.",colnames(bdE4)), with=FALSE]

coord_2011_east <- NULL
coord_2011_north <- NULL
for (i in 1:nrow(tt2d))
{
  x <- tt2d[i,]
  closest_2011 <- which.min(as.Date("2011-03-27")-as.Date(as.numeric(x),origin="1970-01-01"))
  
  coord_2011_east <- rbind(coord_2011_east,c(bdE4$f.eid[i],as.numeric(tt_east[i,closest_2011, with=FALSE])))
  coord_2011_north <- rbind(coord_2011_north,c(bdE4$f.eid[i],as.numeric(tt_north[i,closest_2011, with=FALSE])))
  print(i)
}
coord_2011_eastdf <- data.frame(coord_2011_east)
colnames(coord_2011_eastdf) <- c("f.eid","east")

coord_2011_northdf <- data.frame(coord_2011_north)
colnames(coord_2011_northdf) <- c("f.eid","north")

bdE5 <- merge(bdE4[,!grepl("f.22702.0.|f.22704.0.|f.22700.0.",colnames(bdE4)),with=FALSE],coord_2011_eastdf, by="f.eid")
bdE6 <- merge(bdE5,coord_2011_northdf, by="f.eid")

write.table(bdE6,file="pheno.tsv", row.names=F, col.names=T, quote=F, sep="\t")

### ASSIGN AREA TO EACH INDIVIDUAL ###
pheno <- read.table("pheno.tsv", header=T)
pheno = pheno[!is.na(pheno$east) & !is.na(pheno$north),]
phenob = pheno %>% select(f.eid, east, north)

sp::coordinates(phenob) = ~east+north
# Reading 2011 Census Geography boundaries
# infuse_msoa_lyr_2011_clipped.zip available here https://ec2-54-77-102-218.eu-west-1.compute.amazonaws.com/dataset/2011-census-geography-boundaries-middle-layer-super-output-areas-and-intermediate-zones-7
ogr = rgdal::readOGR('path/infuse_msoa_lyr_2011_clipped/')
shapefile = sp::spTransform(ogr, sp::CRS("+init=epsg:27700"))
phenob@proj4string = shapefile@proj4string
pisb = sp::over(phenob, shapefile)

pheno_final <- data.frame(pheno,pisb)

# Remove non euro individuals (requires UKBB individual-level data)
remove <- read.table('path/non_euro_ID_list', header=F)


# Keep infam
infam <- read.table('path/ukbb_fam_file', header=F)
pheno_final <- pheno_final[!pheno_final$f.eid %in% c(remove$V1) & pheno_final$f.eid %in% c(infam$V1) ,]

# Some check about if annotations are correct
write.table(data.frame(table(pheno_final$geo_label,pheno_final$f.54.0.0)),file="test.csv", quote=F, row.names=F,col.names=T, sep="\t")

# Keep geocodes that are in the census and remove missing education values
# This file is available in genobias/census/
census <- read.csv("genobias/census/MSOA.EA_census.sex.age.csv")

pheno_finalC <- pheno_final[as.character(pheno_final$geo_code) %in% as.character(unique(census$geo_code)) & !is.na(pheno_final$geo_code),]
pheno_finalC <- pheno_finalC[!is.na(pheno_finalC$edu),]

length(unique(pheno_finalC$geo_code))

# Test pheno_final C differences between males and females
t.test(pheno_finalC$edu[pheno_finalC$f.31.0.0==0],pheno_finalC$edu[pheno_finalC$f.31.0.0==1])


### FINAL:UKBB BY SEX AND BY AGE ###
ukbb_new <- data.frame(aggregate(pheno_finalC$edu,list(pheno_finalC$geo_code,pheno_finalC$f.31.0.0,pheno_finalC$age_2011_cat),mean))
colnames(ukbb_new) <- c("geo_code","sex","age","edu")

ukbb_newN <- data.frame(aggregate(pheno_finalC$edu,list(pheno_finalC$geo_code,pheno_finalC$f.31.0.0,pheno_finalC$age_2011_cat),length))
colnames(ukbb_newN) <- c("geo_code","sex","age","eduN")

ukbb_new <- data.frame(ukbb_new,eduN=ukbb_newN$eduN)

ukbb_design <- svydesign(id=~1, strata=~sex + age + geo_code, weights = ~eduN, variable=~edu+sex , data=ukbb_new, nest=TRUE)

svymean(~edu,ukbb_design)
svyby(~edu, ~sex, ukbb_design, svymean)


census <- read.csv("genobias/census/MSOA.EA_census.sex.age.csv")
census <- census[census$geo_code %in% ukbb_new$geo_code,]
censusN <- read.csv("genobias/census/MSOA.EA_census.sex.age.including_frq.csv")
censusN <- censusN[,c("geo_code","census_EA.MSOA.m.35_49.frq","census_EA.MSOA.f.35_49.frq","census_EA.MSOA.m.50_64.frq","census_EA.MSOA.f.50_64.frq","census_EA.MSOA.m.65_plus.frq","census_EA.MSOA.f.65_plus.frq")]
colnames(censusN) <- c("geo_code","census_EA.MSOA.m.35_49","census_EA.MSOA.f.35_49","census_EA.MSOA.m.50_64","census_EA.MSOA.f.50_64","census_EA.MSOA.m.65_plus","census_EA.MSOA.f.65_plus")

census_melt <- melt(census[,c("geo_code","census_EA.MSOA.m.35_49","census_EA.MSOA.f.35_49","census_EA.MSOA.m.50_64","census_EA.MSOA.f.50_64","census_EA.MSOA.m.65_plus","census_EA.MSOA.f.65_plus")])
census_melt$sex <- ifelse(grepl(".m.",census_melt$variable),1,0)
census_melt$age <- ifelse(grepl("35_49",census_melt$variable),"35_49",
                          ifelse(grepl("50_64",census_melt$variable),"50_64","65+"))

census_meltN <- melt(censusN[,c("geo_code","census_EA.MSOA.m.35_49","census_EA.MSOA.f.35_49","census_EA.MSOA.m.50_64","census_EA.MSOA.f.50_64","census_EA.MSOA.m.65_plus","census_EA.MSOA.f.65_plus")])

census_melt <- merge(census_melt,census_meltN,by=c("geo_code","variable"))

census_design <- svydesign(id=~1, strata=~sex + age + geo_code, weights=~value.y, variable=~value.x+sex , data=census_melt, nest=TRUE)

svymean(~value.x,census_design)
svyby(~value.x, ~sex, census_design, svymean)


syukb <- svyby(~edu, ~sex, ukbb_design, svymean)
sycens <- svyby(~value.x, ~sex, census_design, svymean)

syukb$edu/sycens$value.x

df <- data.frame(edu = c(syukb$edu,sycens$value.x), sex=c(syukb$sex,sycens$sex), eduse=c(syukb$se,NA,NA), group=c("UKBB","UKBB","Census","Census"))
df$sexlab <- ifelse(df$sex==1,"Males","Females")


#       edu sex      eduse  group  sexlab
#1 3.951293   0 0.01055599   UKBB Females
#2 4.131304   1 0.01002236   UKBB   Males
#3 3.200241   0         NA Census Females
#4 3.573665   1         NA Census   Males


pdf("genobias/census/figure.pdf", height=4,width=4)
ggplot(df, aes(x=sexlab, y=edu, col=group)) + 
  geom_errorbar(aes(ymin=edu-(1.96*eduse), ymax=edu+(1.96*eduse)), width=.1) +
  geom_point(size=3) +
  ylab("Educational attainment level") +
  xlab(NULL) +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none") +
  theme_light() + scale_color_manual("",values=c("red","blue"),labels=c("Census","UK Biobank")) + guides(col = guide_legend(reverse = TRUE))
dev.off()


### PLOT MAP ON ENGLAND WITH COORDINATES ###

library(broom)
phenob2 = pheno %>% select(f.eid, east, north)
phenob2_sub <- phenob2[phenob2$f.eid %in% pheno_finalC$f.eid,]
sp::coordinates(phenob2_sub) = ~east+north
ogr = rgdal::readOGR('path/infuse_msoa_lyr_2011_clipped/')
england_wales <- tidy(ogr)
geo_codedf <- data.frame(geo_code=as.character(ogr$geo_code),id=as.character(0:(length(ogr$geo_code)-1)),stringsAsFactors=F)

england_walesm <- merge(england_wales,geo_codedf,by="id")
england_walesm$ind <- factor(ifelse(england_walesm$geo_code %in% unique(pheno_finalC$geo_code),1,0))

# Plot it
pdf("genobias/census/map.pdf", height=16,width=16)
ggplot() +
  geom_polygon(data = england_walesm, aes( x = long, y = lat, group = group,  fill=ind), color="black", size=0.05) +
  theme_void() + geom_point(aes(x=phenob2_sub@coords[,1],y=phenob2_sub@coords[,2]),col="blue",size=0.01) + scale_fill_manual("",values=c("white","red"), guide=FALSE)
dev.off()



