setwd("~/Dropbox/gp_2022")

library(dartR)
library(adegenet)
library(pegas)
library(vcfR)
library(dplyr)
library(hierfstat)
library(poppr)

#note I deleted individual GOP794 from the vcf file: 
#vcf_GP_filtered_strict_complete.recodefixed.recode.vcf
#as it was mislabelled

#read vcf file
vcf <- read.vcfR("./data/vcf_GP_filtered_strict_complete.recodefixed.recode1_vcf.recode.vcf")

#convert vcf to genlight
gop_gl <- vcfR2genlight(vcf)

#read in covariate file
cov <- read.csv("./data/gp_database_covariates_Aug2021.csv")

#make covariate file only include individuals that are retained in the vcf file.
#important to make sure the order of individuals in ind.metrics matches the order of individuals in main
#genlight object
cov_reduced <- cov[match(gop_gl$ind.names, cov$id), ]

#add in all the ind.metrics to dataframe

gop_gl@other$ind.metrics <- data

gop_gl@other$ind.metrics <- as.data.frame(matrix(nrow=nrow(cov_reduced),ncol=10))
colnames(gop_gl@other$ind.metrics) <- c( "id", "pop", "sex", "length", 
                                             "date", "oto_age", "natal_region", "stocked", "lat", "long")
gop_gl$pop <- as.factor(cov_reduced$pop)
gop_gl$other$ind.metrics$id <- cov_reduced$id
gop_gl$other$ind.metrics$pop <- as.factor(cov_reduced$pop)
gop_gl$other$ind.metrics$sex <- as.factor(cov_reduced$sex.gonads)
gop_gl$other$ind.metrics$length <- cov_reduced$length
gop_gl$other$ind.metrics$date <- cov_reduced$date_capture
gop_gl$other$ind.metrics$oto_age <- cov_reduced$oto_age
gop_gl$other$ind.metrics$natal_region <-  as.factor(cov_reduced$natal_region)
gop_gl$other$ind.metrics$stocked <-  as.factor(cov_reduced$stocked)
gop_gl$other$ind.metrics$lat <-  cov_reduced$lat
gop_gl$other$ind.metrics$long <-  cov_reduced$long

gop_gl$other$loc.metrics <- as.data.frame(matrix(nrow=(length(vcf@fix[,1]))), ncol=1)
gop_gl$other$loc.metrics$AlleleID <- gop_gl$loc.names
gop_gl$other$loc.metrics$Chrom <- as.factor(vcf@fix[,1])
gop_gl$other$loc.metrics <- gop_gl$other$loc.metrics[,-1]
gop_gl$chromosome <- as.factor(vcf@fix[,1])

#check locnames and loc.metrics match
all.equal(colnames(as.matrix(gop_gl)), gop_gl@other$loc.metrics$AlleleID)

#fix pop codes
gop_gl$other$ind.metrics$pop <- factor(
  gop_gl$other$ind.metrics$pop,
  levels = c("BROKEN", "BROKEN ", "CAMPASPE", "DUMARESQ", "EDWARD-WAKOOL", 
             "GINGHAM", "GOULBURN", "GUNBOWER", "GWYDIR", "L_MURRAY",
             "LODDON_PYRAMID", "LOWER_DARLING", "MACINTYRE", "MACQUARIE",
             "MULLAROO_LINDSAY",  "MURRAY-Barmah", "MURRAY-dsYarrawonga",
             "MURRAY-mid",  "MURRUMBIDGEE",  "OVENS", "SEVERN",   "UPPER_DARLING",  "WIMMERA"),
  labels = c("BROKEN", "BROKEN", "CAMPASPE", "DUMARESQ", "EDWARD-WAKOOL", 
             "GINGHAM", "GOULBURN", "GUNBOWER", "GWYDIR", "L_MURRAY",
             "LODDON_PYRAMID", "LOWER_DARLING", "MACINTYRE", "MACQUARIE",
             "MULLAROO_LINDSAY",  "MURRAY-Barmah", "MURRAY-dsYarrawonga",
             "MURRAY-mid",  "MURRUMBIDGEE",  "OVENS", "SEVERN",   "UPPER_DARLING",  "WIMMERA")
)

gop_gl@pop <- gop_gl$other$ind.metrics$pop

#had to convert to gi and back again to get the dartR to recognise the format
gop_gi <- gl2gi(gop_gl)
#convert back to gl
gop_gl_1 <- gi2gl(gop_gi)

#save(gop_gl_1, file="./data/gop_gl_1.rdata")