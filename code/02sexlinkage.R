# install the package
#devtools::install_github('jdyen/dartR@sexlinkage', dependencies = TRUE)
#manually ran the function gl.sexlinkage

# load the package
library(dartR)
 
#load data
load("./data/gop_gl_1.rdata")

#need to replace blanks with NA in gop_gl_1@other$ind.metrics$sex
gop_gl_1@other$ind.metrics$sex[gop_gl_1@other$ind.metrics$sex == ""] <- NA 
gop_gl_1@other$ind.metrics$sex <- as.factor(gop_gl_1@other$ind.metrics$sex)

#t.het = is the proportion of heterozygotes allowed in the sex that is expected to be homozygous
#t.hom = is the proportion of homozygotes allowed in the sex that is expected to be heterozygous
#t.abs = relevant to Y-linkage test - is the threshold for mis-sexed individuals; t.abs = 1 means a
#single individual in a given sex is treated as that sex being absent
sex_test <- gl.sexlinkage(gop_gl_1, t.het = 0.1, t.hom = 0.1, t.abs = 5,
                                  verbose = FALSE, na.rm = TRUE)

#test for xy homologs
#finds no loci consistent with XY fixed differences or Y linkage
print(sex_test, system = 'xy')
summary(sex_test, include = c('locus', 'count'), system = 'xy')


#increase t.hom threshold to get X-linked loci
sex_test1 <- gl.sexlinkage(gop_gl_1, t.het = 0.01, t.hom = 0.95, t.abs = 1,
                                  verbose = FALSE, na.rm = TRUE)

#print results
print(sex_test1, system = 'xy')
#finds 57 loci consistent with X-linkage
summary(sex_test1, include = c('locus', 'count'), system = 'xy')

#get list of X-linked loci
list_Xlinked <- sex_test1$fem_het_male_hom
#write.csv(list_Xlinked, "./outputs/list_Xlinked_summary_counts.csv")

#keep loci with dartR evidence of sex-linkage
dartR_sex_Xlinked_loci <- as.factor(rownames(list_Xlinked))
dartR_sex_Xlinked_loci_to_keep <-  match(dartR_sex_Xlinked_loci, as.factor(locNames(gop_gl_1)))
gl_gp_dartR_sex_Xlinked_loci_to_keep <- gop_gl_1[, dartR_sex_Xlinked_loci_to_keep]
gl_gp_dartR_sex_Xlinked_loci_to_keep@other$loc.metrics <- gl_gp_dartR_sex_Xlinked_loci_to_keep@other$loc.metrics[dartR_sex_Xlinked_loci_to_keep, ]
#check locnames and locmetrics match
all.equal(colnames(as.matrix(gl_gp_dartR_sex_Xlinked_loci_to_keep)), 
          gl_gp_dartR_sex_Xlinked_loci_to_keep@other$loc.metrics$AlleleID)
#write.csv(gl_gp_dartR_sex_Xlinked_loci_to_keep, "./data/dartRxlinked.csv")

#save(gl_gp_dartR_sex_Xlinked_loci_to_keep, file="./data/gl_gp_dartR_sex_Xlinked_loci_to_keep.rdata")

#basic pca
pc <- gl.pcoa(gl_gp_dartR_sex_Xlinked_loci_to_keep, nfactors=3)
gl_gp_dartR_sex_Xlinked_loci_to_keep@pop <- gl_gp_dartR_sex_Xlinked_loci_to_keep@other$ind.metrics$sex
gl.pcoa.plot(pc, gl_gp_dartR_sex_Xlinked_loci_to_keep, labels="pop", xaxis=1, yaxis=2)

gi_gp_dartR_sex_Xlinked_loci_to_keep <- gl2gi(gl_gp_dartR_sex_Xlinked_loci_to_keep)
#minor allele freq
maf_gop <- minorAllele(gi_gp_dartR_sex_Xlinked_loci_to_keep)
hist(maf_gop)
#per locus heterozygosity
locusstats_gop <- locus_table(gi_gp_dartR_sex_Xlinked_loci_to_keep)
hist(locusstats_gop[,3])

#remove the 57 putative X-linked loci identified by dartR
dartR_sex_Xlinked_loci_to_remove <-  match(dartR_sex_Xlinked_loci, as.factor(locNames(gop_gl_1)))
gop_gl_1_no_sex_linked <- gop_gl_1[, -dartR_sex_Xlinked_loci_to_remove]
gop_gl_1_no_sex_linked@other$loc.metrics <- gop_gl_1_no_sex_linked@other$loc.metrics[-dartR_sex_Xlinked_loci_to_remove, ]

#check locnames and locmetrics match
all.equal(colnames(as.matrix(gop_gl_1_no_sex_linked)), gop_gl_1_no_sex_linked@other$loc.metrics$AlleleID)

#save(gop_gl_1_no_sex_linked, file="./data/gop_gl_1_no_sex_linked.rdata")


