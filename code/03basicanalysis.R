# set working directory
setwd("~/Dropbox/gp_2022")

library(dartR)
library(adegenet)
library(pegas)
library(vcfR)
library(dplyr)
library(hierfstat)
library(poppr)
library(tidyverse)
library(car)  
library(tqdm)

#load data
load("./data/gop_gl_1_no_sex_linked.rdata")

gop_gi <- gl2gi(gop_gl_1_no_sex_linked)

#minor allele freq
maf_gop <- minorAllele(gop_gi)
hist(maf_gop)
#per locus heterozygosity
locusstats_gop <- locus_table(gop_gi)
hist(locusstats_gop[,3])

#create structure/colony file
# define helper function/s
# make a function to convert two-row-per-individual structure format to one-row-per-individual structure format
#  Assumes: "data" is loaded as a matrix or data.frame and the first column is individual IDs
convert_str <- function(data) {
  
  out <- matrix(NA, nrow = (nrow(data) / 2), ncol = (2 * (ncol(data) - 1)))
  out[, seq(1, ncol(out), by = 2)] <- as.matrix(data[seq(1, nrow(data), by = 2), -1])
  out[, seq(2, ncol(out), by = 2)] <- as.matrix(data[seq(2, nrow(data), by = 2), -1])
  rownames(out) <- data[seq(1, nrow(data), by = 2), 1]
  
  out
  
}

str_data <- gl2structure(gop_gl_1_no_sex_linked,
                         exportMarkerNames = FALSE, outpath = "./data/",
                         outfile = "gop_gl_1_no_sex_linked.txt")

# read in two-row format data and convert to one-row format data
str_2row_format <- read.table("./data/gop_gl_1_no_sex_linked.txt", sep = "\t", header = FALSE)

#convert to one row per individual format
str_1row_format <- convert_str(str_2row_format)
#write.table(str_1row_format, "./data/gop_gl_1_no_sex_linked_with-9s.str")

#convert -9 values to 0s for colony
str_1row_with0s <- ifelse(str_1row_format == -9, 0, str_1row_format)
#write.table(str_1row_with0s, "./data/gop_gl_1_no_sex_linked_with0s.str")

#for Htest
str_1row_withNAs <- ifelse(str_1row_format == -9, "NA", str_1row_format)

#HL in R
# source helper functions
source("./code/ALF.R")
source("./code/GENHETv3.1.R")

data_HL <- str_1row_withNAs
data_HL_numeric <- apply(data_HL, 2, function(x) as.numeric(x))
rownames(data_HL_numeric) <- rownames(data_HL)
data_HL_numeric <- cbind(seq_len(nrow(data_HL_numeric)), data_HL_numeric)

locusname <- scan("./data/locusnames.txt",what="character",sep="\t")
#alfreqtest <- ALF(geno=data_HL,locname=locusname)
Htest <- GENHET(dat=data_HL_numeric,estimfreq="T",locname=locusname)
##commented out to avoid overwriting
#write.csv(Htest, "./outputs/Htest_output_2905loci.csv")
htest_results <- read.csv("./data/Htest_output_2905loci_cov.csv", header=TRUE)

pdf("outputs/PHt_2905loci.pdf", width = 10, height = 8)
par(mfrow = c(1, 1), mar = c(10, 5.1, 2, 1.1))

boxplot(PHt~pop,data=htest_results, main="", 
        xlab="", ylab="Individual heterozygosity", col="gray",
        outline = TRUE, boxwex = 0.5, frame=F, las=2, 
        cex.lab=1, cex.axis=0.9)

dev.off()

# pop genetic diversity metrics
gop_gl_1_no_sex_linked_noOVENS <- gl.drop.pop(gop_gl_1_no_sex_linked, pop.list="OVENS")
gop_gi_noOVENS <- gl2gi(gop_gl_1_no_sex_linked_noOVENS)

basicstats <- basic.stats(gop_gi_noOVENS,diploid=TRUE,digits=4)
N <- basicstats$n.ind.samp
Ho <-basicstats$Ho
Hs <-basicstats$Hs
Fis <-basicstats$Fis
overall <-basicstats$overall
#The following will give you mean Ho, Hs and Fis values per pop
popmeanHo <- colMeans(Ho, na.rm=TRUE)
popmeanHs <- colMeans(Hs, na.rm=TRUE)
popmeanFish <- colMeans(Fis, na.rm=TRUE)

#write.csv(basicstats$Hs, "./data/Hs_table.csv")
#write.csv(basicstats$Ho, "./data/Ho_table.csv")

# Calculate mean expected heterozygosity for each population
mean_heterozygosity <- colMeans(basicstats$Hs, na.rm = TRUE)
basicstats_data <- read.csv("./data/Hs_table.csv", row.names = 1)

# Perform Kruskal-Wallis test
kruskal_result <- kruskal.test(heterozygosity ~ population, data = locus_data)

# Print Kruskal-Wallis test results
print("Kruskal-Wallis test results:")
print(kruskal_result)

# Perform pairwise Wilcoxon rank sum tests with Bonferroni correction
pairwise_wilcox <- pairwise.wilcox.test(locus_data$heterozygosity, 
                                        locus_data$population, 
                                        p.adjust.method = "bonferroni")

# Print the first few rows of the pairwise Wilcoxon test results
print("\
Pairwise Wilcoxon rank sum test results (first few rows):")
print(head(pairwise_wilcox$p.value))

# Create a heatmap of p-values from pairwise comparisons
library(pheatmap)

# Convert p-values to a matrix
p_value_matrix <- pairwise_wilcox$p.value
p_value_matrix[is.na(p_value_matrix)] <- 1  # Replace NA with 1 (no significance)

# Create a color palette for the heatmap
color_palette <- colorRampPalette(c("red", "yellow", "green"))(100)

# Save the heatmap as a PNG file
png("pairwise_comparison_heatmap.png", width = 800, height = 600)
pheatmap(p_value_matrix, 
         color = color_palette,
         main = "Pairwise Wilcoxon Test P-values",
         fontsize = 8,
         angle_col = 45)
dev.off()

# Calculate effect size (epsilon squared) for Kruskal-Wallis test
n <- nrow(locus_data)
k <- length(unique(locus_data$population))
epsilon_squared <- (kruskal_result$statistic - k + 1) / (n - k)

print("\
Effect size (Epsilon squared) for Kruskal-Wallis test:")
print(epsilon_squared)

#Kruskal-Wallis chi-squared = 34.062, df = 20, p-value = 0.02571
#statistically significant differences in heterozygosity among populations:
#heatmap visually represents the p-values from the pairwise Wilcoxon rank sum tests, 
#with colors indicating statistical significance. 


# Calculate mean observed heterozygosity for each population
basicstats_data_ho <- read.csv("./data/Ho_table.csv", row.names = 1)

# Perform Kruskal-Wallis test
kruskal_result <- kruskal.test(heterozygosity ~ population, data = locus_data)

# Print Kruskal-Wallis test results
print("Kruskal-Wallis test results:")
print(kruskal_result)

# Perform pairwise Wilcoxon rank sum tests with Bonferroni correction
pairwise_wilcox <- pairwise.wilcox.test(locus_data$heterozygosity, 
                                        locus_data$population, 
                                        p.adjust.method = "bonferroni")

# Print the first few rows of the pairwise Wilcoxon test results
print("\
Pairwise Wilcoxon rank sum test results (first few rows):")
print(head(pairwise_wilcox$p.value))

# Create a heatmap of p-values from pairwise comparisons
library(pheatmap)

# Convert p-values to a matrix
p_value_matrix <- pairwise_wilcox$p.value
p_value_matrix[is.na(p_value_matrix)] <- 1  # Replace NA with 1 (no significance)

# Create a color palette for the heatmap
color_palette <- colorRampPalette(c("red", "yellow", "green"))(100)

# Save the heatmap as a PNG file
png("pairwise_comparison_heatmap_ho.png", width = 800, height = 600)
pheatmap(p_value_matrix, 
         color = color_palette,
         main = "Pairwise Wilcoxon Test P-values",
         fontsize = 8,
         angle_col = 45)
dev.off()

# Calculate effect size (epsilon squared) for Kruskal-Wallis test
n <- nrow(locus_data)
k <- length(unique(locus_data$population))
epsilon_squared <- (kruskal_result$statistic - k + 1) / (n - k)

print("\
Effect size (Epsilon squared) for Kruskal-Wallis test:")
print(epsilon_squared)

#Kruskal-Wallis chi-squared 0.0002305922 




#pop genetic diversity estimates removing highly heterozygous individuals

AR <- allelic.richness(gop_gi_noOVENS,min.n=NULL,diploid=TRUE)
ar <- AR$Ar
popmeanAR <- colMeans(ar, na.rm=TRUE)
#put it all in a table
table <- data.frame(Ho = popmeanHo, Hs = popmeanHs, Fis = popmeanFish, AR = popmeanAR )
#write table
#write.table(table, "./outputs/outputs_div_stats_FINAL_ALL.txt")

#calculate sample sizes
tapply(rep(1, length(gop_gi@pop)), as.character(gop_gi@pop), sum)                  

##remove het individuals
gop_gl_1_no_sex_linked_noOVENSnoHET <- gl.drop.ind(gop_gl_1_no_sex_linked_noOVENS, ind.list = c("GOP759",
                                                                                                "GOP757",
                                                                                                "GOP746",
                                                                                                "GOP735",
                                                                                                "GOP424",
                                                                                                "GOP730",
                                                                                                "GOP390",
                                                                                                "GOP661",
                                                                                                "GOP674",
                                                                                                "GOP742",
                                                                                                "GOP656",
                                                                                                "GOP003",
                                                                                                "GOP417",
                                                                                                "GOP756",
                                                                                                "GOP693",
                                                                                                "GOP740",
                                                                                                "GOP745",
                                                                                                "GOP403",
                                                                                                "GOP786",
                                                                                                "GOP691",
                                                                                                "GOP739",
                                                                                                "GOP676",
                                                                                                "GOP727",
                                                                                                "GOP736",
                                                                                                "GOP758",
                                                                                                "GOP670",
                                                                                                "GOP396",
                                                                                                "GOP749",
                                                                                                "GOP713",
                                                                                                "GOP655",
                                                                                                "GOP414",
                                                                                                "GOP525",
                                                                                                "GOP400"))


gop_gi_noOVENSnoHET <- gl2gi(gop_gl_1_no_sex_linked_noOVENSnoHET)

basicstats <- basic.stats(gop_gi_noOVENSnoHET,diploid=TRUE,digits=4)
N <- basicstats$n.ind.samp
Ho <-basicstats$Ho
Hs <-basicstats$Hs
Fis <-basicstats$Fis
overall <-basicstats$overall
#The following will give you mean Ho, Hs and Fis values per pop
popmeanHo <- colMeans(Ho, na.rm=TRUE)
popmeanHs <- colMeans(Hs, na.rm=TRUE)
popmeanFish <- colMeans(Fis, na.rm=TRUE)

AR <- allelic.richness(gop_gi_noOVENSnoHET,min.n=NULL,diploid=TRUE)
ar <- AR$Ar
popmeanAR <- colMeans(ar, na.rm=TRUE)
#put it all in a table
table <- data.frame(Ho = popmeanHo, Hs = popmeanHs, Fis = popmeanFish, AR = popmeanAR )
#write table
#write.table(table, "./outputs/outputs_div_stats_FINAL_ALL_noHets.txt")




#effective population size estimates using dartR ldne function
#need to download neestimator and give filepath

nes_gp <- gl.LDNe(gop_gl_1_no_sex_linked, 
                   critical = c(0.01,0.05), singleton.rm = T, 
                   neest.path = "./code/ZipFolder_32_bit_191125", 
                   mating = "random")


nes_gp_no_hets <- gl.LDNe(gop_gl_1_no_sex_linked_noOVENSnoHET, 
                  critical = c(0.01,0.05), singleton.rm = T, 
                  neest.path = "./code/ZipFolder_32_bit_191125", 
                  mating = "random")


library(hierfstat)

gop_wild_by_age  <- import2genind("./data/hierfstat_format_4agegroups_wild.DAT")
gop_stocked_by_age  <- import2genind("./data/hierfstat_format_4agegroups_stocked.DAT")

gop_stocked_by_age_gl <- gi2gl(gop_stocked_by_age)
nes_stocked_aged_gp <- gl.LDNe(gop_stocked_by_age_gl, 
                  critical = c(0.01,0.05), singleton.rm = T, 
                  neest.path = "./code/ZipFolder_32_bit_191125", 
                  mating = "random")

gop_wild_by_age_gl <- gi2gl(gop_wild_by_age)
nes_wild_aged_gp <- gl.LDNe(gop_wild_by_age_gl, 
                               critical = c(0.01,0.05), singleton.rm = T, 
                               neest.path = "./code/ZipFolder_32_bit_191125", 
                               mating = "random")


#Ne corrections for ldne estimates
##nb is the effective number of breeders in one reproductive cycle 
#(or time period)- reflects short-term Ne relevant
#to inbreeding
#ne is the effective population size per generation - more relevant to long-term evolutionary
#processes
AL = 23
alpha = 3

#example
raw_nb <- 970

#first correction for physical linkage
nb_adj1 = raw_nb / (0.098 + 0.219 * log(24))
nb_adj1 <- 5718
#second correction for mixed age sample - 2 trait version
nb2 = nb_adj1 / (1.103 - 0.245* log10(AL/alpha))
ne2 = nb2 / (0.485 + 0.758 * log10(AL/alpha))
ne2

