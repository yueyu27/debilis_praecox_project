#####################################
# Project: debilis_praecox_project
#
# Code 8: Inversion Comparison
#
# by: Yue Yu
#
#####################################

#####################################
#
#
# Part 1: Adaptive V.S. Non-Adaptive INV  -- INV LENGTH
#
#
#####################################

# Here we define the following:

cd /Users/yueyu/Desktop/INV_Chracter/Raw_INV_info

# - Adaptive inv
# The ones detected by Lostruct - as they are polymorphic in the pop, and are under balancing selection
# Saved as (adapt 12 INVs): lostruct_12_INV.txt


# - Non-adaptive inv
# The ones detected by comparing ref genome (DEB-PRA) (MINUS) overlap with Lostruct
# Saved as (all 81 INVs):  minimap_all_81_INV.txt
# Saved as (non adapt 74 INVs): minimap_nonadapt_74_INV.txt


#-------------------------
#
# GOAL 1: Check Length
#
#-------------------------

# Local R studio
# R version 4.5.2

setwd("/Users/yueyu/Desktop/INV_Chracter/Raw_INV_info")

adp <- read.delim("lostruct_11_INV.txt", header = T)
non_adp <- read.delim("minimap_nonadapt_75_INV.txt", header = T)

non_adp_PRAH1H2 <- read.delim("minimap_pra1and2_nonadapt_53_INV.txt", header = T)
non_adp_DEBH1H2 <- read.delim("minimap_deb1and2_all_and_nonadapt_17_INV.txt", header = T)


adp_length <- adp$LENGTH
non_adp_length <- c(non_adp$PRA_H2_LENGTH, non_adp_PRAH1H2$PRA_H2_LENGTH, non_adp_DEBH1H2$DEB_H2_LENGTH)


summary(adp_length)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max
#   8142965 15246439 24773886 23561832 29625890 43960821 

summary(non_adp_length)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#   507191   960280  2153884  5455781  6545103 96797362


# Run Comparative analyses:
# Check for normality: If data is normal, use Welch's t-test.
# Check for skew/outliers: If data is heavily skewed or has outliers, use Wilcoxon.


# Run a Welch t-test   <----- statistically significant (p = 9.774e-05 (***) < 0.01)
# - assumes normal distribution (MY DATA VIOLATES THIS)
# - allow unequal sample sizes
# - allow unequal variances 
# t.test(adp_length,non_adp_length, var.equal = FALSE)



# Also run a Wilcoxon rank-sum test  <----- statistically significant (p < 0.01)
# - allow non normal distribution
# - robust to outliers
# - 

wilcox.test(adp_length, non_adp_length)

# Wilcoxon rank sum test with continuity correction

data:  adp_length and non_adp_length
W = 1668, p-value = 9.099e-08. # DOUBLE CHECK THIS VALUE
alternative hypothesis: true location shift is not equal to 0



# Plotting results
library(ggplot2)

# prepare data
df <- data.frame(
  length = c(adp_length, non_adp_length),
  group = factor(c(rep("Adaptive", length(adp_length)),
                   rep("Non-Adaptive", length(non_adp_length))),
                 levels = c("Non-Adaptive", "Adaptive"))
)

# compute p-value (Welch t-test)
pval <- wilcox.test(adp_length, non_adp_length)$p.value
#  9.098948e-08

# convert to significance stars
stars <- ifelse(pval < 0.001, "***",
         ifelse(pval < 0.01, "**",
         ifelse(pval < 0.05, "*", "ns")))

# y position for annotation
y_max <- max(df$length)

ggplot(df, aes(x = group, y = length)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  
  # significance bar + star
  annotate("segment", x = 1, xend = 2, y = y_max*1.05, yend = y_max*1.05) +
  annotate("text", x = 1.5, y = y_max*1.1, label = stars, size = 6) +
  
  labs(x = NULL, y = "Inversion Length") +
  
  theme_classic() +
  theme(
    axis.text = element_text(size = 18, color = "black"),
    axis.title = element_text(size = 20, color = "black"),
    plot.title = element_text(hjust = 0.5)
  )



#####################################
#
#
# Part 2: Adaptive V.S. Non-Adaptive INV -- GENE
#
#
#####################################


# Goal: count the number of genes within all the windows


#--- Sundance:
cd /AsteroidScratch/Annotation_YueYu/INV_character_GENE_COUNT/raw

-rw-r--r-- 1 yue_yu  602 May 22 15:16 lostruct_11_INV.txt
-rw-r--r-- 1 yue_yu 6.4K Apr 23 09:56 minimap_all_81_INV.txt
-rw-r--r-- 1 yue_yu 5.6K May 22 15:16 minimap_nonadapt_75_INV.txt

# Split to process between species
# -- Adaptive INV ---
grep "deb" lostruct_11_INV.txt | cut -f3-5 > deb_adapt_inv.bed
grep "pra" lostruct_11_INV.txt | cut -f3-5 > pra_adapt_inv.bed

# -- Non-Adaptive INV ---
#sed '1d' minimap_nonadapt_75_INV.txt | cut -f5-7 > deb_NONadapt_inv.bed 
# DO NOT RUN THIS FOR DEB, ONLY USE PRA TO REPRESENT NUM OF GENES FOR THIS SET OF INVERSIONS!!

sed '1d' minimap_nonadapt_75_INV.txt | cut -f1-3 > pra_NONadapt_inv.bed


# -- Non-Adaptive INV (from H1H2 compare)---
sed '1d' minimap_deb1and2_all_and_nonadapt_17_INV.txt | cut -f5-7 > deb_NONadapt_inv.bed
# 17 = 17 INV (ALL COORD BASED ON DEB H2)

sed '1d' minimap_pra1and2_nonadapt_53_INV.txt | cut -f5-7 >> pra_NONadapt_inv.bed
# 75 + 53 = 128 INV (ALL COORD BASED ON PRA H2)




cd /AsteroidScratch/Annotation_YueYu/INV_character_GENE_COUNT

cp /AsteroidScratch/Annotation_YueYu/INV_character_GENE_COUNT/raw/*.bed .
# ======
# DEB H2 (Adaptive -- 466 genes)
# ======
GENE="/AsteroidScratch/for_Eric/annotation_result/DEB_H2/braker.gff3"
grep -P "\tgene\t" $GENE | awk 'BEGIN{OFS="\t"}{print $1, $4, $5}' > DEB_H2_gene_all.bed
sed -i 's/^Chr/DEB_&/' DEB_H2_gene_all.bed # ADD "DEB_" to all chromosomes 

bedtools intersect -a deb_adapt_inv.bed -b DEB_H2_gene_all.bed > deb_adapt_inv_overlap_gene.bed
# The output will have columns: Chrom, Start, End
# 466 genes

# Add gene count per INV
bedtools intersect -a deb_adapt_inv.bed -b DEB_H2_gene_all.bed -c > deb_adapt_inv_overlap_gene_withCOUNT.bed
awk 'BEGIN{OFS="\t"} {len=($3-$2)/1000000; print $0, $4/len}' deb_adapt_inv_overlap_gene_withCOUNT.bed > deb_adapt_inv_overlap_gene_withCOUNT_per1Mbp.bed



# ======
# DEB H2 (NON-Adaptive -- 609 genes)
# ======
bedtools intersect -a deb_NONadapt_inv.bed -b DEB_H2_gene_all.bed > deb_NONadapt_inv_overlap_gene.bed
# The output will have columns: Chrom, Start, End
# 609 genes

# Add gene count per INV
bedtools intersect -a deb_NONadapt_inv.bed -b DEB_H2_gene_all.bed -c > deb_NONadapt_inv_overlap_gene_withCOUNT.bed


#awk 'BEGIN{OFS="\t"}{len=($3-$2)/1000000; print $0, $4/len}' \
#deb_NONadapt_inv_overlap_gene_withCOUNT.bed > deb_NONadapt_inv_overlap_gene_withCOUNT_per1Mbp.bed






# ======
# PRA H2 (Adaptive -- 3305 genes)
# ======
GENE="/AsteroidScratch/for_Eric/annotation_result/PRA_H2/praecox_H2_renamed_F.gff3" 
grep -P "\tgene\t" $GENE | awk 'BEGIN{OFS="\t"}{print $1, $4, $5}' > PRA_H2_gene_all.bed
sed -i 's/^Chr/PRA_chr/' PRA_H2_gene_all.bed # Change "PRA_chr" to all chromosomes 
awk 'BEGIN{OFS="\t"}{
  if ($1 ~ /^Chr[0-9]+$/) {
    num = substr($1,4)                  # extract number after "Chr"
    $1 = sprintf("PRA_chr%02d", num)    # zero-pad to 2 digits
  }
  print
}' PRA_H2_gene_all.bed > PRA_H2_gene_all_fixed.bed

cut -f1 PRA_H2_gene_all_fixed.bed | sort | uniq # Check finally correct wow


bedtools intersect -a pra_adapt_inv.bed -b PRA_H2_gene_all_fixed.bed > pra_adapt_inv_overlap_gene.bed
# The output will have columns: Chrom, Start, End
# 3305 genes

# Add gene count per INV
bedtools intersect -a pra_adapt_inv.bed -b PRA_H2_gene_all_fixed.bed -c > pra_adapt_inv_overlap_gene_withCOUNT.bed

awk 'BEGIN{OFS="\t"}{len=($3-$2)/1000000; print $0, $4/len}' pra_adapt_inv_overlap_gene_withCOUNT.bed > pra_adapt_inv_overlap_gene_withCOUNT_per1Mbp.bed


# ======
# PRA H2 (NON-Adaptive -- 8221 genes)
# ======
bedtools intersect -a pra_NONadapt_inv.bed -b PRA_H2_gene_all_fixed.bed > pra_NONadapt_inv_overlap_gene.bed
# The output will have columns: Chrom, Start, End
# 8221 genes

# Add gene count per INV
bedtools intersect -a pra_NONadapt_inv.bed -b PRA_H2_gene_all_fixed.bed -c > pra_NONadapt_inv_overlap_gene_withCOUNT.bed

#awk 'BEGIN{OFS="\t"}{len=($3-$2)/1000000; print $0, $4/len}' \
#pra_NON_adapt_inv_overlap_gene_withCOUNT.bed > pra_NON_adapt_inv_overlap_gene_withCOUNT_per1Mbp.bed


# ============= R plotting ============

# R version 4.3.1

# -- Adaptive (merge) 
adp_deb <- read.delim("deb_adapt_inv_overlap_gene_withCOUNT.bed", header = F)
adp_pra <- read.delim("pra_adapt_inv_overlap_gene_withCOUNT.bed", header = F)

adp <- rbind(adp_deb,adp_pra)
colnames(adp) <- c("CHR","START","END","GENE_COUNT")
# 11 adaptive

# -- Non Adaptive (keep PRA)
# Only use PRA for non adaptive, as they are corresponding regions
non_adp_deb <- read.delim("deb_NONadapt_inv_overlap_gene_withCOUNT.bed", header = F)
non_adp_pra <- read.delim("pra_NONadapt_inv_overlap_gene_withCOUNT.bed", header = F)

non_adp <- rbind(non_adp_deb,non_adp_pra)
colnames(non_adp) <- c("CHR","START","END","GENE_COUNT")
# 145 non adaptiove



########################
#
#     Gene count 
#
########################

# -- Extract comparison 
adp_genecount <- adp$GENE_COUNT
non_adp_genecount <- non_adp$GENE_COUNT

summary(adp_genecount)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max
#     98.0   201.5   307.0   342.8   505.5   606.0

summary(non_adp_genecount)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max
#    0.00   12.00   27.00   60.96   71.00  651.00




wilcox.test(adp_genecount, non_adp_genecount) ############    Significant p == 6.416e-07

# Output
	Wilcoxon rank sum test with continuity correction

data:  adp_genecount and non_adp_genecount
W = 1693.5, p-value = 6.416e-07
alternative hypothesis: true location shift is not equal to 0



# Plotting results
library(ggplot2)

# prepare data
df <- data.frame(
  count = c(adp_genecount, non_adp_genecount),
  group = factor(c(rep("Adaptive", length(adp_genecount)),
                   rep("Non-Adaptive", length(non_adp_genecount))),
                 levels = c("Non-Adaptive", "Adaptive"))
)

# compute p-value (Welch t-test)
pval <- wilcox.test(adp_genecount, non_adp_genecount)$p.value

# convert to significance stars
stars <- ifelse(pval < 0.001, "***",
         ifelse(pval < 0.01, "**",
         ifelse(pval < 0.05, "*", "ns")))

# y position for annotation
y_max <- max(df$count)

p <- ggplot(df, aes(x = group, y = count)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  
  # significance bar + star
  annotate("segment", x = 1, xend = 2, y = y_max*1.05, yend = y_max*1.05) +
  annotate("text", x = 1.5, y = y_max*1.1, label = stars, size = 6) +
  
  labs(x = NULL, y = "Gene Count") +
  
  theme_classic() +
  theme(
    axis.text = element_text(size = 18, color = "black"),
    axis.title = element_text(size = 20, color = "black"),
    plot.title = element_text(hjust = 0.5)
  )

png(filename = "gene_count.png",  width = 1000, height = 800, res = 200)
print(p)
dev.off()


# publication PDF (Completed)
pdf(file = "gene_count.pdf", height = 4.5, width = 6.5)
print(p)
dev.off()




#####################################
#
#
# Part 1+2: CORR length <--> gene count
#
#
#####################################

# Move INV gene count to local laptop

# ===== GENE COUNT =========

cd /Users/yueyu/Desktop/INV_Chracter/Gene_count
scp yue_yu@sundance.zoology.ubc.ca:"/AsteroidScratch/Annotation_YueYu/INV_character_GENE_COUNT/*gene_withCOUNT.bed" .

ls -thor
total 32
-rw-r--r--@ 1 yueyu   541B Jun  4 17:50 deb_NONadapt_inv_overlap_gene_withCOUNT.bed
-rw-r--r--@ 1 yueyu    65B Jun  4 17:50 deb_adapt_inv_overlap_gene_withCOUNT.bed
-rw-r--r--@ 1 yueyu   4.0K Jun  4 17:50 pra_NONadapt_inv_overlap_gene_withCOUNT.bed
-rw-r--r--@ 1 yueyu   293B Jun  4 17:50 pra_adapt_inv_overlap_gene_withCOUNT.bed

# ===== INV LENGTH =========
cd /Users/yueyu/Desktop/INV_Chracter/Raw_INV_info

lostruct_11_INV.txt
minimap_nonadapt_75_INV.txt
minimap_pra1and2_nonadapt_53_INV.txt
minimap_deb1and2_all_and_nonadapt_17_INV.txt


#-------------------------
#
# STEP 1: Load both into R
#
#-------------------------

# Local R studio
# R version 4.5.2

setwd("/Users/yueyu/Desktop/INV_Chracter/Raw_INV_info")

adp <- read.delim("lostruct_11_INV.txt", header = T)
non_adp <- read.delim("minimap_nonadapt_75_INV.txt", header = T)

non_adp_PRAH1H2 <- read.delim("minimap_pra1and2_nonadapt_53_INV.txt", header = T)
non_adp_DEBH1H2 <- read.delim("minimap_deb1and2_all_and_nonadapt_17_INV.txt", header = T)

yes_length <- adp[,c(3,4,5,6)] # USE THIS
no_length <- non_adp[,c(1,2,3,4)]
no_length <- rbind(no_length,non_adp_PRAH1H2[,c(5,6,7,8)])
colnames(no_length) <- c("CHR","START","END","LENGTH")
colnames(non_adp_DEBH1H2)[5:8] <- c("CHR","START","END","LENGTH")
no_length <- rbind(no_length,non_adp_DEBH1H2[,c(5,6,7,8)])

rm(adp,non_adp,non_adp_PRAH1H2,non_adp_DEBH1H2)



setwd("/Users/yueyu/Desktop/INV_Chracter/Gene_count")

adp_deb <- read.delim("deb_adapt_inv_overlap_gene_withCOUNT.bed", header = F)
adp_pra <- read.delim("pra_adapt_inv_overlap_gene_withCOUNT.bed", header = F)

yes_gene <- rbind(adp_deb,adp_pra)
colnames(yes_gene) <- c("CHR","START","END","GENE_COUNT")
# 11 adaptive

# -- Non Adaptive (keep PRA)
# Only use PRA for non adaptive, as they are corresponding regions
non_adp_deb <- read.delim("deb_NONadapt_inv_overlap_gene_withCOUNT.bed", header = F)
non_adp_pra <- read.delim("pra_NONadapt_inv_overlap_gene_withCOUNT.bed", header = F)

no_gene <- rbind(non_adp_deb,non_adp_pra)
colnames(no_gene) <- c("CHR","START","END","GENE_COUNT")
# 145 non adaptiove

rm(adp_deb,adp_pra,non_adp_deb,non_adp_pra)



#-------------------------
#
# STEP 2: Merge into one dataframe
#
#-------------------------
YES <- merge(yes_length,yes_gene, by = c("CHR","START","END"))
NO <- merge(no_length,no_gene, by = c("CHR","START","END"))


YES_AND_NO <- rbind(YES,NO)

cor.test(YES_AND_NO$LENGTH, YES_AND_NO$GENE_COUNT)





#####################################
#
#
# Part 3 - Use gene density <- infer -> recombination rate in DEB H2 & PRA H2
#
#
#####################################

# ======= Step 1: Prepare window.bed file (1Mbp) + gene density

cd /AsteroidScratch/Syri_YueYu_tmp/Recomb_genedensity_predict

# ======
#  DEB
# ======
# Previously prepared for Circo plot, move here
cp /AsteroidScratch/Annotation_YueYu/COUNT_FOR_CIRCO/gene_density_1Mbp_DEB_H2.bed .

# ======
#  PRA
# ======
cp /AsteroidScratch/Annotation_YueYu/COUNT_FOR_CIRCO/gene_density_1Mbp_PRA_H2.bed .


ls -thor
total 204K
-rw-rw-r-- 1 yue_yu 102K May  5 18:55 gene_density_1Mbp_DEB_H2.bed
-rw-rw-r-- 1 yue_yu  97K May  5 18:55 gene_density_1Mbp_PRA_H2.bed




# ======= Step 3: Use HA412 correlation to predict rrate for new ref 
# The formula used for prediction for DEB and PRA:
# rrate= −0.1791 + 0.0421×gene_count

R

deb <- read.delim("gene_density_1Mbp_DEB_H2.bed", header = FALSE)
colnames(deb) <- c("CHR","START","END","gene_count")

pra <- read.delim("gene_density_1Mbp_PRA_H2.bed", header = FALSE)
colnames(pra) <- c("CHR","START","END","gene_count")

deb$rrate <- (-0.1791) + 0.0421 * deb$gene_count
pra$rrate <- (-0.1791) + 0.0421 * pra$gene_count

# Save optional







# ======= Step 4: Plot predicted rrcombination rate for each chromosome

library(dplyr)
library(ggplot2)

# midpoint in Mb
deb_plot <- deb %>%
  mutate(mid_mb = ((START + END) / 2) / 1e6)

# loop through chromosomes
for(chr in unique(deb_plot$CHR)) {

  p <- deb_plot %>%
    filter(CHR == chr) %>%
    ggplot(aes(x = mid_mb, y = rrate)) +
    
    geom_line(linewidth = 0.7) +
    
    theme_classic(base_size = 14) +
    
    labs(
      title = chr,
      x = "Position (Mb)",
      y = "Recombination rate"
    ) +
    
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text = element_text(color = "black")
    )

	png(filename = paste0("rrate_deb/",chr,"_rrate.png"),  width = 1600, height = 800, res = 200)
  	print(p)
	dev.off()
}



# midpoint in Mb
pra_plot <- pra %>%
  mutate(mid_mb = ((START + END) / 2) / 1e6)

# loop through chromosomes
for(chr in unique(pra_plot$CHR)) {

  p <- pra_plot %>%
    filter(CHR == chr) %>%
    ggplot(aes(x = mid_mb, y = rrate)) +
    
    geom_line(linewidth = 0.7) +
    
    theme_classic(base_size = 14) +
    
    labs(
      title = chr,
      x = "Position (Mb)",
      y = "Recombination rate"
    ) +
    
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text = element_text(color = "black")
    )

	png(filename = paste0("rrate_pra/",chr,"_rrate.png"),  width = 1600, height = 800, res = 200)
  	print(p)
	dev.off()
}







# ======= Step 5: Inversion overlapping


cd /AsteroidScratch/Annotation_YueYu/INV_character_GENE_COUNT/raw

ls -thor
# lostruct_11_INV.txt
# minimap_nonadapt_75_INV.txt



# inside R again, continue from code above in step 4:

floor_to <- function(x, unit=1000000) {
  (x %/% unit) * unit
}


ceiling_to <- function(x, unit = 1000000) {
  ((x + 999999) %/% unit) * unit
}



# ====== Make it into a loop ====

# ---- DEB + ADAPTIVE
inv <- read.delim("/AsteroidScratch/Annotation_YueYu/INV_character_GENE_COUNT/raw/deb_adapt_inv.bed", header = FALSE)
colnames(inv) <- c("CHR","START","END")

inv$rrate <- "NA"

for (i in 1:nrow(inv)) {

	chr <- inv$CHR[i]
	start <- floor_to(inv$START[i])
	end <- ceiling_to(inv$END[i])

	rate_within_INV <- deb_plot %>% filter(CHR == chr & START >= start & END <= end )	
	inv$rrate[i] <- round(mean(rate_within_INV$rrate), digits = 5)

	print(i)
}

deb_adap_inv <- inv



# ---- PRA + ADAPTIVE
inv <- read.delim("/AsteroidScratch/Annotation_YueYu/INV_character_GENE_COUNT/raw/pra_adapt_inv.bed", header = FALSE)
colnames(inv) <- c("CHR","START","END")

inv$rrate <- "NA"

for (i in 1:nrow(inv)) {

	chr <- inv$CHR[i]
	start <- floor_to(inv$START[i])
	end <- ceiling_to(inv$END[i])

  rate_within_INV <- pra_plot %>% filter(CHR == chr & START >= start & END <= end ) 
  inv$rrate[i] <- round(mean(rate_within_INV$rrate), digits = 5)

	print(i)
}

pra_adap_inv <- inv



# ---- DEB + NON !!! ADAPTIVE
inv <- read.delim("/AsteroidScratch/Annotation_YueYu/INV_character_GENE_COUNT/raw/deb_NONadapt_inv.bed", header = FALSE)
colnames(inv) <- c("CHR","START","END")

inv$rrate <- "NA"

for (i in 1:nrow(inv)) {

	chr <- inv$CHR[i]
	start <- floor_to(inv$START[i])
	end <- ceiling_to(inv$END[i])

  rate_within_INV <- deb_plot %>% filter(CHR == chr & START >= start & END <= end ) 
  inv$rrate[i] <- round(mean(rate_within_INV$rrate), digits = 5)

	print(i)
}

deb_NONadap_inv <- inv


# ---- PRA + NON !!! ADAPTIVE (ONLY BETWEEN H1 H2)
inv <- read.delim("/AsteroidScratch/Syri_YueYu_tmp/Divergence_pi_2026May04/INV_POS/pra_NONadapt_53_inv_chr_rename.bed", header = FALSE)
colnames(inv) <- c("CHR","START","END")
inv$CHR <- sub("Chr", "PRA_chr", inv$CHR) # Change chr name to match filter

dim(inv)
# 53 3
# should be 53 rows instead of 127, exclude signal from between species comparison INVs

inv$rrate <- "NA"

for (i in 1:nrow(inv)) {

	chr <- inv$CHR[i]
	start <- floor_to(inv$START[i])
	end <- ceiling_to(inv$END[i])

  rate_within_INV <- pra_plot %>% filter(CHR == chr & START >= start & END <= end ) 
  inv$rrate[i] <- round(mean(rate_within_INV$rrate), digits = 5)

	print(i)
}

pra_NONadap_inv <- inv





# ---- POOL ALL 
NONadapt <- rbind(deb_NONadap_inv,pra_NONadap_inv) # 70 rows = 17 + 53 (ONLY WITHIN SPECIES)
adapt <- rbind(deb_adap_inv, pra_adap_inv)        # 11 rows

NONadapt_rrate <- as.numeric(NONadapt$rrate)
adapt_rrate <- as.numeric(adapt$rrate)

summary(NONadapt_rrate)
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# -0.1580  0.1265  0.3005  0.3811  0.4524  1.7996
summary(adapt_rrate)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.02672 0.28715 0.46697 0.50923 0.59340 1.46280

wilcox.test(adapt_rrate, NONadapt_rrate)

# Wilcoxon rank sum test with continuity correction

data:  adapt_rrate and NONadapt_rrate
W = 500, p-value = 0.1144
alternative hypothesis: true location shift is not equal to 0


# Plotting results
library(ggplot2)

# prepare data
df <- data.frame(
  rrate = c(adapt_rrate, NONadapt_rrate),
  group = factor(c(rep("Adaptive", length(adapt_rrate)),
                   rep("Non-Adaptive", length(NONadapt_rrate))),
                 levels = c("Non-Adaptive", "Adaptive"))
)

# compute p-value (Welch t-test)
pval <- wilcox.test(adapt_rrate, NONadapt_rrate)$p.value
#  0.105868

# convert to significance stars
stars <- ifelse(pval < 0.001, "***",
         ifelse(pval < 0.01, "**",
         ifelse(pval < 0.05, "*", "ns")))

# y position for annotation
y_max <- max(df$rrate)

p <- ggplot(df, aes(x = group, y = rrate)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  
  # significance bar + star
  annotate("segment", x = 1, xend = 2, y = y_max*1.05, yend = y_max*1.05) +
  annotate("text", x = 1.5, y = y_max*1.1, label = stars, size = 6) +
  
  labs(x = NULL, y = "Recombination Rate") +
  
  theme_classic() +
  theme(
    axis.text = element_text(size = 18, color = "black"),
    axis.title = element_text(size = 20, color = "black"),
    plot.title = element_text(hjust = 0.5)
  )

png(filename = "rrate_within_INV_2026June06.png",  width = 1000, height = 800, res = 200)
print(p)
dev.off()



# publication PDF (Completed)
pdf(file = "rrate_within_INV_2026June06.pdf", height = 4.5, width = 6.5)
print(p)
dev.off()





#####################################
#
#
# Part 4: Divergence -- Pi
#
#
#####################################

# ---------------
#
# Step 1: Find all INV positions based on H2
#
# ---------------

#--- Sundance:
cd /AsteroidScratch/Annotation_YueYu/INV_character_GENE_COUNT/raw

# Split to process between species
# -- Adaptive INV ---
grep "deb" lostruct_11_INV.txt | cut -f3-5 > deb_adapt_inv.bed
grep "pra" lostruct_11_INV.txt | cut -f3-5 > pra_adapt_inv.bed


# -- Non-Adaptive INV (from H1H2 compare)---
sed '1d' minimap_deb1and2_all_and_nonadapt_17_INV.txt | cut -f5-7 | sed 's/DEB_Chr/Chr/g' >> /AsteroidScratch/Syri_YueYu_tmp/Divergence_pi_2026May04/INV_POS/deb_NONadapt_17_inv_chr_rename.bed
# 17 INV (ALL COORD BASED ON DEB H2)

sed '1d' minimap_pra1and2_nonadapt_53_INV.txt | cut -f5-7 | sed 's/PRA_chr/Chr/g' >> /AsteroidScratch/Syri_YueYu_tmp/Divergence_pi_2026May04/INV_POS/pra_NONadapt_53_inv_chr_rename.bed
# 53 INV (ALL COORD BASED ON PRA H2)




#--- Move to new DIR:
cd /AsteroidScratch/Syri_YueYu_tmp/Divergence_pi_2026May04/INV_POS

# POS in VCF : "Chr01 ... Chr 17"
# POS in INV POS : "DEB_Chr" OR "PRA_chr" --> Change to match above

ls -thor
-rw-rw-r-- 1 yue_yu   75 May  5 10:20 deb_adapt_inv_chr_rename.bed
-rw-rw-r-- 1 yue_yu  222 May  5 10:20 pra_adapt_inv_chr_rename.bed
-rw-rw-r-- 1 yue_yu  425 May  5 15:08 deb_NONadapt_17_inv_chr_rename.bed
-rw-rw-r-- 1 yue_yu 1.3K May  5 15:08 pra_NONadapt_53_inv_chr_rename.bed






# ---------------
#
# Step 2: minimap -> bam -> syri (VCF)
#
# ---------------

# Sundance 

#=========
#  DEB
#=========

cd /AsteroidScratch/Syri_YueYu_tmp/Divergence_pi_2026May04/DEB

mkdir DEB_H1
cd DEB_H1
sed 's/DEB_Chr/Chr/g' /AsteroidScratch/Annotation_YueYu/LIFTOFF/data/debilis1.fasta > debilis1_renamedCHR.fasta
awk '
/^>/ {
    id=$1
    gsub(/^>/,"",id)

    # only keep main chromosomes
    if (id ~ /^Chr[0-9]+$/) {
        keep=1
        out=id".fasta"
    } else {
        keep=0
    }
}

keep { print > out }
' debilis1_renamedCHR.fasta

# DEB 2 already done, stored at: /AsteroidScratch/Syri_YueYu_tmp/Recomb_calc_per_CHROM/DEB_H2


# --- DEB H2 (ref) to DEB H1 (qry) ---

cd /AsteroidScratch/Syri_YueYu_tmp/Divergence_pi_2026May04/DEB

for chr in Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11; do
  nohup bash -c "
    /AsteroidScratch/Syri_YueYu_tmp/minimap2-2.28/minimap2 -ax asm5 --eqx -t 6 /AsteroidScratch/Syri_YueYu_tmp/Recomb_calc_per_CHROM/DEB_H2/${chr}.fasta DEB_H1/${chr}.fasta | \
    samtools sort -o result_deb/DEB_hap2_hap1_${chr}.bam
  " > ${chr}.log 2>&1 &
done

# Completed within 30min




for chr in Chr12 Chr13 Chr14 Chr15 Chr16 Chr17; do
  nohup bash -c "
    /AsteroidScratch/Syri_YueYu_tmp/minimap2-2.28/minimap2 -ax asm5 --eqx -t 6 /AsteroidScratch/Syri_YueYu_tmp/Recomb_calc_per_CHROM/DEB_H2/${chr}.fasta DEB_H1/${chr}.fasta | \
    samtools sort -o result_deb/DEB_hap2_hap1_${chr}.bam
  " > ${chr}.log 2>&1 &
done



# BAM -> SYRI VCF

for chr in Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 Chr12 Chr13 Chr14 Chr15 Chr16 Chr17; do
  nohup bash -c "
  source /DATA/home/yue_yu/miniconda3/etc/profile.d/conda.sh
  conda activate syri_plotsr 
  syri -c result_deb/DEB_hap2_hap1_${chr}.bam -r /AsteroidScratch/Syri_YueYu_tmp/Recomb_calc_per_CHROM/DEB_H2/${chr}.fasta -q DEB_H1/${chr}.fasta -F B --dir syri_deb/ --prefix DEB_hap2_hap1_${chr}_
  " > syri_deb_${chr}.log 2>&1 &
done

# SYRI VCF -> Merge all  (NNED TO DO!!)
for f in DEB_hap2_hap1_Chr*_syri.vcf; do
    bgzip -c $f > ${f}.gz
    tabix -p vcf ${f}.gz
done

bcftools concat -a -Oz -o DEB_hap2_hap1_merged.vcf.gz DEB_hap2_hap1_Chr*_syri.vcf.gz
tabix -p vcf DEB_hap2_hap1_merged.vcf.gz
bcftools index DEB_hap2_hap1_merged.vcf.gz








#=========
#  PRA
#=========

cd /AsteroidScratch/Syri_YueYu_tmp/Divergence_pi_2026May04/PRA

mkdir PRA_H1
cd PRA_H1
sed 's/PRA_chr/Chr/g' /AsteroidScratch/Annotation_YueYu/LIFTOFF/data/praecox1.fasta > praecox1_renamedCHR.fasta


awk '
/^>/ {
    id=$1
    gsub(/^>/,"",id)

    # only keep main chromosomes
    if (id ~ /^Chr[0-9]+$/) {
        keep=1
        out=id".fasta"
    } else {
        keep=0
    }
}

keep { print > out }
' praecox1_renamedCHR.fasta

# PRA 2 already done, stored at: /AsteroidScratch/Syri_YueYu_tmp/Recomb_calc_per_CHROM/PRA_H2


# --- PRA H2 (ref) to PRA H1 (qry) ---


cd /AsteroidScratch/Syri_YueYu_tmp/Divergence_pi_2026May04/PRA

for chr in Chr01 Chr02 Chr03 Chr04 Chr05 Chr06; do
  nohup bash -c "
    /AsteroidScratch/Syri_YueYu_tmp/minimap2-2.28/minimap2 -ax asm5 --eqx -t 6 /AsteroidScratch/Syri_YueYu_tmp/Recomb_calc_per_CHROM/PRA_H2/${chr}.fasta PRA_H1/${chr}.fasta | \
    samtools sort -o result_pra/PRA_hap2_hap1_${chr}.bam
  " > pra_${chr}.log 2>&1 &
done
# COMPLETED





for chr in Chr07 Chr08 Chr09 Chr10 Chr11 Chr12 Chr13 Chr14 Chr15 Chr16 Chr17; do
  nohup bash -c "
    /AsteroidScratch/Syri_YueYu_tmp/minimap2-2.28/minimap2 -ax asm5 --eqx -t 6 /AsteroidScratch/Syri_YueYu_tmp/Recomb_calc_per_CHROM/PRA_H2/${chr}.fasta PRA_H1/${chr}.fasta | \
    samtools sort -o result_pra/PRA_hap2_hap1_${chr}.bam
  " > pra_${chr}.log 2>&1 &
done

# RUNNING 19:38 -- COMPLETED 20:33





# BAM -> SYRI VCF

for chr in Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 Chr12 Chr13 Chr14 Chr15 Chr16 Chr17; do
  nohup bash -c "
  source /DATA/home/yue_yu/miniconda3/etc/profile.d/conda.sh
  conda activate syri_plotsr 
  syri -c result_pra/PRA_hap2_hap1_${chr}.bam -r /AsteroidScratch/Syri_YueYu_tmp/Recomb_calc_per_CHROM/PRA_H2/${chr}.fasta -q PRA_H1/${chr}.fasta -F B --dir syri_pra/ --prefix PRA_hap2_hap1_${chr}_
  " > syri_pra_${chr}.log 2>&1 &
done


# SYRI VCF -> Merge all 
for f in PRA_hap2_hap1_Chr*_syri.vcf; do
    bgzip -c $f > ${f}.gz
    tabix -p vcf ${f}.gz
done

bcftools concat -a -Oz -o PRA_hap2_hap1_merged.vcf.gz PRA_hap2_hap1_Chr*_syri.vcf.gz
tabix -p vcf PRA_hap2_hap1_merged.vcf.gz
bcftools index PRA_hap2_hap1_merged.vcf.gz







# ---------------
#
# Step 3: Extract poistion from VCF
#
# ---------------

# POS in VCF : "Chr01 ... Chr 17"
# POS in INV POS : "DEB_Chr" OR "PRA_chr" --> Changed to match above (DONE)

cd /AsteroidScratch/Syri_YueYu_tmp/Divergence_pi_2026May04/INV_POS

sed 's/DEB_//' deb_adapt_inv.bed > deb_adapt_inv_chr_rename.bed
sed 's/PRA_c/C/' pra_adapt_inv.bed > pra_adapt_inv_chr_rename.bed

sed 's/DEB_//' deb_NONadapt_inv.bed > deb_NONadapt_17_inv_chr_rename.bed
sed 's/PRA_c/C/' pra_NONadapt_inv.bed > pra_NONadapt_53_inv_chr_rename.bed


-rw-rw-r-- 1 yue_yu   75 May  5 10:20 deb_adapt_inv_chr_rename.bed
-rw-rw-r-- 1 yue_yu  222 May  5 10:20 pra_adapt_inv_chr_rename.bed
-rw-rw-r-- 1 yue_yu  425 May  5 15:08 deb_NONadapt_17_inv_chr_rename.bed
-rw-rw-r-- 1 yue_yu 1.3K May  5 15:08 pra_NONadapt_53_inv_chr_rename.bed


#=========
#  DEB
#=========

cd /AsteroidScratch/Syri_YueYu_tmp/Divergence_pi_2026May04/DEB
mkdir -p sub_VCF_adapt sub_VCF_NON_adapt


#  --- adapt inversion subset 
INV="/AsteroidScratch/Syri_YueYu_tmp/Divergence_pi_2026May04/INV_POS/deb_adapt_inv_chr_rename.bed"
INPUTVCF="DEB_hap2_hap1_merged.vcf.gz"

cd /AsteroidScratch/Syri_YueYu_tmp/Divergence_pi_2026May04/DEB/syri_deb

while read chr start end; do
    out="/AsteroidScratch/Syri_YueYu_tmp/Divergence_pi_2026May04/DEB/sub_VCF_adapt/${chr}_${start}_${end}.vcf.gz"
    bcftools view -r ${chr}:${start}-${end} $INPUTVCF -Oz -o $out
    tabix -p vcf $out
done < $INV


#  --- NON adapt inversion subset 
INV="/AsteroidScratch/Syri_YueYu_tmp/Divergence_pi_2026May04/INV_POS/deb_NONadapt_17_inv_chr_rename.bed"
INPUTVCF="DEB_hap2_hap1_merged.vcf.gz"

cd /AsteroidScratch/Syri_YueYu_tmp/Divergence_pi_2026May04/DEB/syri_deb

while read chr start end; do
    out="/AsteroidScratch/Syri_YueYu_tmp/Divergence_pi_2026May04/DEB/sub_VCF_NON_adapt/${chr}_${start}_${end}.vcf.gz"
    bcftools view -r ${chr}:${start}-${end} $INPUTVCF -Oz -o $out
    tabix -p vcf $out
done < $INV






#=========
#  PRA
#=========
cd /AsteroidScratch/Syri_YueYu_tmp/Divergence_pi_2026May04/PRA
mkdir -p sub_VCF_adapt sub_VCF_NON_adapt

#  --- adapt inversion subset 
INV="/AsteroidScratch/Syri_YueYu_tmp/Divergence_pi_2026May04/INV_POS/pra_adapt_inv_chr_rename.bed"
INPUTVCF="PRA_hap2_hap1_merged.vcf.gz"

cd /AsteroidScratch/Syri_YueYu_tmp/Divergence_pi_2026May04/PRA/syri_pra

while read chr start end; do
    out="/AsteroidScratch/Syri_YueYu_tmp/Divergence_pi_2026May04/PRA/sub_VCF_adapt/${chr}_${start}_${end}.vcf.gz"
    bcftools view -r ${chr}:${start}-${end} $INPUTVCF -Oz -o $out
    tabix -p vcf $out
done < $INV


#  --- NON adapt inversion subset 
INV="/AsteroidScratch/Syri_YueYu_tmp/Divergence_pi_2026May04/INV_POS/pra_NONadapt_53_inv_chr_rename.bed"
INPUTVCF="PRA_hap2_hap1_merged.vcf.gz"

cd /AsteroidScratch/Syri_YueYu_tmp/Divergence_pi_2026May04/PRA/syri_pra

while read chr start end; do
    out="/AsteroidScratch/Syri_YueYu_tmp/Divergence_pi_2026May04/PRA/sub_VCF_NON_adapt/${chr}_${start}_${end}.vcf.gz"
    bcftools view -r ${chr}:${start}-${end} $INPUTVCF -Oz -o $out
    tabix -p vcf $out
done < $INV





# ---------------
#
# Step 4: custom code count how many different nucleotide (SNP) / all nucleotide (SYN) in that INV
#
# ---------------

# CURRENT FORMULA (2026 May 05)

π SNP = number of SNPs/ aligned length (<SYN> + <INV>)

# INV opposite direction aligned region
# SYN: same dir aligned region


#=========
#  DEB
#=========

#  --- adapt inversion subset
cd /AsteroidScratch/Syri_YueYu_tmp/Divergence_pi_2026May04/DEB/sub_VCF_adapt

for f in *.vcf.gz; do

    base=$(basename $f .vcf.gz)

    zcat $f | awk -v name="$base" '
    BEGIN {
        snp=0;
        aligned_bp=0;
    }

    !/^#/ {
        ref=$4;
        alt=$5;

        # extract END
        match($8, /END=([0-9]+)/, arr);
        end=arr[1];
        len = end - $2 + 1;

        # --- denominator: aligned regions ---
        if (alt == "<SYN>" || alt == "<INV>") {
            aligned_bp += len;
        }

        # --- numerator: SNPs only ---
        if (alt !~ /^</ && length(ref)==1 && length(alt)==1) {
            snp++;
        }
    }

    END {
        if (aligned_bp > 0) {
            pi = snp / aligned_bp;
        } else {
            pi = 0;
        }

        printf "%s\tSNPs:%d\tAligned_bp:%d\tPi_SNP:%0.6e\n", name, snp, aligned_bp, pi;
    }'

done


#  --- NON adapt inversion subset
cd /AsteroidScratch/Syri_YueYu_tmp/Divergence_pi_2026May04/DEB/sub_VCF_NON_adapt
# Same code as above


#=========
#  PRA
#=========

#  --- adapt inversion subset
cd /AsteroidScratch/Syri_YueYu_tmp/Divergence_pi_2026May04/PRA/sub_VCF_adapt
# Same code as above


#  --- NON adapt inversion subset
cd /AsteroidScratch/Syri_YueYu_tmp/Divergence_pi_2026May04/PRA/sub_VCF_NON_adapt
# Same code as above

# Saved in Excel sheet on local laptop ready for plotting


# ---------------
#
# Step 5: Extract -> T test
#
# ---------------

# If including all aligned region, there is one sig outlier on PRA_chr01_94533311_98968251, but this had a short aligned bp

# For following testing filter and exclude the following:
# - all aligned bp < 100,000bp
# - SNP count = 0

# 7 regions were removed from non adaptive (both DEB and PRA pooled)
# 0 region were removed from adaptive


# Local laptop

setwd("/Users/yueyu/Desktop/INV_Chracter")

adp <- read.delim("Adaptive_pi.txt", header = T)
non_adp <- read.delim("NON_Adaptive_pi.txt", header = T)

adp_pi<- adp$Pi_SNP
non_adp_pi <- as.numeric(non_adp$Pi_SNP)

summary(adp_pi)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.004368 0.005498 0.006776 0.006856 0.008330 0.011193  

summary(non_adp_pi, na.rm = TRUE)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.000016 0.000676 0.001079 0.001381 0.001621 0.009655  

# Plotting results
library(ggplot2)

# prepare data
df <- data.frame(
  pi = c(adp_pi, non_adp_pi),
  group = factor(c(rep("Adaptive", length(adp_pi)),
                   rep("Non-Adaptive", length(non_adp_pi))),
                 levels = c("Non-Adaptive", "Adaptive"))
)

# compute p-value (Welch t-test)
pval <- wilcox.test(adp_pi, non_adp_pi)$p.value
#  2.706253e-06


# convert to significance starspva
stars <- ifelse(pval < 0.001, "***",
         ifelse(pval < 0.01, "**",
         ifelse(pval < 0.05, "*", "ns")))

# y position for annotation
y_max <- max(df$pi)

ggplot(df, aes(x = group, y = pi)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  
  # significance bar + star
  annotate("segment", x = 1, xend = 2, y = y_max*1.05, yend = y_max*1.05) +
  annotate("text", x = 1.5, y = y_max*1.1, label = stars, size = 6) +
  
  labs(x = NULL, y = "Divergence (pi)") +
  
  theme_classic() +
  theme(
    axis.text = element_text(size = 18, color = "black"),
    axis.title = element_text(size = 20, color = "black"),
    plot.title = element_text(hjust = 0.5)
  )




#####################################
#
#
# Part 5: Centromere position
#
#
#####################################

# Pericentric (INV span centromere) --> less recombination, less gene flow -> more divergent -> more polymorphic 
# Paracentric (INV outside centro). --> more recombination -> less divergent


# Run on local R studio
# R version 4.5.2

# --- Step 1: load INV positions
setwd("/Users/yueyu/Desktop/INV_Chracter/Raw_INV_info")

library(dplyr)

adp <- read.delim("lostruct_11_INV.txt", header = T)
non_adp <- read.delim("minimap_nonadapt_75_INV.txt", header = T)

non_adp_PRAH1H2 <- read.delim("minimap_pra1and2_nonadapt_53_INV.txt", header = T)
non_adp_DEBH1H2 <- read.delim("minimap_deb1and2_all_and_nonadapt_17_INV.txt", header = T)

adp <- adp[,c(3,4,5)]

non_adpt_pra <- non_adp[,c(1,2,3)]
non_adpt_pra <- rbind(non_adpt_pra, non_adp_PRAH1H2[,c(5,6,7)])
colnames(non_adpt_pra) <- c("CHR","START","END")

non_adpt_deb <- non_adp_DEBH1H2[,c(5,6,7)]
colnames(non_adpt_deb) <- c("CHR","START","END")


non_adpt_pool <- rbind(non_adpt_deb,non_adpt_pra)

# --- Step 2: load Centromere position

# DEB Hap 2
/Users/yueyu/Desktop/circo_2026March11/centromere/deb_h2_centromere_ready_circo.txt
# PRA Hap 2
/Users/yueyu/Desktop/circo_2026March11/centromere/pra_h2_centromere_ready_circo.txt


cetro_DEB <- read.delim("/Users/yueyu/Desktop/circo_2026March11/centromere/deb_h2_centromere_ready_circo.txt", header = FALSE)
cetro_PRA <- read.delim("/Users/yueyu/Desktop/circo_2026March11/centromere/pra_h2_centromere_ready_circo_CHR_changed.txt", header = FALSE)


colnames(cetro_DEB) <- c("CHR","START","END")
colnames(cetro_PRA) <- c("CHR","START","END")

cetro_PRA$CHR <- gsub("_H2","",cetro_PRA$CHR )

cetro_Pooled <- rbind(cetro_DEB,cetro_PRA)


# --- Step 3: Overlap


# --- Adaptive pericentric

df <- adp

df$pericentric <- "NA"

for (i in 1:nrow(df)) {

	chr <- df$CHR[i]
	start <- df$START[i]
	end <- df$END[i]

	a <- cetro_Pooled %>% filter(CHR == chr & START > start & START < end & END > start & END < end ) 
	
	if(nrow(a) > 0){

	df$pericentric[i] <- "YES"

		} else{
			 df$pericentric[i] <- "NO" }

	print(i)
}


df$length <- df$END - df$START

adapt_df <- df
          CHR     START       END pericentric   length
1   DEB_Chr10  98912144 122472701          NO 23560557
2   DEB_Chr12  16080928  27989816          NO 11908888
3   PRA_chr11  14700032  58660853          NO 43960821
4   PRA_chr17 150634096 158777061          NO  8142965
5   PRA_chr17 104981731 134411976          NO 29430245
6   PRA_chr08  53424150  72008140          NO 18583990
7   PRA_chr08  17075757  42752285          NO 25676528
8   PRA_chr05  57631761  92638308         YES 35006547
9  PRA_chr17  105179356 135000890          NO 29821534
10  PRA_chr08  56732032  65046222          NO  8314190
11  PRA_chr11  16772775  41546661          NO 24773886




# --- NON Adaptive pericentric

df <- non_adpt_pool

df$pericentric <- "NA"


for (i in 1:nrow(df)) {

	chr <- df$CHR[i]
	start <- df$START[i]
	end <- df$END[i]

	a <- cetro_Pooled %>% filter(CHR == chr & START > start & START < end & END > start & END < end ) 
	
	if(nrow(a) > 0){

	df$pericentric[i] <- "YES"

		} else{
			 df$pericentric[i] <- "NO" }

	print(i)
}


df$length <- df$END - df$START

NONadapt_df <- df

          CHR     START       END pericentric   length
1   DEB_Chr09  51758708  64362381         YES 12603673
2   DEB_Chr03  66773446  74566342          NO  7792896
3   DEB_Chr15    875106   3871324          NO  2996218
4   DEB_Chr04   3538315   6316414          NO  2778099
5   DEB_Chr07 205525604 208851832          NO  3326228
6   DEB_Chr11 203667379 205398661          NO  1731282
7   DEB_Chr02 163616830 165120067          NO  1503237
8   DEB_Chr17 123545173 125058147          NO  1512974
9   DEB_Chr10 156637299 157534062          NO   896763
10  DEB_Chr05  84139687  85351646          NO  1211959
11  DEB_Chr12 191766494 192709306          NO   942812
12  DEB_Chr03 151600460 152518442          NO   917982
13  DEB_Chr04 212949765 213863459          NO   913694
14  DEB_Chr04 202032500 202778725          NO   746225
15  DEB_Chr06 190079703 190814744          NO   735041
16  DEB_Chr15  73606299  74269695          NO   663396
17  DEB_Chr06 166520086 167073591          NO   553505
18  PRA_chr12  91222312 188019674         YES 96797362
19  PRA_chr08  95760342 138734423         YES 42974081
20  PRA_chr07  19091258  49532268         YES 30441010
21  PRA_chr02  51720243  78879917         YES 27159674
22  PRA_chr04  47653060  73820251          NO 26167191
23  PRA_chr13 186733118 212503654          NO 25770536
24  PRA_chr13  42068946  66183613          NO 24114667
25  PRA_chr13 125326668 148095206         YES 22768538
26  PRA_chr05 133248880 153891714          NO 20642834
27  PRA_chr10  46106924  62324115          NO 16217191
28  PRA_chr17  91768030 107317857          NO 15549827
29  PRA_chr16  22136938  36484670          NO 14347732
30  PRA_chr03  36123814  49598159          NO 13474345
31  PRA_chr01  15852629  28686802          NO 12834173
32  PRA_chr16 136800539 148362028          NO 11561489
33  PRA_chr13 100139599 111563920          NO 11424321
34  PRA_chr09  49761166  60199104         YES 10437938
35  PRA_chr06 193567641 203715645          NO 10148004
36  PRA_chr06 103945299 113670178          NO  9724879
37  PRA_chr16   6081805  15204166          NO  9122361
38  PRA_chr14 100713021 109305575          NO  8592554
39  PRA_chr07  85450125  93143557          NO  7693432
40  PRA_chr03  25902616  33387718          NO  7485102
41  PRA_chr06 180535267 187914763          NO  7379496
42  PRA_chr07   9950248  16999383          NO  7049135
43  PRA_chr13  91914087  98716979         YES  6802892
44  PRA_chr01  87382987  93972976          NO  6589989
45  PRA_chr01   9125193  15670296          NO  6545103
46  PRA_chr06 132442220 138778440         YES  6336220
47  PRA_chr04  28783926  34863611          NO  6079685
48  PRA_chr06  96521089 102256529          NO  5735440
49  PRA_chr06 188129072 193530740          NO  5401668
50  PRA_chr17  30417858  34460548         YES  4042690
51  PRA_chr11   2359702   6083182          NO  3723480
52  PRA_chr07 126584317 130068026          NO  3483709
53  PRA_chr09  65892747  69346575          NO  3453828
54  PRA_chr03  54925280  57685213          NO  2759933
55  PRA_chr03 128675604 131407941          NO  2732337
56  PRA_chr09   3278248   5881451          NO  2603203
57  PRA_chr14  78877391  81450687          NO  2573296
58  PRA_chr11  64180032  66538963          NO  2358931
59  PRA_chr12  11673452  14025449          NO  2351997
60  PRA_chr04 191278947 193615305          NO  2336358
61  PRA_chr05  94805421  97126217          NO  2320796
62  PRA_chr07 162992101 165304046          NO  2311945
63  PRA_chr12 208755041 210934129          NO  2179088
64  PRA_chr11 102570504 104557636          NO  1987132
65  PRA_chr08 169644358 171607615          NO  1963257
66  PRA_chr14  43398806  44903309          NO  1504503
67  PRA_chr02  87318738  88805310          NO  1486572
68  PRA_chr16 159263122 160610767          NO  1347645
69  PRA_chr03 181333404 182671936          NO  1338532
70  PRA_chr16  83520484  84834561          NO  1314077
71  PRA_chr16  82228655  83368827          NO  1140172
72  PRA_chr06  80455663  81564445          NO  1108782
73  PRA_chr06 121004889 122053844          NO  1048955
74  PRA_chr11 124808546 125840371          NO  1031825
75  PRA_chr07 134725385 135730971          NO  1005586
76  PRA_chr10  91656719  92635435          NO   978716
77  PRA_chr15 104673946 105634226          NO   960280
78  PRA_chr10  44099164  45059419          NO   960255
79  PRA_chr16 181440509 182314829          NO   874320
80  PRA_chr15  39467981  40323205          NO   855224
81  PRA_chr06  64567578  65329308          NO   761730
82  PRA_chr05  97549801  98290639          NO   740838
83  PRA_chr10 154883454 155578589          NO   695135
84  PRA_chr07 155474330 156150842          NO   676512
85  PRA_chr07  54221009  54884717          NO   663708
86  PRA_chr17 162652653 163312830          NO   660177
87  PRA_chr07  93556324  94182581          NO   626257
88  PRA_chr09  19806181  20399106          NO   592925
89  PRA_chr05  98655152  99235526          NO   580374
90  PRA_chr10   6034571   6614687          NO   580116
91  PRA_chr14  16387840  16947815          NO   559975
92  PRA_chr14  55082292  55589483          NO   507191
93  PRA_chr16  22171879  36501814          NO 14329935
94  PRA_chr13  89056578  99566255         YES 10509677
95  PRA_chr06 103934593 113802882          NO  9868289
96  PRA_chr14 100713021 108975803          NO  8262782
97  PRA_chr07  64154687  70930521          NO  6775834
98  PRA_chr06 131846155 138568298         YES  6722143
99  PRA_chr16  15315101  21907286          NO  6592185
100 PRA_chr06  97108441 102255767          NO  5147326
101 PRA_chr01  94533311  98968251          NO  4434940
102 PRA_chr13 101050199 105138228          NO  4088029
103 PRA_chr04  83461063  87470552          NO  4009489
104 PRA_chr04  56537206  60484555          NO  3947349
105 PRA_chr11   2480666   6080360          NO  3599694
106 PRA_chr15 120143941 123548428          NO  3404487
107 PRA_chr06 119390854 122484748          NO  3093894
108 PRA_chr01  91350570  94325084          NO  2974514
109 PRA_chr05 175373292 178240763          NO  2867471
110 PRA_chr13 124689482 127524316          NO  2834834
111 PRA_chr12 133553267 136244930          NO  2691663
112 PRA_chr05  94445022  97103529          NO  2658507
113 PRA_chr14  78876327  81401133          NO  2524806
114 PRA_chr12 208774911 210928795          NO  2153884
115 PRA_chr11  64181514  66269183          NO  2087669
116 PRA_chr11 102778652 104689912          NO  1911260
117 PRA_chr12 102220248 103975992          NO  1755744
118 PRA_chr12 180886716 182620609          NO  1733893
119 PRA_chr15  96572043  98169753          NO  1597710
120 PRA_chr07 134487093 136040043          NO  1552950
121 PRA_chr16  99536602 100924999          NO  1388397
122 PRA_chr06 194977625 196336101          NO  1358476
123 PRA_chr16 159173306 160460498          NO  1287192
124 PRA_chr16  45575527  46817834          NO  1242307
125 PRA_chr14 121101932 122327450          NO  1225518
126 PRA_chr11  85345119  86552341         YES  1207222
127 PRA_chr07  45389174  46583753          NO  1194579
128 PRA_chr16  83697211  84890767          NO  1193556
129 PRA_chr02  87021975  88100397          NO  1078422
130 PRA_chr16  82200767  83261719          NO  1060952
131 PRA_chr16 181287705 182315248          NO  1027543
132 PRA_chr12 138343500 139343396         YES   999896
133 PRA_chr10  91630141  92536430          NO   906289
134 PRA_chr01 101279972 102181104          NO   901132
135 PRA_chr13 162060953 162945525          NO   884572
136 PRA_chr07  93438572  94288433          NO   849861
137 PRA_chr10  44236423  45058535          NO   822112
138 PRA_chr06  64526989  65329652          NO   802663
139 PRA_chr17   7716821   8512098          NO   795277
140 PRA_chr16 178111325 178892531          NO   781206
141 PRA_chr07  83061081  83801754          NO   740673
142 PRA_chr07  43587394  44276252          NO   688858
143 PRA_chr14  43839357  44357729          NO   518372
144 PRA_chr11  90104829  90617466          NO   512637
145 PRA_chr17 162787480 163299877          NO   512397







# --- Step 4 (Version A): Compare & Plot (Simply binary comparison)

> df$pericentric
  [1] "YES" "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO" 
 [16] "NO"  "NO"  "YES" "YES" "YES" "YES" "NO"  "NO"  "NO"  "YES" "NO"  "NO"  "NO"  "NO"  "NO" 
 [31] "NO"  "NO"  "NO"  "YES" "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "YES" "NO"  "NO" 
 [46] "YES" "NO"  "NO"  "NO"  "YES" "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO" 
 [61] "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO" 
 [76] "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO" 
 [91] "NO"  "NO"  "YES" "NO"  "NO"  "NO"  "YES" "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO" 
[106] "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO" 
[121] "NO"  "NO"  "NO"  "NO"  "YES" "NO"  "NO"  "NO"  "NO"  "NO"  "YES" "NO"  "NO"  "NO"  "NO" 
[136] "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO" 
> adapt_df$pericentric
 [1] "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "NO"  "YES" "NO"  "NO"  "NO"


tab <- table(
  group = c(rep("df", length(df$pericentric)),
            rep("adp", length(adapt_df$pericentric))),
  
  pericentric = c(df$pericentric,
                  adapt_df$pericentric)
)

tab
# Contingency table
     pericentric
group  NO YES
  adp  10   1
  df  131  14

fisher.test(tab)

#		Fisher's Exact Test for Count Data

data:  tab
p-value = 1
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
  0.1327434 49.6500407
sample estimates:
odds ratio 
  1.068269 

# CONCLUSION: no significance differnce between two groups





# --- Step 4 (Version B): Compare & Plot (binary comparison + length) --> logistic regression


# Add group labels
adapt_df$type <- "adapt"
NONadapt_df$type <- "NONadapt"

# Combine data
df <- rbind(adapt_df, NONadapt_df)

# Convert variables
df$pericentric_bin <- ifelse(df$pericentric == "YES", 1, 0)
df$type <- factor(df$type, levels = c("NONadapt", "adapt"))

# Optional: log-transform length because genomic lengths are highly skewed
df$log_length <- log10(df$length)

# Logistic regression
model <- glm(pericentric_bin ~ type + log_length,
             data = df,
             family = binomial)

summary(model)


Coefficients:
            Estimate Std. Error z value Pr(>|z|)    
(Intercept) -20.9405     4.7794  -4.381 1.18e-05 ***
typeadapt    -2.0291     1.1676  -1.738   0.0822 .  
# adaptive regions have no significantly lower odds of being pericentric than NONadapt regions, after accounting for interval length

log_length    2.8032     0.6924   4.048 5.16e-05 ***
# longer intervals are much more likely to overlap pericentromeric regions

---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 98.763  on 155  degrees of freedom
Residual deviance: 75.880  on 153  degrees of freedom
AIC: 81.88
