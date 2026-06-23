#####################################
# Project: debilis_praecox_project
#
# Code 7: GEA
#
# by: Yue Yu
#
#####################################

# Using cluster C as an example. Same code can be applied to the remaining clusters

# --------------------------
#   Part 1: AF per pop
# --------------------------
# -- Step 1: subset VCF by populations --> freq file
module load bcftools

cd /home/yueyu/scratch/Texas_clusterC
mkdir LFMM

cd LFMM
VCF="/home/yueyu/scratch/Texas_clusterC/SNP_filter/GBS_TEXAS_SNP_INFO_GENO_CLUSTER_C_BI_NOCALL03_FILTERED_126samples_CLEANED_AF003.vcf.gz"  

bcftools query -l $VCF > C_126_sample_names.txt
# 126 sample names

cat C_126_sample_names.txt | awk -F'_' '{print$1}' | sort | uniq > C_13_pop_names.txt
# 13 populations
# PH 2
# PR 1-3
# DC 2,3,5,6
# DS 1-5


# --  Make a ind file for each population -> match VCF file names
for i in {1..13}; do
	POP_NAME=$(cat C_13_pop_names.txt | head -n $i | tail -1)
	cat C_126_sample_names.txt | grep $POP_NAME > ${POP_NAME}_ind_names.txt
done



# -- VCF subset
module load vcftools

mkdir freq

for i in {1..13}; do
	POP_NAME=$(cat C_13_pop_names.txt | head -n $i | tail -1)
	vcftools --gzvcf $VCF --keep ${POP_NAME}_ind_names.txt --freq --out freq/${POP_NAME}_freq
done

cd freq
rm *.log


# -- Step 2: use R to calculate ALT ALLELE FREQ for each pop (batch process)
cd /home/yueyu/scratch/Texas_clusterC/LFMM/freq

module load StdEnv/2023
module load r/4.4.0

R

library(dplyr)
library(tidyr)

# List all .frq files ending with _freq.frq in the working directory
frq_files <- list.files(pattern = "_freq\\.frq$")

# Initialize empty list to store data frames
maf_list <- list()

# Loop over files
for (file in frq_files) {
  
  # Create label: e.g., pop1 from pop1_freq.frq
  label <- sub("_freq\\.frq$", "", file)
  
  # Read file, skip header
  df <- read.table(file, header = FALSE, skip = 1, sep = "\t", stringsAsFactors = FALSE)
  
  # Columns:
  # V1 = CHROM
  # V2 = POS
  # V3 = N_ALLELES
  # V4 = N_CHR
  # V5 = REF allele:frequency
  # V6 = ALT allele:frequency - THIS WILL BE EXTRACTED AND USED 
  
  # the allele columns --> long
  df_long <- df %>%
    pivot_longer(cols = 5:ncol(df), names_to = "AlleleCol", values_to = "AlleleFreq") %>%
    separate(AlleleFreq, into = c("Allele", "Freq"), sep = ":", convert = TRUE) %>%
    group_by(CHROM = V1, POS = V2) %>%
    summarize(!!paste0("MAF_", label) := min(Freq, na.rm = TRUE)) %>%
    ungroup()
  
  maf_list[[label]] <- df_long 

  # 算出的是alt等位基因的频率
}


# Full join all data frames by CHROM and POS
maf_merged <- Reduce(function(x, y) full_join(x, y, by = c("CHROM", "POS")), maf_list)

dim(maf_merged)
# 36512    15

# Write output
write.table(maf_merged, "C_ALT_AF_20260511.txt", row.names = FALSE, quote = FALSE, sep = "\t")


# --------------------------
#  Part 1.1: MAF avergaed among all pops
# --------------------------
MAF <- read.table("C_ALT_AF_20260511.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(MAF)

# Replace Inf/-Inf with NA in columns 3 to end
MAF[ , 3:ncol(MAF)][!is.finite(as.matrix(MAF[ , 3:ncol(MAF)]))] <- NA

# Now calculate row-wise mean across columns 3 to end
MAF$MAF_mean <- round(rowMeans(MAF[ , 3:ncol(MAF)], na.rm = TRUE), 5)
# 36512    16

summary(MAF$MAF_mean)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.00000 0.04565 0.08278 0.11627 0.17308 0.45192

write.table(MAF, "C_ALT_AF_withAVERAGE_20260511.txt", row.names = FALSE, quote = FALSE, sep = "\t")



# --------------------------
#   Part 2: RUN LFMM -- Notes all in "GBS/PART12_LFMM"
# --------------------------
tmux new-session -s lfmm
tmux attach-session -t lfmm

cd /home/yueyu/scratch/Texas_clusterC/LFMM/freq

module load r

R
# R version 4.5.0

library(lfmm)
library(dplyr)
library(ggplot2)


# ======= Step 1: Load AF data
af <- read.delim("C_ALT_AF_withAVERAGE_20260511.txt", header=TRUE)
colnames(af) <- gsub("MAF_","",colnames(af))

af_maf5 <- af[af$mean >= 0.05,]
dim(af_maf5)
# 26,000    16

sample_ids <- colnames(af[,-c(1,2,ncol(af))])
sample_ids
#  [1] "DC2" "DC3" "DC5" "DC6" "DS1" "DS2" "DS3" "DS4" "DS5" "PH2" "PR1" "PR2"
# [13] "PR3"

# ======= Step 2: Load Climate data
# Load and align environmental data
clim <- read.table("/home/yueyu/scratch/Climate_for_GBS/Coord_forCLIMATE_2026Feb_Normal_1961_1990_2026Feb13.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
names <- readLines("/home/yueyu/scratch/Texas_clusterC/LFMM/C_13_pop_names.txt")

# Subset whe:qre V2 is in names
clim_sub <- clim[clim$ID2 %in% names, ]
dim(clim_sub)
clim_ready <- clim_sub[,-1]
clim_ready 

clim_ready[1:5,1:10]
dim(clim_ready)

clim_ready$ID2
 [1] "PH2" "PR1" "PR2" "PR3" "DC2" "DC3" "DC5" "DC6" "DS1" "DS2" "DS3" "DS4"
[13] "DS5"

# ======= Step 3: Match Sample Ordering
# Vector to match
desired_order <- sample_ids

# Reorder clim_sub to match desired_order
clim_ready_2 <- clim_ready[match(desired_order, clim_ready$ID2), ]
clim_ready_2$ID2
 [1] "DC2" "DC3" "DC5" "DC6" "DS1" "DS2" "DS3" "DS4" "DS5" "PH2" "PR1" "PR2"
[13] "PR3"
# Now matching!!


# ======= Step 4: Best K determined (Cluster C)
# K determined by Admixture best K
# Test K (+/-1) around optimal K, may result in bets candidate detection
# Note that with K=1, LFMM is essentially running a simple linear regression (i.e., no latent factors).


# ======= CLIMATE ANNUAL + SEASONAL IN ONE LOOP ===========
clim_ready_3 <- clim_ready_2[,c(1,2,3,5,6,7,8,9,10,11,12,13,14,15,16,18,19,23,24,25,27,29,30:45,50:57,66:69,74:85)]
dim(clim_ready_3) 
# 13 pops  &  62 = 19 ANNUAL CLIM + 40 SEASONAL CLIM + 3 IDs

colnames(clim_ready_3)
 [1] "ID2"     "lat"     "long"    "MAT"     "MWMT"    "MCMT"    "TD"
 [8] "MAP"     "MSP"     "AHM"     "SHM"     "bFFP"    "eFFP"    "FFP"
[15] "CMD"     "DD_0"    "DD5"     "EMT"     "EXT"     "Eref"    "NFFD"
[22] "RH"      "Tmax_wt" "Tmax_sp" "Tmax_sm" "Tmax_at" "Tmin_wt" "Tmin_sp"
[29] "Tmin_sm" "Tmin_at" "Tave_wt" "Tave_sp" "Tave_sm" "Tave_at" "PPT_wt"
[36] "PPT_sp"  "PPT_sm"  "PPT_at"  "DD_0_wt" "DD_0_sp" "DD_0_sm" "DD_0_at"
[43] "DD5_wt"  "DD5_sp"  "DD5_sm"  "DD5_at"  "NFFD_wt" "NFFD_sp" "NFFD_sm"
[50] "NFFD_at" "Eref_wt" "Eref_sp" "Eref_sm" "Eref_at" "CMD_wt"  "CMD_sp"
[57] "CMD_sm"  "CMD_at"  "RH_wt"   "RH_sp"   "RH_sm"   "RH_at"




# ======= Step 5:Run LFMM
sum(is.na(af_maf5)) # 806 missing data with "NA"
af_maf5_clean <- na.omit(af_maf5)
rownames(af_maf5_clean) <- paste0(af_maf5_clean$CHROM,"_",af_maf5_clean$POS)



# ==== Transform AF before running
af_maf5_clean_2 <- t(af_maf5_clean[,-c(1,2,ncol(af_maf5_clean))])
af_maf5_clean_2[1:5,1:10]
dim(af_maf5_clean_2)
# 13 25245 ( > 26,000 - 806 data with NA, some row with more than 1 NA counted )

# ==== Set parameters 
bon_thershold <- -log10(0.05/ncol(af_maf5_clean_2));bon_thershold  #5.703205
best_K <- 6 


# ==== Testing best K & MAT

best_K <- 1 # git = 1.96
best_K <- 2 # git = 1.91  *** best K indicated by git score <--- used this for final LFMM result
best_K <- 3 # git = 2.03
best_K <- 4 # git = 2.25
best_K <- 5 # git = 2.4
best_K <- 6 # git = 2.7.  *** best K indictaed by ADMIXTURE
best_K <- 7 # git = 3.04
best_K <- 8 # git = 3.73



# ==== Start Loop: only loop through all CLIM, not through all LOCI

# empty df for outliers
outlier_df <- data.frame(
  CHROM = character(),
  WINDOW = character(),
  logP = numeric(),
  var = character(),
  stringsAsFactors = FALSE)

# empty df for outliers
to_plot_p <- data.frame( snp = colnames(af_maf5_clean_2), stringsAsFactors = FALSE)

c <- 23 # Tmax_wt

# Run in loop per climate variable
for (c in 4:62) {

      clim <- colnames(clim_ready_3)[c]
      print(clim)

      # --------- Step 1: Run LFMM  ---------
      Y <- af_maf5_clean_2
      X <- clim_ready_3[,c]

      # LFMM estimation
      macro.lfmm <- lfmm_ridge(Y = Y, X = X, K = best_K)   

      # Calc P value
      macro.pvalues <- lfmm_test(Y = Y, X = X, lfmm = macro.lfmm, calibrate = 'gif')
      print(macro.pvalues$gif)
      
      # outliers.macro.Fetch <- colnames(af_maf5_clean_2)[which(macro.pvalues$calibrated.pvalue < bon_thershold_beforelog)]

    if (macro.pvalues$gif > 0) {

     # --------- Step 2: Log Transform  ---------

        to_plot_p$logP <- -log10(macro.pvalues$calibrated.pvalue)

        to_plot_p$CHROM <- sapply(strsplit(to_plot_p$snp, "_"), `[`, 2)
        to_plot_p$CHROM <- gsub("chr0","",to_plot_p$CHROM) # Change based on CHR format
        to_plot_p$CHROM <- gsub("chr","",to_plot_p$CHROM)
        to_plot_p$CHROM <- as.numeric(to_plot_p$CHROM)

        to_plot_p$WINDOW <- sapply(strsplit(to_plot_p$snp, "_"), `[`, 3)
        to_plot_p$WINDOW <- as.numeric(to_plot_p$WINDOW)

      # --------- Step 3: Save outliers in df ---------
        lets_save <- to_plot_p[to_plot_p$logP >= bon_thershold,c("CHROM","WINDOW","logP")]
        lets_save$var <- rep(clim, nrow(lets_save))

        if(nrow(lets_save) > 0){outlier_df <- rbind(outlier_df, lets_save)} 


      # --------- Step 4: Plot GEA Manplot  ---------

        # ADD CUMULATIVE POSITIONS FOR PLOTTING
        to_plot_p_cum <- to_plot_p %>%
          
          # Compute chromosome size
          group_by(CHROM) %>%
          summarise(chr_len=max(WINDOW)) %>%
          
          # Calculate cumulative position of each chromosome
          mutate(tot=cumsum(chr_len)-chr_len) %>%
          select(-chr_len) %>%
          
          # Add this info to the initial dataset
          left_join(to_plot_p, ., by=c("CHROM"="CHROM")) %>%
          
          # Add a cumulative position of each SNP
          arrange(CHROM, WINDOW) %>%
          mutate(WINDOW_cum=WINDOW+tot)

        lets_plot <- to_plot_p_cum[,c("CHROM","WINDOW_cum","logP")]


        # Assuming CHR is a character variable
        lets_plot$point_color <- ifelse(lets_plot$logP >= bon_thershold, "black", 
                            ifelse(lets_plot$CHROM %% 2 == 1, "cornflowerblue", "darkorange1"))
        #Need to be a factor to be plotted!!
        lets_plot$point_color <- factor(lets_plot$point_color,levels = c("cornflowerblue", "darkorange1", "black"))


        # -- define X axis (after CHR is numeric)
        axisdf = lets_plot %>% group_by(CHROM) %>% summarize(center=( max(WINDOW_cum) + min(WINDOW_cum) ) / 2 )

        plot <- ggplot(lets_plot, aes(x = WINDOW_cum, y = logP, color = point_color)) +
          geom_point(alpha = 0.8, size = 2.7) +
          scale_color_manual(values = c("cornflowerblue", "darkorange1", "black")) +
          geom_hline(yintercept = bon_thershold, linetype = "dashed", color = "black", linewidth = 1.6) +
          scale_x_continuous(labels = axisdf$CHROM, breaks = axisdf$center) +
          scale_y_continuous(limits = c(0, NA)) +
          theme_bw() +
          theme(
            legend.position = "none",
           panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            strip.background = element_blank(),
            axis.text = element_text(size = 20, colour = "black"),  #axis text and title colour black
            axis.title = element_text(size = 24, colour = "black"),
            axis.title.x = element_text(margin = margin(t = 15), colour = "black"), # Adjust t (top margin) as needed
            axis.title.y = element_text(margin = margin(r = 15), colour = "black")  # Adjust r (right margin) as needed
          ) +
          labs(
            y = paste0("-log10(P) for ", clim),
            x = "Chromosome")
         
        png(paste0("plots_2026May11/",clim,"_LFMM_20260511.png"),width = 2000, height = 500, res = 100)
        print(plot)
        dev.off()

        #pdf(file = paste0("plots_2026May11/",clim,"_LFMM_20260511.pdf"),width = 20, height = 5)
        #print(plot)
        #dev.off()
        # Done 2026 May 26 for PDF plots

    } # End of if gif > 0 loop     

  } #  End of c loop


# Output 1: Manhattan plot of all P value, already saved
# Output 2: outlier_df (recorded all outliers --> Need to save)
head(outlier_df)
dim(outlier_df)
# 5195    4

table(outlier_df$var)
write.table(outlier_df, file = "plots_2026May11/outlier_df_20260511.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
# SAVED 2026 May 11




# ---- Find overlap with within-pop inversions
# INV position
inv <- data.frame(
  CHROM = c(17, 8, 8),
  START = c(104981731, 53424150, 17075757 ),
  END = c(134411976, 72008140, 42752285),
  NAME = c("pra.17.02", "pra.08.01","pra.08.02") )

inv
1    17 104981731 134411976 pra.17.02
2     8  53424150  72008140 pra.08.01
3     8  17075757  42752285 pra.08.02


# Read in GEA signals
# outlier_df <- read.delim("plots_2026May11/outlier_df_20260511.txt")
outlier_df$INV_NAME <- NA


# Empty df to store annotated iNV info
a <- data.frame(
  CHROM = character(),
  WINDOW = character(),
  logP = numeric(),
  var = character(),
  INV_NAME = character(),
  stringsAsFactors = FALSE)

# Run in loop through all inversions
i <- 1

for (i in 1:nrow(inv)) {
  outlier_df[outlier_df$CHROM == inv$CHROM[i] & outlier_df$WINDOW <= inv$END[i] & outlier_df$WINDOW >= inv$START[i], 5] <- inv[i,4]
  a <- rbind(a, outlier_df[outlier_df$CHROM == inv$CHROM[i] & outlier_df$WINDOW <= inv$END[i] & outlier_df$WINDOW >= inv$START[i],])

}

a
dim(a)
#  188   5

table(a$INV_NAME)

# pra.08.01 pra.08.02 pra.17.02
#        2        68        50

write.table(a, file = "plots_2026May11/outlier_df_inside_INV_20260522.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")













































# ----------------------------
#
# ------- SOIL -------------- # ALL pop have SOIL data
#
# ----------------------------

# Load climate all
soil <- read.table("/home/yueyu/scratch/Soil_for_GBS/SSURGO_extracted_soil_muagg_and_chorizon_2026Feb24_forGEA.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Load pops
names <- readLines("/home/yueyu/scratch/Texas_clusterC/LFMM/C_13_pop_names.txt")
# 13 pops

# Subset whe:qre V2 is in names
soil_sub <- soil[soil$ID2 %in% names, ]
dim(soil_sub)

# Subset col (samples) in af_maf5_clean
c <- c("CHROM","POS",soil_sub$ID2,"mean")


af <- read.delim("C_ALT_AF_withAVERAGE_20260511.txt", header=TRUE)
colnames(af) <- gsub("MAF_","",colnames(af))

af_maf5 <- af[af$mean >= 0.05,]
dim(af_maf5)

sum(is.na(af_maf5)) # missing data with "NA"
af_maf5_clean <- na.omit(af_maf5)
rownames(af_maf5_clean) <- paste0(af_maf5_clean$CHROM,"_",af_maf5_clean$POS)

af_maf5_clean <- af_maf5_clean[,colnames(af_maf5_clean) %in% c]
dim(af_maf5_clean)
# 25245    16 = 13 pop + IDs

colnames(af_maf5_clean)
 [1] "CHROM" "POS"   "DC2"   "DC3"   "DC5"   "DC6"   "DS1"   "DS2"   "DS3"
[10] "DS4"   "DS5"   "PH2"   "PR1"   "PR2"   "PR3"   "mean"

# Clim ready: remove unneeded col + reorder POP names to match GT
soil_ready <- soil_sub[,-c(1,5,6,7,8,9,10)] # Remove unused columns

soil_ready[1:5,1:10]
dim(soil_ready)

# Vector to match
desired_order <- colnames(af_maf5_clean)[3:15]

# Reorder clim_sub to match desired_order
soil_ready <- soil_ready[match(desired_order, soil_ready$ID2), ]
soil_ready$ID2
#  [1] "DC2" "DC3" "DC5" "DC6" "DS1" "DS2" "DS3" "DS4" "DS5" "PH2" "PR1" "PR2"
# [13] "PR3"



# ======= Step 5:Run LFMM

# ==== Transform AF before running
af_maf5_clean_2 <- t(af_maf5_clean[,-c(1,2,ncol(af_maf5_clean))])
dim(af_maf5_clean_2)
#  13 25245


# ==== Set parameters 
bon_thershold <- -log10(0.05/ncol(af_maf5_clean_2));bon_thershold #5.703205
best_K <- 2



# ==== Start Loop: only loop through all CLIM, not through all LOCI

# empty df for outliers
outlier_df <- data.frame(
  CHROM = character(),
  WINDOW = character(),
  logP = numeric(),
  var = character(),
  stringsAsFactors = FALSE)

# empty df for outliers
to_plot_p <- data.frame( snp = colnames(af_maf5_clean_2), stringsAsFactors = FALSE)

c <- 4 # claytotal_r_0_180cm

# Run in loop per variable

for (c in 4:18) {

      clim <- colnames(soil_ready)[c]
      print(clim)

      # --------- Step 1: Run LFMM  ---------
      Y <- af_maf5_clean_2
      X <- soil_ready[,c]

      # If any climate variable had NAs, remove that pop
      if(any(is.na(X)) == TRUE) {
          
          na_idx <- is.na(X)   
          #which_na <- which(is.na(X))

          # remove NA from X
          X <- X[!na_idx]

          # remove corresponding entries in Y
          Y <- Y[!na_idx,]
      } # End of NA check in soil variable


      # LFMM estimation
      macro.lfmm <- lfmm_ridge(Y = Y, X = X, K = best_K)   

      # Calc P value
      macro.pvalues <- lfmm_test(Y = Y, X = X, lfmm = macro.lfmm, calibrate = 'gif')
      print(macro.pvalues$gif)
      
      # outliers.macro.Fetch <- colnames(af_maf5_clean_2)[which(macro.pvalues$calibrated.pvalue < bon_thershold_beforelog)]

    if (macro.pvalues$gif > 0) {

     # --------- Step 2: Log Transform  ---------

        to_plot_p$logP <- -log10(macro.pvalues$calibrated.pvalue)

        to_plot_p$CHROM <- sapply(strsplit(to_plot_p$snp, "_"), `[`, 2)
        to_plot_p$CHROM <- gsub("chr0","",to_plot_p$CHROM)
        to_plot_p$CHROM <- gsub("chr","",to_plot_p$CHROM)
        to_plot_p$CHROM <- as.numeric(to_plot_p$CHROM)

        to_plot_p$WINDOW <- sapply(strsplit(to_plot_p$snp, "_"), `[`, 3)
        to_plot_p$WINDOW <- as.numeric(to_plot_p$WINDOW)

      # --------- Step 3: Save outliers in df ---------
        lets_save <- to_plot_p[to_plot_p$logP >= bon_thershold,c("CHROM","WINDOW","logP")]
        lets_save$var <- rep(clim, nrow(lets_save))

        if(nrow(lets_save) > 0){outlier_df <- rbind(outlier_df, lets_save)} 


      # --------- Step 4: Plot GEA Manplot  ---------

        # ADD CUMULATIVE POSITIONS FOR PLOTTING
        to_plot_p_cum <- to_plot_p %>%
          
          # Compute chromosome size
          group_by(CHROM) %>%
          summarise(chr_len=max(WINDOW)) %>%
          
          # Calculate cumulative position of each chromosome
          mutate(tot=cumsum(chr_len)-chr_len) %>%
          select(-chr_len) %>%
          
          # Add this info to the initial dataset
          left_join(to_plot_p, ., by=c("CHROM"="CHROM")) %>%
          
          # Add a cumulative position of each SNP
          arrange(CHROM, WINDOW) %>%
          mutate(WINDOW_cum=WINDOW+tot)

        lets_plot <- to_plot_p_cum[,c("CHROM","WINDOW_cum","logP")]


        # Assuming CHR is a character variable
        lets_plot$point_color <- ifelse(lets_plot$logP >= bon_thershold, "black", 
                            ifelse(lets_plot$CHROM %% 2 == 1, "cornflowerblue", "darkorange1"))
        #Need to be a factor to be plotted!!
        lets_plot$point_color <- factor(lets_plot$point_color,levels = c("cornflowerblue", "darkorange1", "black"))


        # -- define X axis (after CHR is numeric)
        axisdf = lets_plot %>% group_by(CHROM) %>% summarize(center=( max(WINDOW_cum) + min(WINDOW_cum) ) / 2 )

        plot <- ggplot(lets_plot, aes(x = WINDOW_cum, y = logP, color = point_color)) +
          geom_point(alpha = 0.8, size = 2.7) +
          scale_color_manual(values = c("cornflowerblue", "darkorange1", "black")) +
          geom_hline(yintercept = bon_thershold, linetype = "dashed", color = "black", linewidth = 1.6) +
          scale_x_continuous(labels = axisdf$CHROM, breaks = axisdf$center) +
          scale_y_continuous(limits = c(0, NA)) +
          theme_bw() +
          theme(
            legend.position = "none",
           panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            strip.background = element_blank(),
            axis.text = element_text(size = 20, colour = "black"),  #axis text and title colour black
            axis.title = element_text(size = 24, colour = "black"),
            axis.title.x = element_text(margin = margin(t = 15), colour = "black"), # Adjust t (top margin) as needed
            axis.title.y = element_text(margin = margin(r = 15), colour = "black")  # Adjust r (right margin) as needed
          ) +
          labs(
            y = paste0("-log10(P) for ", clim),
            x = "Chromosome")
         
        png(paste0("plots_2026May11/",clim,"_LFMM_20260511.png"),width = 2000, height = 500, res = 100)
        print(plot)
        dev.off()

    } # End of if gif > 0 loop     

  } #  End of c loop

# Output 1: Manhattan plot of all P value, already saved
# Output 2: outlier_df (recorded all outliers --> Need to save)
head(outlier_df)
dim(outlier_df)
# 537   4

table(outlier_df$var)
write.table(outlier_df, file = "plots_2026May11/outlier_df_SOIL_20260511.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")


# ---- Find overlap with within-pop inversions
# INV position
inv <- data.frame(
  CHROM = c(17, 8, 8),
  START = c(104981731, 53424150, 17075757 ),
  END = c(134411976, 72008140, 42752285),
  NAME = c("pra.17.02", "pra.08.01","pra.08.02") )

inv


# Read in GEA signals
#outlier_df <- read.delim("plots_2026May11/outlier_df_SOIL_20260511.txt")
outlier_df$INV_NAME <- NA


# Empty df to store annotated iNV info
a <- data.frame(
  CHROM = character(),
  WINDOW = character(),
  logP = numeric(),
  var = character(),
  INV_NAME = character(),
  stringsAsFactors = FALSE)

# Run in loop through all inversions
i <- 1

for (i in 1:nrow(inv)) {
  outlier_df[outlier_df$CHROM == inv$CHROM[i] & outlier_df$WINDOW <= inv$END[i] & outlier_df$WINDOW >= inv$START[i], 5] <- inv[i,4]
  a <- rbind(a, outlier_df[outlier_df$CHROM == inv$CHROM[i] & outlier_df$WINDOW <= inv$END[i] & outlier_df$WINDOW >= inv$START[i],])

}

a
dim(a)
#  8   5

table(a$INV_NAME)
#  pra.08.01 pra.08.02 pra.17.02
#        1         3         4

write.table(a, file = "plots_2026May11/outlier_df_SOIL_inside_INV_20260522.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")



# END