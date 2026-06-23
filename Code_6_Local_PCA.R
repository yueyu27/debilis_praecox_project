#####################################
# Project: debilis_praecox_project
#
# Code 6: Local PCA
#
# by: Yue Yu
#
#####################################

# Using cluster B as an example. Same code can be applied to the remaining clusters


# ========
#
# Cluster B -- 88 sample + 44,214 SNPs
#
# ======= 

# Saved in scratch at (fir): 
cd /home/yueyu/scratch/Texas_clusterB/SNP_filter

#-rw-r----- 1 yueyu 93K Apr 15 17:04 GBS_TEXAS_SNP_INFO_GENO_CLUSTER_B_BI_NOCALL03_FILTERED_88samples_CLEANED_AF003.vcf.gz.tbi
#-rw-r----- 1 yueyu 18M Apr 15 17:04 GBS_TEXAS_SNP_INFO_GENO_CLUSTER_B_BI_NOCALL03_FILTERED_88samples_CLEANED_AF003.vcf.gz



# -------------------
#
# Step 1: VCF -> BCF (.bcf) -> index it
#
# -------------------

module load StdEnv/2023  
module load gcc/12.3
module load bcftools/1.19


cd /home/yueyu/scratch/Texas_clusterB/local_PCA

# ------ 88 sample + 44,214 SNPmv *s, Keep LD SNPs, used for local PCA ------
VCF="/home/yueyu/scratch/Texas_clusterB/SNP_filter/GBS_TEXAS_SNP_INFO_GENO_CLUSTER_B_BI_NOCALL03_FILTERED_88samples_CLEANED_AF003.vcf.gz"
BCF="cluster_B_20260416.bcf"
bcftools view -O b -o $BCF $VCF
bcftools index $BCF
bcftools view -v snps cluster_B_20260416.bcf | grep -v "^#" | wc -l
# DONE




# -------------------
#
# Step 2: Run lostruct
#
# -------------------
# Fir
tmux new-session -s clusterB_login1
tmux attach-session -t clusterB_login1 

cd /home/yueyu/scratch/Texas_clusterB/local_PCA

module load StdEnv/2023
module load gcc/12.3
module load bcftools/1.19    # Need to load bcftools before starting R
module load r/4.4.0

R


# Ready to use 
library(lostruct)
library(dplyr)
library(tidyr)
library(gridExtra)
library(grid)
library(stringr)
library(readr)
library(ggplot2)
select=dplyr::select



# eigen_windows() requires
# Row: variant
# Col: sample

# -- Step 1: Define parameters
window_size <- 50 # window size = 50 SNPs
k_kept <- 40 # Number of MDS kept
min_windows <- 2 # Minimum number of windows to call outlier regions
max_distance_between_outliers <- 4 # Max distance between oulier windows
# Chnaged 2026 May 22, but did not run, directly assesed outliers, all pass this criteria
n_permutations <- 1000 # Number of permutation for chromosomal clustering test
min_cor <- 0.8 # Correlation theshold for collapsing MDS
clust_pvalue_threshold <- 0.01 # Should be 0.01 usually

# -- Step 2: Load BCF
bcf.file <- "cluster_B_20260416.bcf"
sites <- vcf_positions(bcf.file);sites
win.fn.snp <- vcf_windower(bcf.file, size=window_size, type="snp", sites=sites)

# -- Step 3: Local PCA
snp.pca <- eigen_windows(win.fn.snp, k=2, mc.cores=10)
dim(snp.pca)
# 875 179 
# 875 = Number of windows (Can check through 50 SNPs/window * 875 windows ~ 43,750 SNPs, less than the input SNPs)
# 179 = 88 PRA sample from BCF * 2 PCs each samples + 3 columns: total, lam_1, lam_2


# -- Step 4: Distance matrix for MDS 
pcdist <- pc_dist(snp.pca, mc.cores=10) 
dim(pcdist)
# 875 875
# This steps takes some time to run (need screen)

pcdist[1:5,1:5]


# -- Distance matrix for MDS (remove windows with NA)
na.wins <- is.na(snp.pca[,1]) 
pcdist <- pcdist[!na.wins, !na.wins]
nan.wins <- pcdist[,1]=="NaN"
pcdist <- pcdist[!nan.wins, !nan.wins]

# -- Step 5: MDS
mds <- cmdscale(pcdist, eig=TRUE, k=k_kept) 
mds.coords <- mds$points
colnames(mds.coords) <- paste("MDS coordinate", 1:ncol(mds.coords))
mds.coords[1:5,1:5]
dim(mds.coords)
# 868  40

# -- Step 6: Build empty win.regions dataframe and populate window coordinates
win.regions <- region(win.fn.snp)()
win.regions <- win.regions[!na.wins,][!nan.wins,]
win.regions %>% mutate(mid=(start+end)/2) -> win.regions


# -- Step 7: Add MDS 1-40 values to each window
for (k in 1:k_kept){
  name = paste("mds", str_pad(k, 2, pad = "0"), sep="")
  win.regions$tmp <- "NA"
  win.regions <- win.regions %>% rename(!!name := tmp)
}

for (i in 1:k_kept){
  j = i + 4
  win.regions[,j] <- mds.coords[,i]
}
win.regions$n <- 1:nrow(win.regions)


# -- Step 8: 
mds_pcs <- colnames(win.regions)[5:(ncol(win.regions)-1)]

mds_clustering <- tibble(mds_coord=character(), direction=character(), clust_pvalue=numeric(), outliers=numeric(), n1_outliers=numeric(), 
                         high_cutoff=numeric(), lower_cutoff=numeric(), chr=character())

# 3 SD from the mean across all windows in that dimension of MDS
n_sd <- 3

for (mds_chosen in mds_pcs){
  print(paste("Processing",mds_chosen))
  win.regions %>%
    mutate(the_mds = .data[[mds_chosen]]) %>%
    summarize(sd_mds=sd(the_mds)) %>% pull() -> sd_mds
    mds_high_cutoff <- sd_mds*4
    mds_low_cutoff <- sd_mds*3


  win.regions %>%
    mutate(the_mds = .data[[mds_chosen]]) %>%
    mutate(sd_mds=sd(the_mds)) %>%
    filter(the_mds > (sd_mds*n_sd)) -> pos_windows


  if (nrow(pos_windows) >= min_windows){
    permutations <- matrix(nrow=n_permutations, ncol=1)
    for (i in 1:n_permutations){
      max_1 <- win.regions %>% sample_n(nrow(pos_windows)) %>% group_by(chrom) %>% summarize(count=n()) %>% arrange(desc(count)) %>% head(n=1) %>% 
        ungroup() %>% summarize(sum=sum(count)) %>% pull(sum)
      permutations[i,1] <- max_1
    }
    pos_windows %>% group_by(chrom) %>% summarize(count=n()) %>% arrange(desc(count)) %>% head(n=1) %>% 
      ungroup() %>% summarize(sum=sum(count)) %>% pull(sum) -> sampled_max_1
    pos_windows %>% group_by(chrom) %>% summarize(count=n()) %>% arrange(desc(count)) %>% head(n=1) %>%
      pull(chrom) -> clustered_chr
    x <- as_tibble(permutations) %>% filter(V1 >= sampled_max_1) %>% nrow() 
    pvalue <- (x+1)/(n_permutations+1)
    tmp <- tibble(mds_coord=as.character(mds_chosen), direction=as.character("pos"), clust_pvalue=as.numeric(pvalue), outliers=as.numeric(nrow(pos_windows)),
                  n1_outliers=as.numeric(sampled_max_1), high_cutoff=as.numeric(mds_high_cutoff), lower_cutoff=as.numeric(mds_low_cutoff), chr=as.character(clustered_chr))
    mds_clustering <- rbind(mds_clustering, tmp)
  }else{
    tmp <- tibble(mds_coord=as.character(mds_chosen), direction=as.character("pos"), clust_pvalue=as.numeric(NA), outliers=as.numeric(nrow(pos_windows)),
                  n1_outliers=as.numeric(NA), high_cutoff=as.numeric(NA), lower_cutoff=as.numeric(NA), chr=as.numeric(NA))
    mds_clustering <- rbind(mds_clustering, tmp)
  }


 win.regions %>%
    mutate(the_mds=.data[[mds_chosen]]) %>%
    mutate(sd_mds=sd(the_mds)) %>%
    filter(the_mds < -(sd_mds*n_sd)) -> neg_windows


  if (nrow(neg_windows) >= min_windows){
    permutations <- matrix( nrow=n_permutations, ncol=1)
    for (i in 1:n_permutations){
      max_1 <- win.regions %>% sample_n(nrow(neg_windows)) %>% group_by(chrom) %>% summarize(count=n()) %>% arrange(desc(count)) %>% head(n=1) %>% 
        ungroup() %>% summarize(sum=sum(count)) %>% pull(sum)
      permutations[i,1] <- max_1
    }
    neg_windows %>% group_by(chrom) %>% summarize(count=n()) %>% arrange(desc(count)) %>% head(n=1) %>% 
      ungroup() %>% summarize(sum=sum(count)) %>% pull(sum) -> sampled_max_1
    neg_windows %>% group_by(chrom) %>% summarize(count=n()) %>% arrange(desc(count)) %>% head(n=1) %>%
      pull(chrom) -> clustered_chr
    x <- as_tibble(permutations) %>%  filter(V1 >= sampled_max_1) %>% nrow()
    pvalue <- (x+1)/(n_permutations+1)
    tmp <- tibble(mds_coord = as.character(mds_chosen),direction = as.character("neg"),clust_pvalue = as.numeric(pvalue), outliers=as.numeric(nrow(neg_windows)),
                  n1_outliers=as.numeric(sampled_max_1), high_cutoff = as.numeric(-mds_high_cutoff),lower_cutoff=as.numeric(-mds_low_cutoff),
                  chr=as.character(clustered_chr))
    mds_clustering <- rbind(mds_clustering, tmp)
  }else{
    tmp <- tibble(mds_coord = as.character(mds_chosen),direction = as.character("neg"),clust_pvalue = as.numeric(NA), outliers=as.numeric(nrow(neg_windows)),
                  n1_outliers=as.numeric(NA), high_cutoff = as.numeric(NA),lower_cutoff=as.numeric(NA),
                  chr=as.numeric(NA))
    mds_clustering <- rbind(mds_clustering, tmp)
  }
}


dim(mds_clustering) 
head(mds_clustering)

mds_clustering$clust_pvalue
mds_clustering$outliers


# Keeps only the MDS coordinates where clustering was statistically significant (p < 0.01).
clust_pvalue_threshold <- 0.05 
mds_clustering %>% filter(clust_pvalue < clust_pvalue_threshold) -> sig_mds_clusters
sig_mds_clusters


# -- Step 9 and 10 in one big loop ---

#目的
# 1. check: current_windows are real outlier clusters — i.e., not isolated single points
# 2. run PCA on SNPs within the outlier region

# -- Step 9: Generate datasets of outlier windows and genotypes
outlier_windows <- tibble(chrom=character(), start=numeric(), end=numeric(), mid=numeric(), 
                          the_mds=numeric(), mds_coord=character(), outlier=character(), 
                          n=numeric())

cluster_genotypes <- tibble(mds_coord=character(), name=character(), PC1=numeric(), genotype=character())

i <- 1

for (i in 1:nrow(sig_mds_clusters)){

  coord <- pull(sig_mds_clusters[i,1])
  direction <- pull(sig_mds_clusters[i,2])
  high_cutoff <- pull(sig_mds_clusters[i,6])
  low_cutoff <- pull(sig_mds_clusters[i,7]) # CUT OFF APPLIED
  cluster_chr <- pull(sig_mds_clusters[i,8])
  coord_direction <- paste(coord, "-", direction, sep="")
  print(paste("Testing",coord_direction))
  
  if (direction == "pos"){
    current_windows <- win.regions %>%
      mutate(the_mds= .data[[coord]]) %>% 
      mutate(outlier=case_when((the_mds > low_cutoff) & (chrom == cluster_chr) ~ "Outlier", 
                               TRUE ~ "Non-outlier")) %>%
      filter(outlier != "Non-outlier") %>%
      select(chrom,start,end,mid,the_mds,outlier,n) %>%
      mutate(mds_coord=coord_direction) %>%
      mutate(ahead_n=n-lag(n), behind_n=abs(n-lead(n))) %>%
      mutate(min_dist=pmin(ahead_n, behind_n, na.rm=T)) %>%
      filter(min_dist <= max_distance_between_outliers) %>% # Removes isolated outliers. Keeps only outlier windows that are near at least one other outlier (e.g., part of a cluster).
      select(-ahead_n,-behind_n,-min_dist)
    windows <- current_windows %>% pull(n)
    outlier_windows <- rbind(outlier_windows, current_windows)
  }else{
    current_windows <- win.regions %>%
      mutate(the_mds= .data[[coord]]) %>% 
      mutate(outlier=case_when((the_mds < low_cutoff) & (chrom == cluster_chr) ~ "Outlier", 
                               TRUE ~ "Non-outlier")) %>%
      filter(outlier != "Non-outlier") %>%
      select(chrom,start,end,mid,the_mds,outlier,n) %>%
      mutate(mds_coord=coord_direction) %>%
      mutate(ahead_n=n-lag(n), behind_n=abs(n-lead(n))) %>%
      mutate(min_dist=pmin(ahead_n, behind_n, na.rm=T)) %>%
      filter(min_dist <= max_distance_between_outliers) %>%
      select(-ahead_n,-behind_n,-min_dist)
    windows <- current_windows %>% pull(n)
    outlier_windows <- rbind(outlier_windows, current_windows)
  } 


# -- Step 10: Determine 0/1/2 cluster on PC1 within outlier window

# 1. PCA from a custom SNP windowing function
samples <- attr(win.fn.snp, "samples")
out <- cov_pca(win.fn.snp(windows), k=2)
 
# 2. Reshape PCA output into a matrix
matrix.out <- t(matrix(out[4:length(out)], ncol=length(samples), byrow=T))

out <- as_tibble(cbind(samples, matrix.out)) %>% 
    rename(name=samples, PC1=V2, PC2=V3) %>% 
    mutate(PC1=as.double(PC1), PC2=as.double(PC2))


# 3. Tries to cluster samples along PC1 into 3 clusters
try_3_clusters <- try(kmeans(matrix.out[,1], 3, centers=c(min(matrix.out[,1]), (min(matrix.out[,1])+max(matrix.out[,1]))/2, max(matrix.out[,1]))))
  
if("try-error" %in% class(try_3_clusters)){
    kmeans_cluster <-kmeans(matrix.out[,1], 2, centers=c(min(matrix.out[,1]), max(matrix.out[,1])))
  }else{
    kmeans_cluster <- kmeans(matrix.out[,1], 3, centers=c(min(matrix.out[,1]), (min(matrix.out[,1])+max(matrix.out[,1]))/2, max(matrix.out[,1])))
  }

# 4. Cluster labels are stored (and converted to 0/1/2 instead of 1/2/3).
  out$cluster <- kmeans_cluster$cluster - 1
  out$cluster <- as.character(out$cluster)
  out$mds_coord <- paste(coord, direction, sep="-")
  genotype.out <- out %>% select(mds_coord,name,PC1,cluster) %>% rename(genotype=cluster)
  cluster_genotypes <- rbind(cluster_genotypes, genotype.out)

}

# If error occur "... XXX File name too long ..." it is becuase windows == 0 and does not run any downsteam PCA step with cov_pca() function
# Do not worry, it is not fatal or an error to the result, the remaining outlier windows is still printed and ran correctly through this big loop


#  -- Output 1: cluster_genotypes
dim(cluster_genotypes)
# 264 rows (from 88 samples * 3 outlier window calculations) 
# all outlier window ran successfully

#  -- Output 2: Outlier windows 
dim(outlier_windows)
# 31  8

outlier_windows



# -- Step 11: Counting
mds_counts <- tibble(mds_coord=character(), n_outliers=numeric())

for (mds in unique(outlier_windows$mds_coord)){
  count_outliers <- outlier_windows %>%
    filter(mds_coord == mds) %>% nrow()
  tmp_tibble <- tibble(mds_coord=as.character(mds), n_outliers=as.numeric(count_outliers))
  mds_counts <- rbind(mds_counts, tmp_tibble)
}

mds_counts
  mds_coord n_outliers
  <chr>          <dbl>
1 mds09-neg         12
2 mds12-neg         11
3 mds13-neg          8


# -- Step 12: Collapsing correlated MDS 
correlated_mds <- tibble(mds1=character(), mds2=character(), correlation=numeric())
total_mds_coords <- unique(outlier_windows$mds_coord)

for (i in 1:(length(total_mds_coords)-1)){
  for (j in (i+1):length(total_mds_coords)){

    mds1 <- total_mds_coords[[i]]
    mds2 <- total_mds_coords[[j]]
    chr1 <- outlier_windows %>% filter(mds_coord == mds1) %>% select(chrom) %>% unique() %>% pull()
    chr2 <- outlier_windows %>% filter(mds_coord == mds2) %>% select(chrom) %>% unique() %>% pull()
    if (chr1 != chr2){next;}
    cluster_genotypes %>% mutate(mds_coord=gsub("_", "-", mds_coord)) %>% filter(mds_coord == mds1 | mds_coord == mds2) %>%
      select(-PC1) %>%
      spread(mds_coord, genotype) %>% select(-name) -> tmp
    x <- tmp %>% pull(1) %>% as.numeric()
    y <- tmp %>% pull(2) %>% as.numeric()
    test_result <- cor.test(x, y, na.rm=T)
    tmp_tibble <- tibble(mds1=as.character(mds1), mds2=as.character(mds2), correlation=as.numeric(abs(test_result$estimate)))
    correlated_mds <- rbind(correlated_mds, tmp_tibble) # This step did not bind? very odd

  }
}

dim(correlated_mds)
correlated_mds
  mds1      mds2      correlation
  <chr>     <chr>           <dbl>
1 mds09-neg mds12-neg       0.957
# both close by on CHR11, should be collapsed



# -- This loop de-duplicates highly correlated MDS windows (on the same Chromosome)
# -- if highly correlated keeping the MDS windows with more outliers (presumably more informative) found in total_mds_coords

for (i in 1:nrow(correlated_mds)){
  if (correlated_mds[i,3] >= min_cor){
    count1 <- mds_counts %>% filter(mds_coord == as.character(correlated_mds[i,1])) %>% pull(n_outliers)
    count2 <- mds_counts %>% filter(mds_coord == as.character(correlated_mds[i,2])) %>% pull(n_outliers)
    if (count1 < count2){
      total_mds_coords[which(total_mds_coords != as.character(correlated_mds[i,1]))] -> total_mds_coords
    }else{
      total_mds_coords[which(total_mds_coords != as.character(correlated_mds[i,2]))] -> total_mds_coords
    }
  }
}


total_mds_coords





# -- Step 13: Plot outliers (putative inversions)

cluster_genotypes <- tibble(mds_coord=character(), name=character(), 
                            PC1=numeric(), genotype=character())

mds_info <- tibble(mds_coord=character(), mds=character(), chromosome=character(), 
                   start=numeric(), end=numeric(), n_outliers=numeric(),
                   PC1_perc=numeric(), PC2_perc=numeric(), 
                   betweenSS_perc=numeric())




# This i is for nrow(sig_mds_clusters)

i <- 1 
i <- 3

#下一步才会进一步筛选是不是被corr已经筛选掉，如果是的话就不会作图

for (i in 1:nrow(sig_mds_clusters)){

# -- Prepare dataframe for plotting

  coord <- pull(sig_mds_clusters[i,1])
  direction <- pull(sig_mds_clusters[i,2])
  if(! paste(coord,"-", direction, sep="") %in% total_mds_coords){
    next;
  }
  high_cutoff <- pull(sig_mds_clusters[i,6])
  low_cutoff <- pull(sig_mds_clusters[i,7])
  cluster_chr <- pull(sig_mds_clusters[i,8])
  coord_direction <- paste(coord, "-", direction, sep="")
  
  print(coord_direction)


  if (direction == "pos"){
    current_windows <- win.regions %>%
      mutate(the_mds=.data[[coord]] ) %>% 
      mutate(outlier=case_when((the_mds > low_cutoff) & (chrom == cluster_chr) ~ "Outlier", 
                               TRUE ~ "Non-outlier")) %>%
      filter(outlier != "Non-outlier") %>%
      select(chrom,start,end,mid,the_mds,outlier,n) %>%
      mutate(mds_coord=coord_direction) %>%
      mutate(ahead_n=n-lag(n), behind_n=abs(n-lead(n))) %>%
      mutate(min_dist=pmin(ahead_n, behind_n, na.rm=T)) %>%
      filter(min_dist < max_distance_between_outliers ) %>%
      select(-ahead_n,-behind_n,-min_dist)
    windows <- current_windows %>% pull(n)
  }else{
    current_windows <- win.regions %>%
      mutate(the_mds=.data[[coord]]) %>% 
      mutate(outlier=case_when((the_mds < low_cutoff) & (chrom == cluster_chr) ~ "Outlier", 
                               TRUE ~ "Non-outlier")) %>%
      filter(outlier != "Non-outlier") %>%
      select(chrom,start,end,mid,the_mds,outlier,n) %>%
      mutate(mds_coord=coord_direction) %>%
      mutate(ahead_n=n-lag(n), behind_n=abs(n-lead(n))) %>%
      mutate(min_dist=pmin(ahead_n, behind_n, na.rm=T)) %>%
      filter(min_dist < max_distance_between_outliers) %>%
      select(-ahead_n,-behind_n,-min_dist)
    windows <- current_windows %>% pull(n)
  }

# Outputs "current_windows" and "windows"






# -- Plot # 1: 显著的mds在每个染色体上的 outlier window 区域
  chromosome <- current_windows %>% head(1) %>% pull(chrom)
  start <- current_windows %>% summarize(s=min(start)) %>% pull()
  end <- current_windows %>% summarize(e=max(end)) %>% pull()
  
  genome_plot <- win.regions %>% mutate(the_mds= .data[[coord]]) %>% 
    mutate(chrom_num=substr(win.regions$chrom, 8,9)) %>% 
    mutate(outlier=case_when(n %in% current_windows$n ~ "Outlier", 
                             TRUE ~ "Non-outlier")) %>%
    ggplot(., aes(x=mid/1000000, y=the_mds, color=outlier)) + geom_point() + theme_bw() +
    facet_wrap(~chrom_num, scales="free_x", nrow=1) +
    scale_color_manual(values=c("grey40","#E41A1C")) +
    xlab("Mbp") + ylab(toupper(coord)) +
    theme(legend.title=element_blank()) + labs(tag = "(a)")

png(filename = paste0("Cluster_B_lostruct_win50_manplot_",coord_direction,"_2026Apr17.png"),  width = 1500, height = 800, res = 200)
print(genome_plot)
dev.off()



# -- Plot # 2: PCA 三个显著分区
  out <- cov_pca(win.fn.snp(windows), k=2)
  PC1_perc <- out[2]/out[1]
  PC2_perc <- out[3]/out[1]
  matrix.out <- t(matrix(out[4:length(out)], ncol=length(samples), byrow=T)) 
  out <- as_tibble(cbind(samples, matrix.out)) %>% 
    rename(name=samples, PC1=V2, PC2=V3) %>% 
    mutate(PC1=as.double(PC1), PC2=as.double(PC2))
  #head(out)

  try_3_clusters <-try(kmeans(matrix.out[,1], 3, centers=c(min(matrix.out[,1]), (min(matrix.out[,1])+max(matrix.out[,1]))/2, max(matrix.out[,1]))))
  
  if("try-error" %in% class(try_3_clusters)){
    kmeans_cluster <-kmeans(matrix.out[,1], 2, centers=c(min(matrix.out[,1]), max(matrix.out[,1])))
  }else{
    kmeans_cluster <- kmeans(matrix.out[,1], 3, centers=c(min(matrix.out[,1]), (min(matrix.out[,1])+max(matrix.out[,1]))/2, max(matrix.out[,1])))
  }


  out$cluster <- kmeans_cluster$cluster - 1
  out$cluster <- as.character(out$cluster)
  out$mds_coord <- paste(coord, direction, sep="_")
  betweenSS_perc <- kmeans_cluster$betweenss / kmeans_cluster$totss


  pca_plot <- out %>%
    ggplot(., aes(x=PC1, y=PC2, col=cluster)) + geom_point() + theme_bw() +
    scale_color_manual(name="Cluster", values=c("red","purple","blue")) +
    xlab(paste("PC",1," (",round(100*PC1_perc,2),"% PVE)",sep="")) +
    ylab(paste("PC",2," (",round(100*PC2_perc,2),"% PVE)",sep="")) +
    theme(plot.margin = unit(c(0.3,0.3,0.3,0.3), "inches")) + labs(tag = "(b)")

png(filename = paste0("Cluster_B_lostruct_win50_PCAplot_",coord_direction,"_2026Apr17.png"),  width = 1000, height = 700, res = 200)
print(pca_plot)
dev.off()
   



# -- Plot # 3: Heterzygosity boxplot

  win.fn.snp(windows) %>% as_tibble() -> snps
  colnames(snps) <- samples

  windows 
  # 851 852 853 854 856 857 858 859

  dim(snps)
  snps[1:6,1:5]

  snps %>% gather("name","genotype",1:ncol(snps)) %>% 
    filter(!is.na(genotype)) %>%
    group_by(name, genotype) %>%
    summarize(count=n()) %>%
    spread(genotype, count) %>%
    summarize(het=`1`/(`0` + `1` + `2`)) -> heterozygosity

  head(heterozygosity)

  het_plot <- inner_join(out, heterozygosity) %>% 
    ggplot(., aes(x=as.character(cluster), y=het, fill=as.character(cluster))) + geom_boxplot() + theme_bw() +
    scale_fill_manual(name="Cluster", values=c("red","purple","blue")) +
    xlab("Cluster") + ylab("Heterozygosity") +
    theme(plot.margin = unit(c(0.3,0.3,0.3,0.3), "inches")) + labs(tag = "(c)")
  
png(filename = paste0("Cluster_B_lostruct_win50_HET_",coord_direction,"_2026Apr17.png"),  width = 1000, height = 700, res = 200)
print(het_plot)
dev.off()




# -------------------
#
# Step 3 : Save output
#
# -------------------

 # Output the information
  genotype.out <- out %>% select(mds_coord,name,PC1,cluster) %>% rename(genotype=cluster) 
  cluster_genotypes <- rbind(cluster_genotypes, genotype.out)
  tmp_info <- tibble(mds_coord=as.character(coord_direction), mds=as.character(toupper(coord)), chromosome=as.character(chromosome), 
                     start=as.numeric(start), end=as.numeric(end), n_outliers=as.numeric(length(windows)),
                     PC1_perc=as.numeric(round(PC1_perc*100,2)), PC2_perc=as.numeric(round(PC2_perc*100,2)),
                     betweenSS_perc=as.numeric(round(betweenSS_perc,4)))
  mds_info <- rbind(mds_info, tmp_info)

    } # END OF THE BIG PLOTTING LOOP OF FINAL all WINDOWS 




library(readr)

write_tsv(cluster_genotypes, "Cluster_B_INV_genotypes_20250417.tsv")
write_tsv(mds_info, "Cluster_B_mds_info_20250417.tsv")
