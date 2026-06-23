#####################################
# Project: debilis_praecox_project
#
# Code 4: GBS on the Global Set
#
# by: Yue Yu
#
#####################################

# Step 1: Trimmomatic       fastq
# Step 2: Multi-QC          fastq
# Step 3.1: bwa-mem2:       fastq -> BAM
# Step 3.2: HaplotypeCaller:  BAM -> g.vcf
# Step 3.3: GenomicDBImport:         g.vcf -> database
# Step 3.4: GenotypeGVCF:                     database -> VCF
# Step 4.1: VCFtools:                                     VCF -> raw SNPs 
# Step 4.2: SNP filtering:										 raw SNPs -> filtered SNPs

# Step 5: ADMIXTURE
# Step 6: PCA



# The following code includes 

# --------------------------------
#   PRA called on PRA-H2 ref (DONE 2025 April)
# --------------------------------
# Previously done 2025 April
# BAM saved in: /home/yueyu/scratch/GBS/BAM_PRA
# GVCF saved in: /home/yueyu/scratch/GBS/G_VCF_GATK_PRA


# --------------------------------
#   DEB_Texas called on PRA-H2 ref (DONE 2025 June)
# --------------------------------
# Previously done 2025 June
# BAM saved in: /home/yueyu/scratch/PET_forTexas/BAM_DEB_callon_PRAHap2
# GVCF saved in: /home/yueyu/scratch/PET_forTexas/G_VCF_GATK_DEB_callon_PRAHap2


# --------------------------------
#   DEB_Florida call on PRA-H2 ref (IN CODE BELOW)
# --------------------------------
cd /home/yueyu/scratch/ALL_GBS_call_on_PRA/BAM_DEB_FLORIDA_callon_PRAHap2

nano run_BWA_PRA_2026Mar29.sh

# -------------- run_BWA_PRA_2026Mar29.sh --------------------
#!/bin/bash
#SBATCH --account=rrg-rieseber-ac
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G  
#SBATCH --time=04:00:00
#SBATCH --array=1-204  # 204 FLORIDA DEB samples

module load StdEnv/2023
module load bwa-mem2/2.2.1
module load gcc/12.3
module load samtools/1.20

cd /home/yueyu/scratch/GBS/Trimmed/subset_DEB_Trimmed

# Define variables
REF="/home/yueyu/scratch/PRA_FASTA/praecox2.fasta"  # PRAECOX

# Extract TEXAS name from all files
i=$(ls *paired_R1.fastq.gz | grep -E '^DD|^DT|^DV'| head -n $SLURM_ARRAY_TASK_ID | tail -n 1)   # TEXAS NAMES FOR 204 SAMPLES
SAMPLE=$(echo $i | cut -d "_" -f 1-2)

R1="${SAMPLE}_paired_R1.fastq.gz"
R2="${SAMPLE}_paired_R2.fastq.gz"

# OutputBAM
OUTPUT_BAM="/home/yueyu/scratch/ALL_GBS_call_on_PRA/BAM_DEB_FLORIDA_callon_PRAHap2/${SAMPLE}.callonPRAHap2.sort.bam"

# BWA-MEM2 & samtools sort
bwa-mem2 mem -t 8 -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:illumina" \
    "$REF" "$R1" "$R2" | \
samtools sort -@ 8 -m 7G -o "$OUTPUT_BAM"

# -- Note:
# -R: add read group (RG) information to the output BAM/SAM file. 
# It is necessary for GATK to run later on.

# -------------- run_BWA_PRA_2026Mar29.sh (END) --------------------
sbatch run_BWA_PRA_2026Mar29.sh


#================================
#  Part 3: Index BAM
#================================
module load samtools/1.20

# -- Index bam file --> bam.bai
cd /home/yueyu/scratch/ALL_GBS_call_on_PRA/BAM_DEB_FLORIDA_callon_PRAHap2

for bam in *.sort.bam; do
    samtools index "$bam"
done


#================================
#  Part 5: Calculate Coverage
#================================
module load StdEnv/2023
module load nextgenmap/0.5.5
module load samtools/1.20

cd /home/yueyu/scratch/ALL_GBS_call_on_PRA/BAM_DEB_FLORIDA_callon_PRAHap2

DEB_COVERAGE="DEB_Florida_GBS_sample_callon_PRAHap2_coverage_results.txt"

for bam_file in *.sort.bam; do
  sample=$(echo $bam_file | cut -d "." -f 1) # Extract sample name by removing the .sort.bam extension
  echo ${sample}

  # Calculate coverage
  samtools depth "${bam_file}" > "${sample}_coverage.txt"
  average_coverage=$(awk '{sum+=$3} END {print sum/NR}' "${sample}_coverage.txt" | bc -l)
  
  # Output the result to the file
  echo "Sample: ${sample}, Overall Coverage: ${average_coverage}" >> "${DEB_COVERAGE}"
done

# DONE - MOVED TO LOCAL COMPUTER TOO
cd /Users/yueyu/Desktop/GBS/Coverage_2026March
scp yueyu@fir.computecanada.ca:"/home/yueyu/scratch/ALL_GBS_call_on_PRA/BAM_DEB_FLORIDA_callon_PRAHap2/DEB_Florida_GBS_sample_callon_PRAHap2_coverage_results.txt" .

# Remove individual coverage files
# rm *coverage.txt


#----------------  run_DEB_HapCaller_2026March29.sh (start)  ----------------
#!/bin/bash
#SBATCH --account=rrg-rieseber-ac
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30G
#SBATCH --array=1-204

module load StdEnv/2020
module load gatk/4.2.4.0

cd /home/yueyu/scratch/ALL_GBS_call_on_PRA/BAM_DEB_FLORIDA_callon_PRAHap2

# Define variables
REF="/home/yueyu/scratch/PRA_FASTA/praecox2.fasta"  # !!!!  PRAECOX !!!! CAUTION!!

i=$(ls *sort.bam | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
SAMPLE=$(echo $i | cut -d "." -f 1) 

gatk --java-options "-Xmx28G" HaplotypeCaller \
	 -R ${REF} \
	 -I ${SAMPLE}.callonPRAHap2.sort.bam \
	 -O /home/yueyu/scratch/ALL_GBS_call_on_PRA/G_VCF_GATK_DEB_Florida_callon_PRAHap2/${SAMPLE}.callonPRAHap2.g.vcf \
	 -ERC GVCF \
	 --max-alternate-alleles 3 \
	 --pcr-indel-model AGGRESSIVE \
	 -G StandardAnnotation -G AS_StandardAnnotation 

#----------------  run_DEB_HapCaller_2026March29.sh (end)  ----------------


# Pool all samples and call SNPs for the global VCF set
# ---------  PREACOX  (130 samples) ------------
# ---------  DEBILIS Texas (117 samples)-------------
# ---------  DEBILIS Florida  (204 samples) -------------

cd /home/yueyu/scratch/ALL_GBS_call_on_PRA/G_VCF_ALL_merged

for i in *.g.vcf; do
    name=$(echo "$i" | cut -d "." -f 1)
    echo -e "$name\t$i" >> sample_name_451_all.txt
done

head sample_name_451_all.txt


# -- Import all samples into DBI per CHROMOSOME 

cd /home/yueyu/scratch/ALL_GBS_call_on_PRA/G_VCF_ALL_merged

for i in $(seq -w 01 17)
do
  cat <<EOL > PRA_Chr${i}_makeDB.sh
#!/bin/bash
#SBATCH --account=rrg-rieseber-ac
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=10G

module load StdEnv/2020
module load gatk/4.2.4.0

cd /home/yueyu/scratch/ALL_GBS_call_on_PRA/G_VCF_ALL_merged

gatk --java-options "-Xmx48g" GenomicsDBImport \\
    --genomicsdb-workspace-path PRA_CHROM$i \\
    --batch-size 50 \\
    --sample-name-map sample_name_451_all.txt \\
    --reader-threads 3 \\
    -L PRA_chr$i
EOL

  chmod +x PRA_Chr${i}_makeDB.sh
done


for i in $(seq -w 01 17)
do
  cat <<EOL > PRA_Chr${i}_genotypeVCF.sh
#!/bin/bash
#SBATCH --account=rrg-rieseber-ac
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G

module load StdEnv/2020
module load gatk/4.2.4.0

cd /home/yueyu/scratch/ALL_GBS_call_on_PRA/G_VCF_ALL_merged

gatk --java-options "-Xmx90g" GenotypeGVCFs \\
     -R /home/yueyu/scratch/PRA_FASTA/praecox2.fasta \\
     -V gendb://PRA_CHROM$i \\
     -O /home/yueyu/scratch/ALL_GBS_call_on_PRA/VCF/PRA_Chr$i.vcf.gz
EOL

  chmod +x PRA_Chr${i}_genotypeVCF.sh

done


#Make final run script
for i in {01..17};
do 
    echo "sbatch PRA_Chr${i}_genotypeVCF.sh" >> run_PRA_GT_byCHR.sh
done


nano run_PRA_vcf_merge.sh

#----------  vcf_merge.sh  ---------- 
#!/bin/bash
#SBATCH --account=rrg-rieseber-ac
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G

module load StdEnv/2020
module load gcc/9.3.0
module load bcftools/1.16
module load tabix/0.2.6


cd /home/yueyu/scratch/ALL_GBS_call_on_PRA/VCF

bcftools concat \
  PRA_Chr01.vcf.gz \
  PRA_Chr02.vcf.gz \
  PRA_Chr03.vcf.gz \
  PRA_Chr04.vcf.gz \
  PRA_Chr05.vcf.gz \
  PRA_Chr06.vcf.gz \
  PRA_Chr07.vcf.gz \
  PRA_Chr08.vcf.gz \
  PRA_Chr09.vcf.gz \
  PRA_Chr10.vcf.gz \
  PRA_Chr11.vcf.gz \
  PRA_Chr12.vcf.gz \
  PRA_Chr13.vcf.gz \
  PRA_Chr14.vcf.gz \
  PRA_Chr15.vcf.gz \
  PRA_Chr16.vcf.gz \
  PRA_Chr17.vcf.gz \
  -O z > /home/yueyu/scratch/ALL_GBS_call_on_PRA/VCF/PRA_451_sample_raw_v2026Mar31.vcf.gz

tabix -p vcf PRA_451_sample_raw_v2026Mar31.vcf.gz

#----------  vcf_merge.sh  (END) ---------- 




# --------------------------------
#
#
#   Step 4.2: SNP filtering
#
#
# --------------------------------

# SNP FILTERING STEPS BELOW

#--------------------------
#  Step 2: Extract SNPs
#--------------------------

cd /home/yueyu/scratch/ALL_GBS_call_on_PRA/Filtered_VCF

nano Filter_Step1.sh
# ----------- Filter_Step1.sh ---------
#!/bin/bash
#SBATCH --account=rrg-rieseber-ac 
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=2:00:00

cd /home/yueyu/scratch/ALL_GBS_call_on_PRA/Filtered_VCF

module load StdEnv/2023
module load gatk/4.6.1.0

unset JAVA_TOOL_OPTIONS

#SNP
gatk --java-options "-Xmx8g" SelectVariants \
  -V /home/yueyu/scratch/ALL_GBS_call_on_PRA/VCF/PRA_451_sample_raw_v2026Mar31.vcf.gz \
  --select-type-to-include SNP \
  --exclude-filtered \
  -O SNP_FILTERED.vcf.gz
# ----------- Filter_Step1.sh (END)---------

# Submitted batch job 30828893 - DONE

module load StdEnv/2020
module load gcc/9.3.0
module load bcftools/1.16

bcftools view -H SNP_FILTERED.vcf.gz | wc -l 
# 4,611,861

#--------------------------
#  Step 3: Plot SNP quality (DONE)
#--------------------------
module load StdEnv/2023  
module load gcc/12.3
module load bcftools/1.19

cd /home/yueyu/scratch/ALL_GBS_call_on_PRA/Filtered_VCF

#Extract info and plot in R 
bcftools query SNP_FILTERED.vcf.gz -f '%ExcessHet\t%FS\t%SOR\t%MQRankSum\t%ReadPosRankSum\t%QD\t%MQ\t%DP\n' > SNP_quality_check/SNP_HET.FS.SOR.MQRS.RPRS.QD.MQ.DP.txt

cd /home/yueyu/scratch/ALL_GBS_call_on_PRA/Filtered_VCF/SNP_quality_check
wc -l SNP_HET.FS.SOR.MQRS.RPRS.QD.MQ.DP.txt
# 4,611,861 (4 million SNPs)


module load StdEnv/2023
module load r/4.4.0

R
# ------------- in R -------------
library(data.table)

# Load DEB SNP table
snps <- fread("SNP_HET.FS.SOR.MQRS.RPRS.QD.MQ.DP.txt", na.strings = c(".", "NA"), colClasses = "numeric", header = F)

str(snps)
colnames(snps) <- c("ExcessHet","FS","SOR","MQRankSum","ReadPosRankSum","QD","MQ","DP")

# !!! IMPORTANT!! change saved PNG names

# -- ExcessHet
png(filename = "ExcessHet_2026March31.png")
dExcessHet <- density(snps$ExcessHet,na.rm=T)
plot(dExcessHet,main="ExcessHet distribution for PRA", xlab="ExcessHet")
dev.off()

# -- FS
# or can plot FS with a log base 10 scale
png(filename = "FS_2026March31.png")
dFS <- density(snps$FS,na.rm=T)
plot(dFS,main="FS distribution for PRA", xlab="FS")
dev.off()

# -- SOR
png(filename = "SOR_2026March31.png")
dSOR <- density(snps$SOR,na.rm=T)
plot(dSOR,main="SOR distribution for PRA", xlab="SOR")
dev.off()

# -- MQRankSum
png(filename = "MQRankSum_2026March31.png")
dMQRankSum <- density(snps$MQRankSum,na.rm=T)
plot(dMQRankSum,main="MQRankSum distribution for PRA", xlab="MQRankSum")
dev.off()

# -- ReadPosRankSum
png(filename = "ReadPosRankSum_2026March31.png")
dReadPosRankSum <- density(snps$ReadPosRankSum,na.rm=T)
plot(dReadPosRankSum,main="ReadPosRankSum distribution for PRA", xlab="ReadPosRankSum")
dev.off()

# -- QD
png(filename = "QD_2026March31.png")
dQD <- density(snps$QD,na.rm=T)
plot(dQD,main="QD distribution for PRA", xlab="QD")
dev.off()

# -- MQ
png(filename = "MQ_2026March31.png")
dMQ <- density(snps$MQ,na.rm=T)
plot(dMQ,main="MQ distribution for PRA", xlab="MQ")
dev.off()

# -- DP
png(filename = "DP_2026March31.png")
dDP <- density(snps$DP,na.rm=T)
plot(dDP,main="DP distribution for PRA", xlab="DP")
dev.off()



#--------------------------
#  Step 3: Hard filter
#--------------------------

cd /home/yueyu/scratch/ALL_GBS_call_on_PRA/Filtered_VCF

nano Filter_Step2.sh
# ----------- Filter_Step2.sh ---------
#!/bin/bash
#SBATCH --account=rrg-rieseber-ac
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=1:00:00

cd /home/yueyu/scratch/ALL_GBS_call_on_PRA/Filtered_VCF

module load StdEnv/2023
module load gatk/4.6.1.0

unset JAVA_TOOL_OPTIONS

gatk --java-options "-Xmx8g" VariantFiltration \
    -V SNP_FILTERED.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSumNeg12.5" \
    -filter "ExcessHet > 54.69" --filter-name "ExcessHet5469" \
    -O SNP_INFO_MARKED.vcf.gz


# ADD CODE TO EXCLUDE THE FILTERED
gatk --java-options "-Xmx8g" SelectVariants \
  -V SNP_INFO_MARKED.vcf.gz \
  --exclude-filtered \
  -O SNP_INFO_FILTERED.vcf.gz

# ----------- Filter_Step2.sh (END)---------

# Submitted batch job 30832542
# only took 7 min



moduel load vcftools
vcftools --gzvcf SNP_INFO_FILTERED.vcf.gz --missing-site --out SNP_INFO_FILTERED_site_missing

awk 'NR>1 {
  if($6 <= 0.1) bin1++ ; 
  else if($6 <= 0.2) bin2++ ; 
  else if($6 <= 0.3) bin3++ ; 
  else if($6 <= 0.4) bin4++ ; 
  else if($6 <= 0.5) bin5++ ; 
  else if($6 <= 0.6) bin6++ ; 
  else if($6 <= 0.7) bin7++ ; 
  else if($6 <= 0.8) bin8++ ; 
  else if($6 <= 0.9) bin9++ ; 
  else bin10++
} END {
  print "0–0.1:",bin1, \
        "\n0.1–0.2:",bin2, \
        "\n0.2–0.3:",bin3, \
        "\n0.3–0.4:",bin4, \
        "\n0.4–0.5:",bin5, \
        "\n0.5–0.6:",bin6, \
        "\n0.6–0.7:",bin7, \
        "\n0.7–0.8:",bin8, \
        "\n0.8–0.9:",bin9, \
        "\n0.9–1.0:",bin10
}' SNP_INFO_FILTERED_site_missing.lmiss

# --  Result -- 
0–0.1: 3724957
0.1–0.2: 1558
0.2–0.3: 527
0.3–0.4: 211
0.4–0.5: 130
0.5–0.6: 62
0.6–0.7: 23
0.7–0.8: 9
0.8–0.9:
0.9–1.0:






# --------------------
# Step 6: BCFTOOLS on SAMPLE field for SAMPLE-DP and SAMPLE-GQ (Zhe code, very good!)
# --------------------

cd /home/yueyu/scratch/ALL_GBS_call_on_PRA/Filtered_VCF

#!! IMPORTANT NOTE !!
# AFTER MUCH STRUGGLES 
# +++++++++++++++ SEPERATOR AS | NOT || , not regular expression, specific rule in BCFTOOLS for this step ++++++++++

nano Filter_Step3.sh
# ----------- Filter_Step3.sh ---------
#!/bin/bash
#SBATCH --account=rrg-rieseber-ac
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=2:00:00

cd /home/yueyu/scratch/ALL_GBS_call_on_PRA/Filtered_VCF

module load StdEnv/2020
module load gcc/9.3.0
module load bcftools/1.16
module load tabix/0.2.6

bcftools filter \
    -e '(GQ<30 | FORMAT/DP<10)' \
    SNP_INFO_FILTERED.vcf.gz \
    -S . \
    -s "PASS" \
    -Oz -o SNP_INFO_GENO_FILTERED.vcf.gz

tabix -p vcf SNP_INFO_GENO_FILTERED.vcf.gz

# ----------- Filter_Step3.sh (END)---------

# Submitted batch job 30834381 - DONE




# --------------------
# Step 7: FILTER MAX MISSINGNESS (SNP level)
# --------------------

cd /home/yueyu/scratch/ALL_GBS_call_on_PRA/Filtered_VCF

nano Filter_Step4.sh
# ----------- Filter_Step4.sh ---------
#!/bin/bash
#SBATCH --account=rrg-rieseber-ac
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=1:00:00

cd /home/yueyu/scratch/ALL_GBS_call_on_PRA/Filtered_VCF

module load StdEnv/2023
module load gatk/4.6.1.0

unset JAVA_TOOL_OPTIONS

gatk --java-options "-Xmx18g" SelectVariants \
            -V SNP_INFO_GENO_FILTERED.vcf.gz \
            --remove-unused-alternates \
            --restrict-alleles-to BIALLELIC \
            --max-nocall-fraction 0.3 \
            -O SNP_INFO_GENO_BI_NOCALL03_FILTERED.vcf.gz

# ----------- Filter_Step4.sh (END)---------

# Submitted batch job 30843542 -- DONE 

module load StdEnv/2020
module load gcc/9.3.0
module load bcftools/1.16

bcftools view -H SNP_INFO_GENO_BI_NOCALL03_FILTERED.vcf.gz | wc -l
# 242,812



# --------------------
# Step 8: REMOVE SAMPLES THAT ARE OF LOW COVERAGE -- YES, do it after filter for missing rate (of the sites)
# -------------------- 

module load StdEnv/2020 
module load vcftools/0.1.16

cd /home/yueyu/scratch/ALL_GBS_call_on_PRA/Filtered_VCF

vcftools --gzvcf SNP_INFO_GENO_BI_NOCALL03_FILTERED.vcf.gz --missing-indv --out sample_missing_ratio
awk '$5 > 0.3' sample_missing_ratio.imiss

# Just one sample that is low coverage
# INDV  N_DATA  N_GENOTYPES_FILTERED  N_MISS  F_MISS
# PH1_10  242812  0 120158  0.49486







# --------------------
# Step 8(ADDITIONAL): REMOVE SAMPLES FROM PREVUOUS KNOWN PCA RESULTS
# -------------------- 

# REDO 2026 May 14, FOR PUBLICATION, MATCH WITH THE OTHER SUBSETS


# Samples to remove (from previous PCA results)
SHOULD REMOVE (35 samples in total): 
Sample DD1_10,  DD4_3, PR6_7, DC2_8
low coverage: PH1_10

Population DC1, DC4 and DC7

nano samples_to_remove.txt
PH1_10
PR6_7
DC2_8
DD4_3
DD1_10


# -- Keep DC 1,4,and 7 (2026 March 31) --
# Reasons behind this:
# Based on previous ADMIXTURE + PCA results, these 3 pops (especially DC 1 and 7) seperate from the Texas and Floria cluster
# DC 4 on PCA is seperated from all, but in ADMIXTURE show cluster share with DS4


# -- remove samples
bcftools view -S ^samples_to_remove.txt -Oz -o SNP_INFO_GENO_BI_NOCALL03_FILTERED_removed5samples.vcf.gz SNP_INFO_GENO_BI_NOCALL03_FILTERED.vcf.gz
tabix -p vcf SNP_INFO_GENO_BI_NOCALL03_FILTERED_removed5samples.vcf.gz


# -- REMOVE DC 1,4,and 7 (2026 April 03) --
bcftools query -l SNP_INFO_GENO_BI_NOCALL03_FILTERED_removed5samples.vcf.gz \
  | grep -E '^(DC1|DC4|DC7)' > samples_to_remove_30.txt

bcftools view -S ^samples_to_remove_30.txt -Oz -o SNP_INFO_GENO_BI_NOCALL03_FILTERED_removed35samples.vcf.gz SNP_INFO_GENO_BI_NOCALL03_FILTERED_removed5samples.vcf.gz
tabix -p vcf SNP_INFO_GENO_BI_NOCALL03_FILTERED_removed35samples.vcf.gz

bcftools query -l SNP_INFO_GENO_BI_NOCALL03_FILTERED_removed35samples.vcf.gz | wc -l
# Removed the number of samples left 416

# ------------------------------------
# Step 9: RECODE INFO + AF filter
# -----------------------------------


nano Filter_Step5.sh
# ----------- Filter_Step5.sh ---------
#!/bin/bash
#SBATCH --account=rrg-rieseber-ac
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --time=2:00:00

cd /home/yueyu/scratch/ALL_GBS_call_on_PRA/Filtered_VCF

module load StdEnv/2020
module load gcc/9.3.0
module load bcftools/1.16
module load tabix/0.2.6

# not needed and recode information (INFO - AN, AC, AF, DP need to recode for downstream)
# Remove other INFO field that are not accurate anymore and can not be re-calculated after filtering (eg. require initial seq info)
# ^ means to keep these INFO fields but drop the rest

bcftools annotate -x "^INFO/AN,INFO/AC,INFO/AF,INFO/DP" SNP_INFO_GENO_BI_NOCALL03_FILTERED_removed35samples.vcf.gz -Ou | bcftools +fill-tags  -O z -o SNP_INFO_GENO_BI_NOCALL03_FILTERED_RECODED_CLEANED.vcf.gz -- -t AN,AC,AF
tabix -p vcf SNP_INFO_GENO_BI_NOCALL03_FILTERED_RECODED_CLEANED.vcf.gz

# -- AF filter
bcftools view \
  -i 'AF>=0.03 && AF<=0.97' \
  -O z \
  -o SNP_INFO_GENO_BI_NOCALL03_FILTERED_RECODED_CLEANED_AF003.vcf.gz \
  SNP_INFO_GENO_BI_NOCALL03_FILTERED_RECODED_CLEANED.vcf.gz

tabix -p vcf SNP_INFO_GENO_BI_NOCALL03_FILTERED_RECODED_CLEANED_AF003.vcf.gz

# ----------- Filter_Step5.sh (END) ---------












# --------------------------------
#
#
#   Step 5: ADMIXTURE
#
#
# --------------------------------

#####################################
#  Step 0: VCF LD filter
#####################################

VCF="/home/yueyu/scratch/ALL_GBS_call_on_PRA/Filtered_VCF/SNP_INFO_GENO_BI_NOCALL03_FILTERED_RECODED_CLEANED_AF003.vcf.gz" 

cd /scratch/yueyu/ALL_GBS_call_on_PRA/Admixture

#Default version: PLINK v1.90b6.21 64-bit (19 Oct 2020)
module load StdEnv/2020
module load plink/1.9b_6.21-x86_64

plink --vcf $VCF --keep-allele-order --make-bed --set-missing-var-ids @:# --allow-extra-chr --double-id --out GBS_all_bfile
#50511 variants and 416 people pass filters and QC.

# Step 1: LD prune
plink --bfile GBS_all_bfile \
      --indep-pairwise 10kb 50 0.2 \
      --allow-extra-chr \
      --out GBS_all_bfile_10Kb_r02

# Step 2: Apply LD pruning
plink --bfile GBS_all_bfile \
      --extract GBS_all_bfile_10Kb_r02.prune.in \
      --make-bed \
      --allow-extra-chr \
      --out GBS_all_LDpruned

# Step 3: Distance prune (reduced by keeping 1 SNP/300 bp)
plink --bfile GBS_all_LDpruned \
      --bp-space 300 \
      --make-bed \
      --allow-extra-chr \
      --recode vcf \
      --keep-allele-order \
      --out GBS_all_LDpruned_dist300
#9594 variants and 416 people pass filters and QC

module load StdEnv/2023
module load bcftools
module load tabix/0.2.6

bgzip GBS_all_LDpruned_dist300.vcf
tabix -p vcf GBS_all_LDpruned_dist300.vcf.gz

# rm *.log *.nosex *.bim *.bed *.fam *.prune* 

bcftools query -l GBS_all_LDpruned_dist300.vcf.gz > GBS_all_LD_pruned_samples.txt
wc -l GBS_all_LD_pruned_samples.txt
# 416 DEB_LD_pruned_samples.txt





#####################################
#  Step 1: VCF -> BED
#####################################

cd /scratch/yueyu/ALL_GBS_call_on_PRA/Admixture

module load nixpkgs/16.09 intel/2018.3
module load tabix/0.2.6

#Change name to integer: PRA_chr01 >> 01
zcat GBS_all_LDpruned_dist300.vcf.gz |\
  sed s/^PRA_chr//g |\
  bgzip > GBS_all_LDpruned_dist300_numericChr.vcf.gz

tabix GBS_all_LDpruned_dist300_numericChr.vcf.gz
 
module load StdEnv/2020
module load plink/1.9b_6.21-x86_64

#Then, run PLINK
plink --make-bed \
    --vcf GBS_all_LDpruned_dist300_numericChr.vcf.gz \
    --out GBS_all_LDpruned_dist300_numericChr \
    --set-missing-var-ids @:# \
    --double-id \
    --allow-extra-chr



#####################################
#  Step 2: Run ADMIXTURE
#####################################

cd /scratch/yueyu/ALL_GBS_call_on_PRA/Admixture

module load StdEnv/2020
module load admixture/1.3.0

for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15; \
do admixture --cv GBS_all_LDpruned_dist300_numericChr.bed $K |\
tee GBS_all.${K}.out; \
done

cat *out | grep CV

CV error (K=1): 0.61037
CV error (K=2): 0.40407
CV error (K=3): 0.33614
CV error (K=4): 0.30604
CV error (K=5): 0.29175
CV error (K=6): 0.28433
CV error (K=7): 0.27726
CV error (K=8): 0.27033
CV error (K=9): 0.26609
CV error (K=10): 0.26062
CV error (K=11): 0.25982
CV error (K=12): 0.25178 ***
CV error (K=13): 0.25850
CV error (K=14): 0.25112
CV error (K=15): 0.24938




# Colors to use for plotting:
#9F5A33 # Dark brown, secondary choice
#CCC5B5 # light grey, secondary choice
#A28B3B
#843435
#70B6E7
#797979 # grey, secondary choice
#DFBA36
#CE774A
#30A595
#6B2C4B
#D2698A
#007817
#94BE9A
#5A80A4


#####################################
#  Step 3: Plot ADMIXTURE in R
#####################################
cd /scratch/yueyu/ALL_GBS_call_on_PRA/Admixture


module load StdEnv/2020
module load r/4.2.1

R

library(tidyverse)
library(ggplot2)

fam <- read.table("GBS_all_LDpruned_dist300_numericChr.fam", header = FALSE)
sample_from_fam <- fam[,1]

samplelist <- sapply(strsplit(sample_from_fam, "_"), function(y) {
  # If the first two pieces are identical, just take the first piece
  if(length(y) >= 2 && y[1] == y[2]) {
    y[1]
  } else {
    paste(y[1:2], collapse="_")
  }
})

samplelist


all_data <- tibble(sample=character(),
                   k=numeric(),
                   Q=character(),
                   value=numeric())


for (k in 1:12){
  data <- read_delim(paste0("GBS_all_LDpruned_dist300_numericChr.",k,".Q"),
                  col_names = paste0("Q",seq(1:k)),
                  delim=" ")

  data$sample <- samplelist
  data$k <- k
  
  #This step converts from wide to long.
  data %>% gather(Q, value, -sample,-k) -> data
  all_data <- rbind(all_data,data)

  }

# -- Plot only bext K (K = 3)

# Subset best K 
all_data2 <- all_data[all_data$k == 3,]

# Define NEW color, use same in Map

cluster_colors <- c(
  "#DFBA36",
  "#6B2C4B",
  "#007817"
)

names(cluster_colors) <- paste0("Q", 1:3)
cluster_colors


# -- Reorder data for plotting (PH -> PP -> PR -> DS -> DC -> DT -> DV -> DD)

# Step 1: Create a subspecies column 
all_data2$subspecies <- gsub("[0-9_].*", "", all_data2$sample)  # subspecies number
group <- as.numeric(gsub(".*?([0-9]+)_.*", "\\1", all_data2$sample))  # e.g., 2 from DD2_1 (subspecies number)
index <- as.numeric(sub(".*_", "", all_data2$sample))                 # e.g., 1 from DD2_1 (IND number)


# Step 2:  Create a numeric score for ordering
subspecies_order <- match(all_data2$subspecies, c("PH","PP","PR","DS","DC","DT","DV","DD"))
order_score <- subspecies_order * 1e6 + group * 1e3 + index

# Step 3: Reorder sample factor levels based on the sort key
all_data2$sample <- factor(all_data2$sample, levels = unique(all_data2$sample[order(order_score)]))



# Plot
plot <- all_data2 %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack") +
  xlab("Sample") + ylab("Ancestry") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 2)) +
  scale_fill_manual(values = cluster_colors,name = "K",labels = 1:3) + 
  facet_wrap(~k,ncol=1)

# High resolution
pdf(file = "GBS_all_admixture_bestK3_with416samples_2026June05.pdf", width = 12, height = 3) 
print(plot)
dev.off()











# --------------------------------
#
#
#   Step 6: PCA
#
#
# --------------------------------


#####################################
#  Step 1: Set working directory
#####################################
# Fir
cd /home/yueyu/scratch/ALL_GBS_call_on_PRA/PCA

# Use VCF that is filtered, same as for Admixture
/scratch/yueyu/ALL_GBS_call_on_PRA/Admixture/GBS_all_LDpruned_dist300.vcf.gz



#####################################
#  Step 2: Run PCA 
#####################################

module load StdEnv/2020
module load r/4.2.1

R

# -- install packages
# install.packages("tidyverse")
# install.packages("ggplot2")
#if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#   BiocManager::install("SNPRelate")

# -- load pkgs
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(SNPRelate)
library(viridis)


#####################################
#  6.2 Load data for plotting
#####################################

setwd("/home/yueyu/scratch/ALL_GBS_call_on_PRA/PCA")

# -- VCF to GDS (IMPORTANT STEP)
#SNPRelate works with a compressed version of a genotype file called a “gds”
snpgdsVCF2GDS("/scratch/yueyu/ALL_GBS_call_on_PRA/Admixture/GBS_all_LDpruned_dist300.vcf.gz",
              "ALL_forPCA_2026May14.gds",
              method="biallelic.only")

# -- Then load GDS file
genofile <- snpgdsOpen("ALL_forPCA_2026May14.gds")


# -- Run the PCA
pca <- snpgdsPCA(genofile, eigen.cnt = 15, autosome.only = F)

    # of samples: 416
    # of SNPs: 9,594


# -- Here's the percent variance explained for each eigenvector
pc.percent <- round(pca$varprop*100,2)

  [1] 20.51  9.36  4.54  3.93  2.36  1.85  1.78  1.38  1.26  1.15  1.13  0.92
 [13]  0.81  0.77  0.72


# -- Make a dataframe of your PCA results
#Sample ID
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))


#Add population info manually 
pop_code <- substr(sample.id, 1, 3) # Extract first three characters
subpop_code <- substr(sample.id, 1, 2) # Extract first two characters
cbind(sample.id, pop_code,subpop_code)


## Make a dataframe with PC1 to PC6
tab <- data.frame(sample = pca$sample.id,
    pop = factor(pop_code)[match(pca$sample.id, sample.id)],
    subpop_code = factor(subpop_code)[match(pca$sample.id, sample.id)],
    PC1 = pca$eigenvect[,1],    # the first eigenvector
    PC2 = pca$eigenvect[,2],    # the second eigenvector
    PC3 = pca$eigenvect[,3],    
    PC4 = pca$eigenvect[,4],    
    PC5 = pca$eigenvect[,5], 
    PC6 = pca$eigenvect[,6],   
    stringsAsFactors = FALSE)

head(tab)

# -- CLEAN UP THE SAMPLE NAMES IF IT IN FORMAT OF "DC1_1_DC1_1" ETC before plotting
split_list <- strsplit(tab$sample, "_")
# Take only the first two parts for each element
tab$sample_id_clean <- sapply(split_list, function(x) paste(x[1:2], collapse = "_"))




# -------- Plot: PC 1 V.S. PC2 - Manuscript ---------
# Add: for each sample, use the highest represented Q number as color to plot in PCA

# Get highest ancestry cluster per sample
max_Q <- all_data2 %>%
  group_by(sample) %>%
  slice_max(order_by = value, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(sample, max_Q = Q)

tab2 <- tab %>%
  left_join(max_Q, by = c("sample_id_clean" = "sample"))

head(tab2)

# Use the same color palette as Admixture used above



# ======== PLOTTING FOR MANUSCTIP ===========

# ------- Plot: PC 1 V.S. PC2 - Subspecies level + NO NO NO POP NAMES ---------

pdf(file = "PC_1and2_2026June04.pdf", width = 8, height = 5)
tab2 %>%
  ggplot(.,aes(x = PC1,y = PC2, color = max_Q)) +
  geom_point(size = 6) +
  theme_classic()+
  theme(legend.text = element_text(color = "black", size = 20),
        axis.text = element_text(size=22),
        axis.title=element_text(size=28)) +
  labs(y="PC2 (9.36%)", x = "PC1 (20.51%)")+
  scale_color_manual(values = cluster_colors)  # Apply the extracted colors
dev.off()

# Completed 2026June04 -- Color same as ADMIXTURE



# ------- Plot: PC 1 V.S. PC2 - Subspecies level + NO NO NO POP NAMES + COLOR SSP ---------

subspecies_colors <- c(
  PH = "#B17CD6",
  PP = "#70B6E7",
  PR = "#5A80A4",
  DS = "#D2698A",
  DC = "#F9FE86",
  DT = "#30A595",
  DV = "#2A3F8F",
  DD = "#CE774A"
)

pdf(file = "PC_1and2_SSP_2026June04.pdf", width = 8, height = 5)
tab2 %>%
  ggplot(aes(x = PC1, y = PC2, color = subpop_code)) +
  geom_point(size = 6, alpha = 0.8) +
  scale_color_manual(values = subspecies_colors) +
  theme_classic() +
  theme(
    legend.text = element_text(size = 20),
    legend.title = element_blank(),
    axis.text = element_text(size = 22),
    axis.title = element_text(size = 28)
  ) +
  labs(y = "PC2 (9.36%)", x = "PC1 (20.51%)")
dev.off()

# Completed 2026June04 -- Color same as SSP

