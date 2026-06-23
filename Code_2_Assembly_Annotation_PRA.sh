#####################################
# Project: debilis_praecox_project
#
# Code 2: Genome Assembly & Annotation for PRA
#
# Date: 2024 Oct
#
# by: Yue Yu
#
#####################################


# Including all the following part:
# -- Part 1: Juicer
# -- Part 2: Run YAHS on Hap1 and Hap2 seperately
# -- Part 3: Minimap
# -- Part 4: Merged 1+2 -> Juicer again
# -- Part 5: Reorder and rename CHR based on HA412
# -- Part 6: Check assembly readiness: BUSCO
# -- Part 7: Final Minimap
# -- Part 8: Final Minimap between DEB hap2 - PRA hap2







#####################################
#
# Part 1: Juicer
# 
#####################################

#-----------
#   Hap 1
#-----------
cd Ready_Hap1
ln -s /home/yueyu/scratch/GA/juicer/CPU/ scripts
ln -s /home/yueyu/scratch/PRA_GA/fastq fastq

nano Pra_Hap01_juicer.sh

#----------------- Pra_Hap01_juicer.sh --------------

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=def-rieseber
#SBATCH --time=3-0
#SBATCH --ntasks-per-node=32
#SBATCH --mem=149G

cd /home/yueyu/scratch/PRA_GA/Ready_Hap1

module load StdEnv/2020 bwa/0.7.17 java/17.0.2 samtools/1.15.1

bash scripts/juicer.sh -D $PWD \
					   -g contigassembly \
					   -s DpnII \
					   -p restriction_sites/pra_hap1_DpnII.chrom.sizes \
					   -y restriction_sites/pra_hap1_DpnII.txt \
					   -z reference/Hpra677_HIC_oldhifiasm_S4.hic.hap1.fasta \
					   -t 32 \
					   -S early

#----------------- Pra_Hap01_juicer.sh (END) --------------


#Format of the directory structure

total 88K
drwxr-x---. 2 yueyu 25K Oct 24 15:33 reference
drwxr-x---. 2 yueyu 25K Oct 24 16:42 restriction_sites
lrwxrwxrwx. 1 yueyu  34 Oct 24 16:48 scripts -> /home/yueyu/scratch/GA/juicer/CPU/
lrwxrwxrwx. 1 yueyu  32 Oct 24 16:48 fastq -> /home/yueyu/scratch/PRA_GA/fastq
-rw-r-----. 1 yueyu 520 Oct 24 16:57 Pra_Hap01_juicer.sh


# Run Juicer!!!
chmod +x Pra_Hap01_juicer.sh
sbatch Pra_Hap01_juicer.sh


#-----------
#   Hap 2
#-----------
cd Ready_Hap2
ln -s /home/yueyu/scratch/GA/juicer/CPU/ scripts
ln -s /home/yueyu/scratch/PRA_GA/fastq fastq

nano Pra_Hap02_juicer.sh

#----------------- Pra_Hap02_juicer.sh --------------

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=def-rieseber
#SBATCH --time=3-0
#SBATCH --ntasks-per-node=32
#SBATCH --mem=149G

cd /home/yueyu/scratch/PRA_GA/Ready_Hap2

module load StdEnv/2020 bwa/0.7.17 java/17.0.2 samtools/1.15.1

bash scripts/juicer.sh -D $PWD \
					   -g contigassembly \
					   -s DpnII \
					   -p restriction_sites/pra_hap2_DpnII.chrom.sizes \
					   -y restriction_sites/pra_hap2_DpnII.txt \
					   -z reference/Hpra677_HIC_oldhifiasm_S4.hic.hap2.fasta \
					   -t 32 \
					   -S early

#----------------- Pra_Hap02_juicer.sh (END) --------------


#Format of the directory structure
ls -thor
total 88K
drwxr-x---. 2 yueyu 25K Oct 24 15:36 reference
drwxr-x---. 2 yueyu 25K Oct 24 16:42 restriction_sites
lrwxrwxrwx. 1 yueyu  34 Oct 24 16:48 scripts -> /home/yueyu/scratch/GA/juicer/CPU/
lrwxrwxrwx. 1 yueyu  32 Oct 24 16:48 fastq -> /home/yueyu/scratch/PRA_GA/fastq
-rwxr-x---. 1 yueyu 520 Oct 24 17:03 Pra_Hap02_juicer.sh


# Run Juicer!!!
chmod +x Pra_Hap02_juicer.sh
sbatch Pra_Hap02_juicer.sh



#####################################
#
# Part 2: Run YAHS on Hap1 and Hap2 seperately
# 
#####################################


# ---------------
#
# Praecox Hap 1
#
# ---------------

cd /home/yueyu/scratch/PRA_GA/Ready_Hap1/aligned

nano prep_for_Yahs.sh
# --------- prep_for_Yahs.sh --------
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=0-3
#SBATCH --mem=20G

cd /home/yueyu/scratch/PRA_GA/Ready_Hap1/aligned

awk '
{
     if ($9 > 0 && $12 > 0) {
    # Calculate end position for read 1 using pos1 and cigar1
    cigar = $10
    pos_start1 = $3
    sum1 = 0
    while (match(cigar, /([0-9]+)([MDNX=])/, arr)) {
        sum1 += arr[1]
        cigar = substr(cigar, RSTART + RLENGTH)
    }
    pos_end1 = pos_start1 + sum1  # End position for read 1

    # Calculate end position for read 2 using pos2 and cigar2
    cigar = $13
    pos_start2 = $7
    sum2 = 0
    while (match(cigar, /([0-9]+)([MDNX=])/, arr)) {
        sum2 += arr[1]
        cigar = substr(cigar, RSTART + RLENGTH)
    }
    pos_end2 = pos_start2 + sum2  # End position for read 2

    # Print interleaved output for read 1 and read 2
    print $2, pos_start1, pos_end1, $15"/1", $9
    print $6, pos_start2, pos_end2, $16"/2", $12
}
}
' merged_nodups.txt > yahs_result/merged_nodups_for_yahs_Hap1.bed

# --------- prep_for_Yahs.sh (END) --------

# Submitted batch job 39379431 (Hap1)
# Done 36 min






# ======= Step 2: download and install yahs (DONE)
# to run it
/home/yueyu/scratch/software/yahs/yahs


# ======= Step 3: run yahs

nano run_Yahs.sh
# --------- run_Yahs.sh --------
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=def-rieseber
#SBATCH --time=0-1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=10G

cd /home/yueyu/scratch/PRA_GA/Ready_Hap1/aligned/yahs_result

CONTIG="/home/yueyu/scratch/PRA_GA/Ready_Hap1/reference/Hpra677_HIC_oldhifiasm_S4.hic.hap1.fasta"
# This location need to have fasta.fai file to accompany the fasta file
HIC_TO_CONTIG="merged_nodups_for_yahs_Hap1.bed" # same

/home/yueyu/scratch/software/yahs/yahs $CONTIG $HIC_TO_CONTIG
# --------- run_Yahs.sh (END) --------



# ======= Step 4:  generate a HiC contact map file > for JuiceBox

# ------- Step 4.1: generate .assembly file

# The following can dirctly run in command line without any resource request
CONTIG_FAI="/home/yueyu/scratch/PRA_GA/Ready_Hap1/reference/Hpra677_HIC_oldhifiasm_S4.hic.hap1.fasta.fai"
/home/yueyu/scratch/software/yahs/juicer pre -a -o out_JBAT yahs.out.bin yahs.out_scaffolds_final.agp $CONTIG_FAI >out_JBAT.log 2>&1

#2>&1: means == 2 (stderr) is being redirected to 1 (stdout), effectively merging both output streams

# ------- Step 4.2: generate .hic file

nano run_generateHic.sh
# -------- run_generateHic.sh --------
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=def-rieseber
#SBATCH --time=0-1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=64G

cd /home/yueyu/scratch/PRA_GA/Ready_Hap1/aligned/yahs_result

module load StdEnv/2020 python/3.11.2 java/17.0.2 lastz/1.04.03

(java -jar -Xmx60G /home/yueyu/scratch/software/juicer_jar/juicer_tools.1.9.9_jcuda.0.8.jar pre out_JBAT.txt out_JBAT.hic.part <(cat out_JBAT.log  | grep PRE_C_SIZE | awk '{print $2" "$3}')) && (mv out_JBAT.hic.part out_JBAT.hic)

# --------- run_generateHic.sh (END) --------



# ------- Step 5: Move to Juicebox to check 
cd ~/Desktop/test_3d_DNA/yahs_pra/yahs_pra_hap1
scp yueyu@narval.computecanada.ca:"/home/yueyu/scratch/PRA_GA/Ready_Hap1/aligned/yahs_result/out_JBAT.hic" .
scp yueyu@narval.computecanada.ca:"/home/yueyu/scratch/PRA_GA/Ready_Hap1/aligned/yahs_result/out_JBAT.assembly" .


# ------- Step 6: Extract first 17 CHR
module load seqkit
seqkit head -n 17 yahs.out_scaffolds_final.fa > Pra_HAP_1_first_17_super_scaffolds.fasta
grep -c "^>" Pra_HAP_1_first_17_super_scaffolds.fasta

# ---------------
#
# Praecox Hap 2 (same as above)
#
# ---------------



#####################################
#
# Part 3: Minimap
# 
#####################################

# Steps to take
# -- Step 1: Run minimap to PAF format
# -- Step 2: Plot PAF result -> to guide futrue steps
# -- Step 3: use PAF result -> Rename CHR to match in H1 and H2
# -- Step 4: use PAF result -> Reverse compliment
# -- Step 5: Run minimap to BAM format (corrected CHR name and RC fasta)
# -- Step 6: Run Syri with bam file
# -- Step 7: Run PlotSR with syri output




#---------------------------------
# Step 1: Run minimap to PAF format
#---------------------------------

cd /home/yueyu/scratch/PRA_GA/YAHS_MINIMAP_PRA/PAF

cat run_pra_minimap_yahsHap1_yahsHap2_20250115_topaf.sh
# ------------------------- (START) ---------------
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=3-0
#SBATCH --cpus-per-task=1
#SBATCH --mem=300G

module load StdEnv/2023
module load minimap2/2.28

pra_hap1="/home/yueyu/scratch/PRA_GA/Ready_Hap1/aligned/yahs_result/Pra_HAP_1_first_17_super_scaffolds.fasta"
pra_hap2="/home/yueyu/scratch/PRA_GA/Ready_Hap2/aligned/yahs_result/Pra_HAP_2_first_17_super_scaffolds.fasta"


cd /home/yueyu/scratch/PRA_GA/YAHS_MINIMAP_PRA
minimap2 -cx asm5 $pra_hap1 $pra_hap2 > pra_minimap_yahsHap1_yahsHap2_20250115.paf


for i in $(ls *paf)
do

name=$(echo $i | cut -d "." -f 1)

awk '$12 >= 30 {print $8, $9, $3, $4, "30", $1, $6, $5 }' $i > "$name""MQ30_tab.txt"

done


# ---------------------------------------------------------------------
echo "Job finished with exit code $? at: `date`"
# ---------------------------------------------------------------------
# -------------------------(END) ---------------

Submitted batch job 39496519
# DONE





#---------------------------------
# -- Step 2: Plot PAF result 
#---------------------------------
module load StdEnv/2023
module load r/4.4.0
cd /home/yueyu/scratch/PRA_GA/YAHS_MINIMAP_PRA

R


nano minimap_plot.R
#======= minimap_plot.R =============

library(dplyr)
library(magrittr)
library(knitr)
library(ggplot2)
library(tidyr)
library(tidyverse)

# Functions

filterMum <- function(df, minl=1000, flanks=1e4){
    coord = df %>% filter(abs(re-rs)>minl) %>% group_by(qid, rid) %>%
        summarize(qsL=min(qs)-flanks, qeL=max(qe)+flanks, rs=median(rs)) %>%
        ungroup %>% arrange(desc(rs)) %>%
        mutate(qid=factor(qid, levels=unique(qid))) %>% select(-rs)
    merge(df, coord) %>% filter(qs>qsL, qe<qeL) %>%
        mutate(qid=factor(qid, levels=levels(coord$qid))) %>% select(-qsL, -qeL)
}



diagMum <- function(df){
    ## Find best qid order
    rid.o = df %>% group_by(qid, rid) %>% summarize(base=sum(abs(qe-qs)),
                                                    rs=weighted.mean(rs, abs(qe-qs))) %>%
        ungroup %>% arrange(desc(base)) %>% group_by(qid) %>% do(head(., 1)) %>%
        ungroup %>% arrange(desc(rid), desc(rs)) %>%
        mutate(qid=factor(qid, levels=unique(qid)))
    ## Find best qid strand
    major.strand = df %>% group_by(qid) %>%
        summarize(major.strand=ifelse(sum(sign(qe-qs)*abs(qe-qs))>0, '+', '-'),
                  maxQ=max(c(qe, qs)))
    merge(df, major.strand) %>% mutate(qs=ifelse(major.strand=='-', maxQ-qs, qs),
                                       qe=ifelse(major.strand=='-', maxQ-qe, qe),
                                       qid=factor(qid, levels=levels(rid.o$qid)))
}



#Load data
f <- "pra_minimap_yahsHap1_yahsHap2_20250115MQ30_tab.txt"
name <- gsub(".txt", "", f)
mumgp <- read.table(f, header=F)
colnames(mumgp) <- c("rs", "re", "qs", "qe", "error", "qid", "rid", "strand")
dim(mumgp)
head(mumgp)

mumgp[c("qs", "qe")] <- t(mapply(\(a, b, c, d, e, f, g, h){
                 if(h == "-") c(d,c) else c(c,d) }, mumgp$rs, mumgp$re, mumgp$qs, mumgp$qe, mumgp$error, mumgp$qid, mumgp$rid, mumgp$strand))

mumgp.filt = filterMum(mumgp, minl=1e4)

mumgp.filt.diag = diagMum(mumgp.filt)

head(mumgp.filt.diag)
dim(mumgp.filt.diag)

#Plotting
P1 <- ggplot(mumgp.filt.diag, aes(x=rs, xend=re, y=qs, yend=qe, colour=strand)) +
  geom_segment(show.legend=FALSE, size=3) + geom_point(alpha=0.09) + theme_bw() +
  facet_grid(qid~rid, scales='free', space='free', switch='both') +
  guides(colour=guide_legend(override.aes=list(alpha=1))) +
  theme(strip.text.y=element_text(angle=120, size=10),
        strip.text.x=element_text(angle=60, size=10),
        strip.background=element_blank(),
        legend.position=c(1,-.03), legend.justification=c(1,1),
        legend.direction='horizontal',
        axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        panel.spacing=unit(0, 'cm')) +
  xlab('reference sequence') +
  ylab('assembly') +
  scale_colour_brewer(palette='Set1')


png(paste0(name, ".png"), width=2000, height= 2000)
print(P1)
dev.off()

#======= minimap_plot.R (END) =============

cd /home/yueyu/scratch/PRA_GA/YAHS_MINIMAP_PRA
scp yueyu@narval.computecanada.ca:"/home/yueyu/scratch/PRA_GA/YAHS_MINIMAP_PRA/pra_minimap_yahsHap1_yahsHap2_20250115MQ30_tab.png" .

# DONE - plotting - saved at local laptop too






#---------------------------------
# -- Step 3: Rename CHR to match in H1 and H2
#---------------------------------

# Make sure chromosome number is the same in both files
cd /home/yueyu/scratch/PRA_GA

#-------- Change CHR names (HAP 1)

nano Change_Hap1_toCHR_names_yahs_PRA.txt

scaffold_1   CHR01
scaffold_10   CHR02
scaffold_11   CHR03
scaffold_12   CHR04
scaffold_13   CHR05
scaffold_14   CHR06
scaffold_15   CHR07
scaffold_16   CHR08
scaffold_17  CHR09
scaffold_2   CHR10
scaffold_3   CHR11
scaffold_4   CHR12
scaffold_5   CHR13
scaffold_6   CHR14
scaffold_7   CHR15
scaffold_8   CHR16
scaffold_9   CHR17

cd /home/yueyu/scratch/PRA_GA/YAHS_MINIMAP_PRA
cp /home/yueyu/scratch/PRA_GA/Ready_Hap1/aligned/yahs_result/Pra_HAP_1_first_17_super_scaffolds.fasta Pra_HAP_1_first_17_super_scaffolds_CHR_RENAMED.fasta

Hap1fasta="Pra_HAP_1_first_17_super_scaffolds_CHR_RENAMED.fasta"

while read -r old_name new_name; do
    sed -i "s/\b$old_name\b/$new_name/g" $Hap1fasta
done < Change_Hap1_toCHR_names_yahs_PRA.txt

# DONE



#-------- Change CHR names (HAP 2)

nano Change_Hap2_toCHR_names_yahs_PRA.txt

scaffold_1   CHR01
scaffold_3   CHR02
scaffold_10   CHR03
scaffold_14   CHR04
scaffold_11   CHR05
scaffold_16   CHR06
scaffold_15   CHR07
scaffold_17   CHR08
scaffold_13  CHR09
scaffold_5   CHR10
scaffold_2   CHR11
scaffold_4   CHR12
scaffold_8   CHR13
scaffold_6   CHR14
scaffold_9   CHR15
scaffold_12   CHR16
scaffold_7   CHR17

cd /home/yueyu/scratch/PRA_GA/YAHS_MINIMAP_PRA
cp /home/yueyu/scratch/PRA_GA/Ready_Hap2/aligned/yahs_result/Pra_HAP_2_first_17_super_scaffolds.fasta Pra_HAP_2_first_17_super_scaffolds_CHR_RENAMED.fasta

Hap2fasta="Pra_HAP_2_first_17_super_scaffolds_CHR_RENAMED.fasta"

while read -r old_name new_name; do
    sed -i "s/\b$old_name\b/$new_name/g" $Hap2fasta
done < Change_Hap2_toCHR_names_yahs_PRA.txt




#---------------------------------
# -- Step 4: Reverse compliment
#---------------------------------


module load StdEnv/2023
module load r/4.4.0
cd /home/yueyu/scratch/PRA_GA/YAHS_MINIMAP_PRA

R

#------ r (start) -------- 
library(Biostrings)

#------ Hap 1 ------
# Define file paths
input_fasta <- "Pra_HAP_1_first_17_super_scaffolds_CHR_RENAMED.fasta"   # Input FASTA file
output_fasta <- "Pra_HAP_1_first_17_super_scaffolds_CHR_RENAMED_RC.fasta" # Output FASTA file

# Read the sequences from the FASTA file
sequences <- readDNAStringSet(input_fasta)

sequences

# Define the headers of sequences you want to reverse complement 
# must be the ones in default RED in minimap PAF plot!!
# RED means need to be Reversed

headers_to_reverse_complement <- c("CHR02","CHR04","CHR05","CHR07","CHR09","CHR12","CHR13","CHR14","CHR16","CHR17")


# Reverse complement the sequences with matching headers
for (header in headers_to_reverse_complement) {
    if (header %in% names(sequences)) {
        sequences[[header]] <- reverseComplement(sequences[[header]])
    }
}

sequences

# Write the modified sequences to a new FASTA file
writeXStringSet(sequences, output_fasta)


#------ Hap 2 (Do not need to adjust) ------
# because hap1 is already reverse cimplimented, Hap2 do not need to change

#------ r (end) --------




#---------------------------------
# -- Step 5: Run minimap to BAM format (corrected CHR name and RC fasta)
#---------------------------------


cd /home/yueyu/scratch/PRA_GA/YAHS_MINIMAP_PRA/BAM_2025feb23

nano run_pra_minimap_yahsHap1_yahsHap2_20250222_tobam.sh

# -------------------
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=1-08:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=300G


module load StdEnv/2023
module load minimap2/2.28

module load gcc/12.3
module load samtools/1.20


cd /home/yueyu/scratch/PRA_GA/YAHS_MINIMAP_PRA

pra_hap1="Pra_HAP_1_first_17_super_scaffolds_CHR_RENAMED_RC.fasta"
pra_hap2="Pra_HAP_2_first_17_super_scaffolds_CHR_RENAMED.fasta"


minimap2 -ax asm5 --eqx $pra_hap1 $pra_hap2 > BAM_2025feb23/pra_minimap_yahs_20250222.sam
# asm5: average divergence is << 5%.
# eqx: Output =/X CIGAR operators for sequence match/mismatch.

cd BAM_2025feb23
samtools view -b pra_minimap_yahs_20250222.sam > pra_minimap_yahs_20250222.bam

# ---------------------------------------------------------------------
echo "Job finished with exit code $? at: `date`"
# ---------------------------------------------------------------------
# -------------------

# Submitted batch job 40678771 - Tue 25 Feb 2025 01:51:54 PM EST
# DONE in 1 day 06 hours




#---------------------------------
# -- Step 6: Run Syri with bam file
#---------------------------------

#narval3
tmux new-session -s syri_narval3
tmux attach-session -t syri_narval3


module load StdEnv/2020 gcc python/3.9 igraph
source ~/ENV/bin/activate

#=======
#  Syri 
#=======

cd /home/yueyu/scratch/PRA_GA/YAHS_MINIMAP_PRA/BAM_2025feb23

pra_hap1="/home/yueyu/scratch/PRA_GA/YAHS_MINIMAP_PRA/Pra_HAP_1_first_17_super_scaffolds_CHR_RENAMED_RC.fasta"
pra_hap2="/home/yueyu/scratch/PRA_GA/YAHS_MINIMAP_PRA/Pra_HAP_2_first_17_super_scaffolds_CHR_RENAMED.fasta"

syri -c pra_minimap_yahs_20250222.bam -r $pra_hap1 -q $pra_hap2 -F B --prefix syri_20250227

# DONE - 14 hours




#---------------------------------
# -- Step 7: Run PlotSR with syri output
#---------------------------------

#=========
#  plotSR
#=========

nano genomes_pra_20250227.txt
# ------- genomes_pra.txt ------
#file   name    tags
/home/yueyu/scratch/PRA_GA/YAHS_MINIMAP_PRA/Pra_HAP_1_first_17_super_scaffolds_CHR_RENAMED_RC.fasta   PRA_Hap_1 lw:1.5
/home/yueyu/scratch/PRA_GA/YAHS_MINIMAP_PRA/Pra_HAP_2_first_17_super_scaffolds_CHR_RENAMED.fasta   PRA_Hap_2 lw:1.5
# ------- genomes_pra.txt ------


plotsr --sr syri_20250227syri.out --genomes genomes_pra_20250227.txt -H 8 -W 5

# genomes_pra.txt file should eb seperated by TAB not SPACE
# -H and -W stand for height and width of output graph.pdf/png/other
# Run time : 5min (for running the entire genome together)


cd ~/Desktop/REF_DEB_PRA
scp yueyu@narval.computecanada.ca:"/home/yueyu/scratch/GA/PAF/deb_minimap_20240920MQ30_tab.png" .
# Done 2025 Feb 28th




#---------------------------------
# -- Step 8: extract INV > 500Kbp 
#---------------------------------

# -- Extract all inversions
awk -F'\t' '$11 == "INV"' syri_20250227syri.out > syri_INV_all.out
wc -l syri_INV_all
# 317 inversions detected (matches result in the summary file)


# -- Extract inversions > 500Kbp
awk -F'\t' '($3 - $2) >= 500000' syri_INV_all.out > syri_INV_over500Kbp.out
wc -l syri_INV_over500Kbp.out
# 55 inversions > 500Kbp

# -- Copy paste 55 rows into EXCEL for manual curation (2025 Feb 27th)



#---------------------------------
# -- Step 9: SYRI summary
#---------------------------------
# Syri output
cat syri_20250227syri.summary








#####################################
#
# Part 4: Merged 1+2 -> Juicer again
# 
#####################################

#=====================
#
#  Prepare H1 + H2 
#
#=====================

# Merge and run Juicer again

#---- Make file directory to run Juicer
cd /home/yueyu/scratch/PRA_GA/YAHS_Merged_Hap1_and_Hap2_PRA


#---- FILE PREP (1/4): soft link to fastq
ln -s /home/yueyu/scratch/PRA_GA/fastq 


#---- FILE PREP (2/4): Scripts
ln -s /home/yueyu/scratch/GA/juicer/CPU/ scripts

scripts/juicer.sh 


#---- FILE PREP (3/4): Merge FASTA files
mkdir reference
cd reference

HAP1="/home/yueyu/scratch/PRA_GA/Ready_Hap1/aligned/yahs_result/Pra_HAP_1_first_17_super_scaffolds.fasta"
HAP2="/home/yueyu/scratch/PRA_GA/Ready_Hap2/aligned/yahs_result/Pra_HAP_2_first_17_super_scaffolds.fasta"

# 加这一步是因为之前yahs组装出来的初步fasta里面的scaffold命名为（scaffold_1）不区分是在hap1或hap2
# 加上这一步后，在merge fasta的时候可以在后续的步骤中看到是哪一个hap中的super scaffold
sed '/^>/s/>/>h1_HiC_/' $HAP1 > Pra_HAP_1_first_17_super_scaffolds_names_added_hap1.fasta
sed '/^>/s/>/>h2_HiC_/' $HAP2 > Pra_HAP_2_first_17_super_scaffolds_names_added_hap2.fasta

# Merge FASTA
cat Pra_HAP_1_first_17_super_scaffolds_names_added_hap1.fasta Pra_HAP_2_first_17_super_scaffolds_names_added_hap2.fasta > PRA_mergedHap1and2_20250227.fasta

grep ">" PRA_mergedHap1and2_20250227.fasta



#---- Index FASTA files

nano index.sh
#----------------- index.sh --------------
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=def-rieseber
#SBATCH --time=05:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=30G

cd /home/yueyu/scratch/PRA_GA/YAHS_Merged_Hap1_and_Hap2_PRA/reference

module load StdEnv/2023
module load bwa/0.7.18

bwa index PRA_mergedHap1and2_20250227.fasta
#Took ~2 hours to do this step
#----------------- index.sh (END) --------------
# Submitted batch job 40792565  - DONE




#---- FILE PREP (4/4): Resitriction sites

mkdir restriction_sites
cd restriction_sites

# generate_site_positions.py (Download)
wget -O generate_site_positions.py https://raw.githubusercontent.com/aidenlab/juicer/main/misc/generate_site_positions.py


nano generate_site_positions.py

#Add to this section the path to your fasta (Modify)
filenames = {
    'hg19': '/seq/reference/Homo_sapiens_assembly19.fasta',
    'PraHap1and2': '../reference/PRA_mergedHap1and2_20250227.fasta', #here I am adding the path to my fasta
    'mm9' : '/seq/reference/Mus_musculus_assembly9.fasta',
    'mm10': '/seq/reference/Mus_musculus_assembly10.fasta',
    'hg18': '/seq/reference/Homo_sapiens_assembly18.fasta',
  }
  
mv generate_site_positions.py generate_site_positions_PRA_Hap1and2.py

#run generate_site_positions.py

nano run_python.sh
#----------------- run_python.sh --------------

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=def-rieseber
#SBATCH --time=3:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=5G

cd /home/yueyu/scratch/PRA_GA/YAHS_Merged_Hap1_and_Hap2_PRA/restriction_sites

module load StdEnv/2023 
module load python/3.12.4

python generate_site_positions_PRA_Hap1and2.py DpnII PraHap1and2

#----------------- run_python.sh (END) --------------
# ~41 min to run
# Submitted batch job 40792606 - DONE



cd restriction_sites
nano Chromosome_sizes.sh

#------- Chromosome_sizes.sh
for i in $(ls *DpnII.txt)
do

name=$(echo $i | cut -d "." -f 1 )

awk 'BEGIN{OFS="\t"}{print $1, $NF}'  $i > "$name"".chrom.sizes"

done

# {print $1, $NF}: For each line in the file $i, this prints the first field ($1) and the last field ($NF).
#------- Chromosome_sizes.sh (END)

bash Chromosome_sizes.sh
wc -l PraHap1and2_DpnII.chrom.sizes
#34 chromosomes: adding up both 17 chomosomes
# DONE



#=====================
#
#  Run Juicer Hap 1 + Hap2
#
#=====================

cd /home/yueyu/scratch/PRA_GA/YAHS_Merged_Hap1_and_Hap2_PRA


# Data Structure
lrwxrwxrwx. 1 yueyu  32 Feb 28 17:33 fastq -> /home/yueyu/scratch/PRA_GA/fastq
drwxr-x---. 2 yueyu 25K Feb 28 19:10 restriction_sites
drwxr-x---. 2 yueyu 25K Feb 28 20:11 reference
lrwxrwxrwx. 1 yueyu  34 Feb 28 21:11 scripts -> /home/yueyu/scratch/GA/juicer/CPU/


nano PRA_BOTH_HAP_juicer_2025Feb.sh
#----------------- PRA_BOTH_HAP_juicer_2025Feb.sh --------------
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=def-rieseber
#SBATCH --time=4-0
#SBATCH --ntasks-per-node=32
#SBATCH --mem=200G

cd /home/yueyu/scratch/PRA_GA/YAHS_Merged_Hap1_and_Hap2_PRA

module load StdEnv/2020 bwa/0.7.17 java/17.0.2 samtools/1.15.1

bash scripts/juicer.sh -D $PWD \
                       -g contigassembly \
                       -s DpnII \
                       -p restriction_sites/PraHap1and2_DpnII.chrom.sizes \
                       -y restriction_sites/PraHap1and2_DpnII.txt \
                       -z reference/PRA_mergedHap1and2_20250227.fasta \
                       -t 32 \
                       -S early

#Check the name: reference or referenceS

#----------------- PRA_BOTH_HAP_juicer_2025Feb.sh (END) --------------

chmod +x PRA_BOTH_HAP_juicer_2025Feb.sh
sbatch PRA_BOTH_HAP_juicer_2025Feb.sh

Submitted batch job 40798186  - DONE 2025 FEB 28TH

seff 40798186
Job ID: 40798186
Cluster: narval
User/Group: yueyu/yueyu
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 32
CPU Utilized: 11-00:00:57
CPU Efficiency: 41.83% of 26-07:07:44 core-walltime
Job Wall-clock time: 19:43:22
Memory Utilized: 185.73 GB
Memory Efficiency: 92.86% of 200.00 GB


# ---- Output in 

cd /home/yueyu/scratch/PRA_GA/YAHS_Merged_Hap1_and_Hap2_PRA/aligned
le merged_nodups.txt




#=====================
#  Run YAHS Hap 1 + Hap2 （ do not do this for hap 1+2)
#  Note: !! DO NOT USE!! ran it on 2025 March 3rd and result in weird assembly as it corrects the assembly again
#=====================


#=====================
#
#  Run 3D-DNA Hap 1 + Hap2
#
#=====================

cd /home/yueyu/scratch/PRA_GA/YAHS_Merged_Hap1_and_Hap2_PRA/aligned
mkdir 3D_DNA_run_H1_and2

nano Run_Hap_BOTH_with_3D_DNA.sh

#----------------- Run_Hap_BOTH_with_3D_DNA.sh --------------

#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=7-0
#SBATCH --cpus-per-task=15
#SBATCH --mem=100G

cd /home/yueyu/scratch/PRA_GA/YAHS_Merged_Hap1_and_Hap2_PRA

module load StdEnv/2020 python/3.10.2 java/17.0.2 lastz/1.04.03

export PATH="/Hap_BOTH_3D_DNA/3d-dna:$PATH"

source /home/yueyu/scratch/GA/3d-dna/3ddna/bin/activate

/home/yueyu/scratch/GA/3d-dna/run-asm-pipeline.sh -r 0 reference/PRA_mergedHap1and2_20250227.fasta aligned/merged_nodups.txt

deactivate
#----------------- Run_Hap_BOTH_with_3D_DNA.sh (END) --------------


#=====================
#
#  Move to Juicebox to check 
#
#=====================

cd /Users/yueyu/Desktop/test_3d_DNA/yahs_pra/CORRECT_3D_DNA_pra_hap1and2_2025March
scp yueyu@narval.computecanada.ca:"/home/yueyu/scratch/PRA_GA/YAHS_Merged_Hap1_and_Hap2_PRA/PRA_mergedHap1and2_20250227.0.hic" .
scp yueyu@narval.computecanada.ca:"/home/yueyu/scratch/PRA_GA/YAHS_Merged_Hap1_and_Hap2_PRA/PRA_mergedHap1and2_20250227.0.assembly" .


#=====================
#
# Manual Curation Hap 1+ Hap2
# Check all INV (>500Kbp)
#
#   Done 2024 March 10th
#
#=====================

# Checked location in
PRA_Syri_2025March10.xlsx

# Checked screenshot in
Ref_PRA_2024Dec.pptx

# Reviewed assembly in
cd /Users/yueyu/Desktop/test_3d_DNA/yahs_pra/CORRECT_3D_DNA_pra_hap1and2_2025March
less -SN PRA_FINAL_2025March10_mergedHap1and2_20250227.0.review.assembly

# Trasport to Narval from Local laptop
scp PRA_FINAL_2025March10_mergedHap1and2_20250227.0.review.assembly yueyu@narval.computecanada.ca:"/home/yueyu/scratch/PRA_GA/YAHS_Merged_Hap1_and_Hap2_PRA/reviewed_assembly/"




#=====================
#
#  review.assembly  --> FASTA (with Moj/Eric R.script)
#
#=====================

cd /home/yueyu/scratch/PRA_GA/YAHS_Merged_Hap1_and_Hap2_PRA/reviewed_assembly

# Final assembly 
#1. MIX_ASSEM <- args[1] #Reviewed assembly file
#2. PREFIX <- args[2]    #a random name to add to output
#3. MIX_FASTA <- args[3] #FASTA used to run Juicer (in the reference file)
#4. SAVE_DIR <- args[4]

#----- Hap 1+2 
MIX_ASSEM="PRA_FINAL_2025March10_mergedHap1and2_20250227.0.review.assembly"
PREFIX="praecox_hap1_and_hap2"
MIX_FASTA="/home/yueyu/scratch/PRA_GA/YAHS_Merged_Hap1_and_Hap2_PRA/reference/PRA_mergedHap1and2_20250227.fasta"

module load StdEnv/2023
module load r/4.4.0
cd /home/yueyu/scratch/PRA_GA/YAHS_Merged_Hap1_and_Hap2_PRA/reviewed_assembly

#Run this in a new screen 
ASM_TO_FASTA_Yue_for_Hap1_and_Hap2_forPRACOX_2025Mar10.R


# ------ Two important outputs from R to fasta script
cd /home/yueyu/scratch/PRA_GA/YAHS_Merged_Hap1_and_Hap2_PRA/reviewed_assembly

# RAW: Hap 1
Praecox_20250310_Yue_hap1_save.fasta
# 18 chromosomes

# RAW: Hap 2
Praecox_20250310_Yue_hap2_save.fasta
# 20 chromosomes


# ------- Extract first 17 CHR
cd /home/yueyu/scratch/PRA_GA/YAHS_Merged_Hap1_and_Hap2_PRA/reviewed_assembly

module load seqkit/2.5.1

seqkit grep -f Hap1_17scaffolds_to_extract.txt Praecox_20250310_Yue_hap1_save.fasta -o Praecox_20250310_Yue_hap1_save_17CHR.fasta
seqkit grep -f Hap2_17scaffolds_to_extract.txt Praecox_20250310_Yue_hap2_save.fasta -o Praecox_20250310_Yue_hap2_save_17CHR.fasta





#####################################
#
# Part 5: Reorder and rename CHR based on HA412
# 
#####################################

#-------------------------------------------------
#
# (Hap 1) Step 1: Run minimap to PAF format (PRA H1 - HA412)
#
#-------------------------------------------------

cd /home/yueyu/scratch/PRA_FASTA_almostdone

nano run_pra_minimap_finalHap1_HA412_20250312_topaf.sh
# ------------------------- (START) ---------------
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=3-0
#SBATCH --cpus-per-task=1
#SBATCH --mem=300G

module load StdEnv/2023
module load minimap2/2.28

HA412="/home/yueyu/projects/def-rieseber/yueyu/HA412_REF/Ha412HOv2.0-20181130.fasta"
pra_hap1="/home/yueyu/scratch/PRA_GA/YAHS_Merged_Hap1_and_Hap2_PRA/reviewed_assembly/Praecox_20250310_Yue_hap1_save_17CHR.fasta"

cd /home/yueyu/scratch/PRA_FASTA_almostdone
minimap2 -cx asm5 $HA412 $pra_hap1 > pra_minimap_Hap1_HA412_20250312.paf


for i in $(ls *paf)
do

name=$(echo $i | cut -d "." -f 1)

awk '$12 >= 30 {print $8, $9, $3, $4, "30", $1, $6, $5 }' $i > "$name""MQ30_tab.txt"

done


# ---------------------------------------------------------------------
echo "Job finished with exit code $? at: `date`"
# ---------------------------------------------------------------------
# -------------------------(END) ---------------

Submitted batch job 41195089
# DONE - 2024 March 13th




#-------------------------------------------------
# (Hap 2) Step 1: Run minimap to PAF format (PRA H2 - HA412)
#-------------------------------------------------

cd /home/yueyu/scratch/PRA_FASTA_almostdone

nano run_pra_minimap_finalHap2_HA412_20250312_topaf.sh
# ------------------------- (START) ---------------
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=3-0
#SBATCH --cpus-per-task=1
#SBATCH --mem=300G

module load StdEnv/2023
module load minimap2/2.28

HA412="/home/yueyu/projects/def-rieseber/yueyu/HA412_REF/Ha412HOv2.0-20181130.fasta"
pra_hap2="/home/yueyu/scratch/PRA_GA/YAHS_Merged_Hap1_and_Hap2_PRA/reviewed_assembly/Praecox_20250310_Yue_hap2_save_17CHR.fasta"

cd /home/yueyu/scratch/PRA_FASTA_almostdone
minimap2 -cx asm5 $HA412 $pra_hap2 > pra_minimap_Hap2_HA412_20250312.paf


for i in $(ls *paf)
do

name=$(echo $i | cut -d "." -f 1)

awk '$12 >= 30 {print $8, $9, $3, $4, "30", $1, $6, $5 }' $i > "$name""MQ30_tab.txt"

done


# ---------------------------------------------------------------------
echo "Job finished with exit code $? at: `date`"
# ---------------------------------------------------------------------
# -------------------------(END) ---------------

# Submitted batch job 41195112 - DONE
# Real time: 244676.796 sec; CPU: 232425.984 sec; Peak RSS: 209.387 GB




#-------------------------------------------------
#
#  Plot minimap result PAF >> to guide 
#  1. Name into HA412 - CHR
#  2. Reverse compliment
#
#-------------------------------------------------
module load StdEnv/2023
module load r/4.4.0
cd /home/yueyu/scratch/PRA_FASTA_almostdone

R


library(dplyr)
library(magrittr)
library(knitr)
library(ggplot2)
library(tidyr)
library(tidyverse)

# Functions

filterMum <- function(df, minl=1000, flanks=1e4){
    coord = df %>% filter(abs(re-rs)>minl) %>% group_by(qid, rid) %>%
        summarize(qsL=min(qs)-flanks, qeL=max(qe)+flanks, rs=median(rs)) %>%
        ungroup %>% arrange(desc(rs)) %>%
        mutate(qid=factor(qid, levels=unique(qid))) %>% select(-rs)
    merge(df, coord) %>% filter(qs>qsL, qe<qeL) %>%
        mutate(qid=factor(qid, levels=levels(coord$qid))) %>% select(-qsL, -qeL)
}



diagMum <- function(df){
    ## Find best qid order
    rid.o = df %>% group_by(qid, rid) %>% summarize(base=sum(abs(qe-qs)),
                                                    rs=weighted.mean(rs, abs(qe-qs))) %>%
        ungroup %>% arrange(desc(base)) %>% group_by(qid) %>% do(head(., 1)) %>%
        ungroup %>% arrange(desc(rid), desc(rs)) %>%
        mutate(qid=factor(qid, levels=unique(qid)))
    ## Find best qid strand
    major.strand = df %>% group_by(qid) %>%
        summarize(major.strand=ifelse(sum(sign(qe-qs)*abs(qe-qs))>0, '+', '-'),
                  maxQ=max(c(qe, qs)))
    merge(df, major.strand) %>% mutate(qs=ifelse(major.strand=='-', maxQ-qs, qs),
                                       qe=ifelse(major.strand=='-', maxQ-qe, qe),
                                       qid=factor(qid, levels=levels(rid.o$qid)))
}


#Load data
#  Hap 1
f <- "pra_minimap_Hap1_HA412_20250312MQ30_tab.txt"

#  Hap 2
#f <- "pra_minimap_Hap2_HA412_20250312MQ30_tab.txt"


name <- gsub(".txt", "", f)
mumgp <- read.table(f, header=F)
colnames(mumgp) <- c("rs", "re", "qs", "qe", "error", "qid", "rid", "strand")
dim(mumgp)
head(mumgp)

mumgp[c("qs", "qe")] <- t(mapply(\(a, b, c, d, e, f, g, h){
                 if(h == "-") c(d,c) else c(c,d) }, mumgp$rs, mumgp$re, mumgp$qs, mumgp$qe, mumgp$error, mumgp$qid, mumgp$rid, mumgp$strand))

mumgp.filt = filterMum(mumgp, minl=1e4)

mumgp.filt.diag = diagMum(mumgp.filt)

head(mumgp.filt.diag)
dim(mumgp.filt.diag)

#Plotting
P1 <- ggplot(mumgp.filt.diag, aes(x=rs, xend=re, y=qs, yend=qe, colour=strand)) +
  geom_segment(show.legend=FALSE, size=3) + geom_point(alpha=0.09) + theme_bw() +
  facet_grid(qid~rid, scales='free', space='free', switch='both') +
  guides(colour=guide_legend(override.aes=list(alpha=1))) +
  theme(strip.text.y=element_text(angle=120, size=10),
        strip.text.x=element_text(angle=60, size=10),
        strip.background=element_blank(),
        legend.position=c(1,-.03), legend.justification=c(1,1),
        legend.direction='horizontal',
        axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        panel.spacing=unit(0, 'cm')) +
  xlab('reference sequence') +
  ylab('assembly') +
  scale_colour_brewer(palette='Set1')


png(paste0(name, ".png"), width=2000, height= 2000)
print(P1)
dev.off()

# -- Done in R





#-------------------------------------------------
#
# (Hap 1) CHR name align to HA412
#
#-------------------------------------------------
 
imgcat pra_minimap_Hap1_HA412_20250312MQ30_tab.png
# Notes: This alignent map is very messy, try to identify chromosomes based on this
# All RED lines with ( X = Y line ) = Reversed (to be revered)



#================================
#  Prepare Step 1: Name Match
#================================


# -- Change to new folder

cd /home/yueyu/scratch/PRA_FINAL_FASTA

#-------- Change Names of scaffold change to HA412 chr - Hap 1
nano Change_Hap1_toCHR_names_PRA_2025March16.txt

h1_HiC_scaffold_10   PRA_chr01
h1_HiC_scaffold_15   PRA_chr02
h1_HiC_scaffold_1   PRA_chr03
h1_HiC_scaffold_12   PRA_chr04
h1_HiC_scaffold_11   PRA_chr05
h1_HiC_scaffold_7   PRA_chr06
h1_HiC_scaffold_4   PRA_chr07
h1_HiC_scaffold_16   PRA_chr08
h1_HiC_scaffold_2  PRA_chr09
h1_HiC_scaffold_17   PRA_chr10
h1_HiC_scaffold_14   PRA_chr11
h1_HiC_scaffold_5   PRA_chr12
h1_HiC_scaffold_13   PRA_chr13
h1_HiC_scaffold_8   PRA_chr14
h1_HiC_scaffold_3   PRA_chr15
h1_HiC_scaffold_6   PRA_chr16
h1_HiC_scaffold_9   PRA_chr17



cp /home/yueyu/scratch/PRA_GA/YAHS_Merged_Hap1_and_Hap2_PRA/reviewed_assembly/Praecox_20250310_Yue_hap1_save_17CHR.fasta ./Praecox_20250310_Yue_hap1_save_17CHR_RENAMED.fasta

Hap1fasta="Praecox_20250310_Yue_hap1_save_17CHR_RENAMED.fasta"

while read -r old_name new_name; do
    sed -i "s/\b$old_name\b/$new_name/g" $Hap1fasta
done < Change_Hap1_toCHR_names_PRA_2025March16.txt

grep ">" Praecox_20250310_Yue_hap1_save_17CHR_RENAMED.fasta

# DONE - 2025 March 16th (for PRA FINAL FASTA version)




#================================
#  Prepare Step 2: Reverse Compliment - HA412 match
#================================

module load StdEnv/2023
module load r/4.4.0
cd /home/yueyu/scratch/PRA_FINAL_FASTA


R

#------ r (start) -------- # Just run in the normal window, do not need 
library(Biostrings)


#------ Hap 1 (TEST 1)------
# Define file paths
input_fasta <- "Praecox_20250310_Yue_hap1_save_17CHR_RENAMED.fasta"   # Input FASTA file
output_fasta <- "Praecox_20250310_Yue_hap1_save_17CHR_RENAMED_RC.fasta" # Output FASTA file

# Read the sequences from the FASTA file
sequences <- readDNAStringSet(input_fasta)

sequences

# Define the headers of sequences you want to reverse complement
# Names to Reverse for PRA_HAP1 
headers_to_reverse_complement <- c("PRA_chr01","PRA_chr03","PRA_chr04","PRA_chr07","PRA_chr11","PRA_chr13","PRA_chr15","PRA_chr17")

# Reverse complement the sequences with matching headers
for (header in headers_to_reverse_complement) {
    if (header %in% names(sequences)) {
        sequences[[header]] <- reverseComplement(sequences[[header]])
    }
}

sequences

# Write the modified sequences to a new FASTA file
writeXStringSet(sequences, output_fasta)

#------ r (end) --------




#-------------------------------------------------
#
# (Hap 2) CHR name align to HA412
#
#-------------------------------------------------
 
imgcat pra_minimap_Hap2_HA412_20250312MQ30_tab.png
# Notes: This alignent map is very messy, try to identify chromosomes based on this
# All RED lines with ( X = Y line ) = Reversed (to be revered)


#================================
#  Prepare Step 1: Name Match
#================================


cd /home/yueyu/scratch/PRA_FINAL_FASTA

#-------- Change Names of scaffold change to HA412 chr - Hap 2
nano Change_Hap2_toCHR_names_PRA_2025March16.txt

h2_HiC_scaffold_10   PRA_chr01
h2_HiC_scaffold_15   PRA_chr02
h2_HiC_scaffold_1   PRA_chr03
h2_HiC_scaffold_12   PRA_chr04
h2_HiC_scaffold_11   PRA_chr05
h2_HiC_scaffold_7   PRA_chr06
h2_HiC_scaffold_4   PRA_chr07
h2_HiC_scaffold_16   PRA_chr08
h2_HiC_scaffold_2  PRA_chr09
h2_HiC_scaffold_17   PRA_chr10
h2_HiC_scaffold_14   PRA_chr11
h2_HiC_scaffold_5   PRA_chr12
h2_HiC_scaffold_13   PRA_chr13
h2_HiC_scaffold_8   PRA_chr14
h2_HiC_scaffold_3   PRA_chr15
h2_HiC_scaffold_6   PRA_chr16
h2_HiC_scaffold_9   PRA_chr17


cp /home/yueyu/scratch/PRA_GA/YAHS_Merged_Hap1_and_Hap2_PRA/reviewed_assembly/Praecox_20250310_Yue_hap2_save_17CHR.fasta ./Praecox_20250310_Yue_hap2_save_17CHR_RENAMED.fasta

Hap2fasta="Praecox_20250310_Yue_hap2_save_17CHR_RENAMED.fasta"

while read -r old_name new_name; do
    sed -i "s/\b$old_name\b/$new_name/g" $Hap2fasta
done < Change_Hap2_toCHR_names_PRA_2025March16.txt

grep ">" Praecox_20250310_Yue_hap2_save_17CHR_RENAMED.fasta

# DONE - 2025 March 16th (for PRA FINAL FASTA version)





#================================
#  Prepare Step 2: Reverse Compliment - HA412 match 
#================================

cd /home/yueyu/scratch/PRA_FINAL_FASTA

module load StdEnv/2023
module load r/4.4.0

R

#------ r (start) -------- # Just run in the normal window, do not need 
library(Biostrings)


#------ Hap 1 (TEST 1)------
# Define file paths
input_fasta <- "Praecox_20250310_Yue_hap2_save_17CHR_RENAMED.fasta"   # Input FASTA file
output_fasta <- "Praecox_20250310_Yue_hap2_save_17CHR_RENAMED_RC.fasta" # Output FASTA file

# Read the sequences from the FASTA file
sequences <- readDNAStringSet(input_fasta)

sequences

# Define the headers of sequences you want to reverse complement
# Names to Reverse for PRA_HAP1 
headers_to_reverse_complement <- c("PRA_chr02","PRA_chr05","PRA_chr06","PRA_chr08","PRA_chr09","PRA_chr10","PRA_chr12","PRA_chr14","PRA_chr16")

# Reverse complement the sequences with matching headers
for (header in headers_to_reverse_complement) {
    if (header %in% names(sequences)) {
        sequences[[header]] <- reverseComplement(sequences[[header]])
    }
}

sequences

# Write the modified sequences to a new FASTA file
writeXStringSet(sequences, output_fasta)

#------ r (end) --------





#-------------------------------------------------
#
# 2nd round of CHR name and Reverse adjustment - after GENESPACE compare to DEBILIS
#
# 2025 April 07th
#
#-------------------------------------------------
 
# Things to change
# PRA15 -> Reverse
# PRA17 -> change name to PRA16 – Reverse
# PRA16 -> change name to PRA17 


#-------------------------------------------------
#
# Hap 1  - adjust to match HA412 and DEB
#
#-------------------------------------------------

cd /home/yueyu/scratch/PRA_FINAL_FASTA


#================================
#  Prepare Step 1: Name Match
#================================


cp Praecox_20250310_Yue_hap1_save_17CHR_RENAMED_RC.fasta praecox1_20250407_before_correction.fasta

Hap1fasta="praecox1_20250407_before_correction.fasta"

# Change PRA_chr17 to PRA_chr10000
# Change PRA_chr16 to PRA_chr17
# Change PRA_chr10000 back to PRA_chr16

sed -i "s/PRA_chr17/PRA_chr10000/g" $Hap1fasta
grep ">" praecox1_20250407_before_correction.fasta

sed -i "s/PRA_chr16/PRA_chr17/g" $Hap1fasta
grep ">" praecox1_20250407_before_correction.fasta

sed -i "s/PRA_chr10000/PRA_chr16/g" $Hap1fasta
grep ">" praecox1_20250407_before_correction.fasta
# DONE - 2025 April 07 (for PRA FINAL FASTA version)




#================================
#  Prepare Step 2: Reverse Compliment
#================================

module load StdEnv/2023
module load r/4.4.0
cd /home/yueyu/scratch/PRA_FINAL_FASTA


R

#------ r (start) -------- # Just run in the normal window, do not need 
library(Biostrings)


#------ Hap 1 ------
# Define file paths
input_fasta <- "praecox1_20250407_before_correction.fasta"   # Input FASTA file
output_fasta <- "praecox1_corrected_20250407.fasta" # Output FASTA file

# Read the sequences from the FASTA file
sequences <- readDNAStringSet(input_fasta)

sequences

# Define the headers of sequences you want to reverse complement
# Names to Reverse for PRA_HAP1 
headers_to_reverse_complement <- c("PRA_chr15","PRA_chr16")

# Reverse complement the sequences with matching headers
for (header in headers_to_reverse_complement) {
    if (header %in% names(sequences)) {
        sequences[[header]] <- reverseComplement(sequences[[header]])
    }
}

sequences

# Write the modified sequences to a new FASTA file
writeXStringSet(sequences, output_fasta)

#------ r (end) --------










#-------------------------------------------------
#
# Hap 2 - adjust to match HA412 and DEB
#
#-------------------------------------------------
 
cd /home/yueyu/scratch/PRA_FINAL_FASTA

#================================
#  Prepare Step 1: Name Match
#================================

cp Praecox_20250310_Yue_hap2_save_17CHR_RENAMED_RC.fasta praecox2_20250407_before_correction.fasta

Hap2fasta="praecox2_20250407_before_correction.fasta"

# Change PRA_chr17 to PRA_chr10000
# Change PRA_chr16 to PRA_chr17
# Change PRA_chr10000 back to PRA_chr16

sed -i "s/PRA_chr17/PRA_chr10000/g" $Hap2fasta
grep ">" praecox2_20250407_before_correction.fasta

sed -i "s/PRA_chr16/PRA_chr17/g" $Hap2fasta
grep ">" praecox2_20250407_before_correction.fasta

sed -i "s/PRA_chr10000/PRA_chr16/g" $Hap2fasta
grep ">" praecox2_20250407_before_correction.fasta
# DONE - 2025 April 07 (for PRA FINAL FASTA version)




#================================
#  Prepare Step 2: Reverse Compliment
#================================

module load StdEnv/2023
module load r/4.4.0
cd /home/yueyu/scratch/PRA_FINAL_FASTA

R

#------ r (start) -------- # Just run in the normal window, do not need 
library(Biostrings)


#------ Hap 2 ------
# Define file paths
input_fasta <- "praecox2_20250407_before_correction.fasta"   # Input FASTA file
output_fasta <- "praecox2_corrected_20250407.fasta" # Output FASTA file

# Read the sequences from the FASTA file
sequences <- readDNAStringSet(input_fasta)
sequences

# Define the headers of sequences you want to reverse complement
# Names to Reverse for PRA_HAP1 
headers_to_reverse_complement <- c("PRA_chr15","PRA_chr16")

# Reverse complement the sequences with matching headers
for (header in headers_to_reverse_complement) {
    if (header %in% names(sequences)) {
        sequences[[header]] <- reverseComplement(sequences[[header]])
    }
}

sequences

# Write the modified sequences to a new FASTA file
writeXStringSet(sequences, output_fasta)

#------ r (end) --------











#####################################
#
# Part 6: Check assembly readiness: BUSCO
# 
#####################################

#=========
#  Hap 1
#=========

nano run_busco_pra_hap1.sh
#----------------- run_busco_pra_hap1.sh --------------
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=50G

module load StdEnv/2020 gcc/9.3.0 python/3.10 augustus/3.5.0 hmmer/3.3.2 blast+/2.13.0 metaeuk/6 prodigal/2.6.3 r/4.3.1 bbmap/38.86
source ~/busco_env/bin/activate

cd /home/yueyu/scratch/PRA_GA/YAHS_Merged_Hap1_and_Hap2_PRA/reviewed_assembly

INPUT="Praecox_20250310_Yue_hap1_save_17CHR.fasta"
OUTPUT="BUSCO_check/PRA_Hap1_BUSCO_2025Mar12"

busco --offline -m genome -c 10 -i $INPUT -o $OUTPUT -l /home/yueyu/scratch/GA/Hdeb2414_merged_Hap1_andHap2/reviewed_assembly/busco_downloads/lineages/eudicots_odb10

#----------------- run_busco_pra_hap1.sh (END) --------------

#=========
#  Hap 2
#=========

nano run_busco_pra_hap2.sh
#----------------- run_busco_pra_hap2.sh --------------
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=50G

module load StdEnv/2020 gcc/9.3.0 python/3.10 augustus/3.5.0 hmmer/3.3.2 blast+/2.13.0 metaeuk/6 prodigal/2.6.3 r/4.3.1 bbmap/38.86
source ~/busco_env/bin/activate

cd /home/yueyu/scratch/PRA_GA/YAHS_Merged_Hap1_and_Hap2_PRA/reviewed_assembly

INPUT="Praecox_20250310_Yue_hap2_save_17CHR.fasta"
OUTPUT="BUSCO_check/PRA_Hap2_BUSCO_2025Mar12"

busco --offline -m genome -c 10 -i $INPUT -o $OUTPUT -l /home/yueyu/scratch/GA/Hdeb2414_merged_Hap1_andHap2/reviewed_assembly/busco_downloads/lineages/eudicots_odb10

#----------------- run_busco_pra_hap2.sh (END) --------------










#####################################
#
# Part 7: Final Minimap
# 
#####################################

#=====================
#
#  FINAL FASTA files (2025 March 16th)
#
#=====================
# Renamed to HA412 and reversed 

cd /AsteroidScratch/Syri_YueYu_tmp

mkdir BAM_2025April_PRA

pra_hap1="/AsteroidScratch/Annotation_YueYu/LIFTOFF/data/praecox1.fasta"
pra_hap2="/AsteroidScratch/Annotation_YueYu/LIFTOFF/data/praecox2.fasta"

# do NOT use default minimap2 on Sundance, install the latest: /AsteroidScratch/Syri_YueYu_tmp/minimap2-2.28/minimap2

/AsteroidScratch/Syri_YueYu_tmp/minimap2-2.28/minimap2 -ax asm5 --eqx $pra_hap1 $pra_hap2 > BAM_2025April_PRA/pra_minimap_final_20250407.sam
# asm5: average divergence is << 5%.
# eqx: Output =/X CIGAR operators for sequence match/mismatch.

cd BAM_2025April_PRA
samtools view -b pra_minimap_final_20250407.sam > pra_minimap_final_20250407.bam

# ---------------------------------------------------------------------
echo "Job finished with exit code $? at: `date`"
# ---------------------------------------------------------------------


#=====================
#
#  Syri + plotSR
#
#=====================


cd /home/yueyu/scratch/PRA_FINAL_FASTA/BAM_2025April_PRA

module load StdEnv/2020 gcc python/3.9 igraph
source ~/ENV/bin/activate

#CHANGE THIS!dd!
pra_hap1="/home/yueyu/scratch/PRA_FINAL_FASTA/praecox1_corrected_20250407.fasta"
pra_hap2="/home/yueyu/scratch/PRA_FINAL_FASTA/praecox2_corrected_20250407.fasta"

syri -c pra_minimap_final_20250407.bam -r $pra_hap1 -q $pra_hap2 -F B --prefix final_PRA_20250407_


nano genomes_pra_20250407.txt
# ------- genomes_pra_20250407.txt ------
#file   name    tags
/home/yueyu/scratch/PRA_FINAL_FASTA/praecox1_corrected_20250407.fasta   PRA_Hap_1       lw:1.5
/home/yueyu/scratch/PRA_FINAL_FASTA/praecox2_corrected_20250407.fasta   PRA_Hap_2       lw:1.5
# ------- genomes_pra_20250407.txt ------
# Must be TAB seperated, not space 

plotsr --sr final_PRA_20250407_syri.out --genomes genomes_pra_20250407.txt -H 8 -W 5

mv plotsr.pdf 2025April11_plotsr_DEB_final.pdf

cd /Users/yueyu/Desktop/REF_DEB_PRA/PRA
scp yueyu@narval.computecanada.ca:"/home/yueyu/scratch/PRA_FINAL_FASTA/BAM_2025April_PRA/FINAL_PRA_SYRI_RESULT_2025April11.pdf" .










#####################################
#
# Part 8: Final Minimap between DEB hap2 - PRA hap2
# 
#####################################

deb_hap2="/AsteroidScratch/Syri_YueYu_tmp/BAM_2026March_DEBhap2_PRAhap2/debilis2_renamedCHR.fasta"
pra_hap2="/AsteroidScratch/Syri_YueYu_tmp/BAM_2026March_DEBhap2_PRAhap2/praecox2_renamedCHR.fasta"

cd BAM_2026March_DEBhap2_PRAhap2
/AsteroidScratch/Syri_YueYu_tmp/minimap2-2.28/minimap2 -ax asm5 --eqx $pra_hap2 $deb_hap2 | samtools sort -o prahap2_debhap2_minimap_final_20260323.bam
# asm5: average divergence is << 5% (DOUBLE CHECK!!!!!)
# eqx: Output =/X CIGAR operators for sequence match/mismatch.
# -ax asm5 Full genome/assembly alignment

# ---------------------------------------------------------------------
echo "Job finished with exit code $? at: `date`"
# ---------------------------------------------------------------------
# ------------------

cd /AsteroidScratch/Syri_YueYu_tmp/BAM_2026March_DEBhap2_PRAhap2

pra_hap2="praecox2_renamedCHR.fasta"
deb_hap2="debilis2_renamedCHR.fasta"

syri -c prahap2_debhap2_minimap_final_20260323.bam -r $pra_hap2 -q $deb_hap2 -F B --prefix final_PRA_hap2_DEB_hap2_2026March23_

# --  plotSR 
nano genomes_2026march23.txt
# ------- genomes_2026march23.txt ------
#file   name    tags
/AsteroidScratch/Syri_YueYu_tmp/BAM_2026March_DEBhap2_PRAhap2/praecox2_renamedCHR.fasta   PRA_Hap_2       lw:1.5
/AsteroidScratch/Syri_YueYu_tmp/BAM_2026March_DEBhap2_PRAhap2/debilis2_renamedCHR.fasta   DEB_Hap_2       lw:1.5
# ------- genomes_2026march23.txt ------
# Must be TAB seperated, not space 


plotsr --sr final_PRA_hap2_DEB_hap2_2026March23_syri.out --genomes genomes_2026march23.txt -H 8 -W 5
mv plotsr.pdf 2026March23_plotsr_final.pdf


# END