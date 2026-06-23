#####################################
# Project: debilis_praecox_project
#
# Code 3: Genome Assembly & Annotation for DEB
#
# by: Yue Yu
#
#####################################


# Including all the following part:
# -- Part 1: Juicer
# -- Part 2: 3D-DNA
# -- Part 3: JuiceBox manual curation
# -- Part 4: Review Assembly -> FASTA
# -- Part 5: Reverse Compliment
# -- Part 6: Minimap + syri + plotSR
# -- Part 7: Merged 1+2 -> Juicer again
# -- Part 8: Reorder and rename CHR based on HA412
# -- Part 9: Check assembly readiness: BUSCO
# -- Part 10: Final Minimap

# -- Part 11: GENESPACE
# -- Part 12: Circo plot 


 
#####################################
#
# Part 1: Juicer
# 
#####################################

#Download Juicer  (correct version!!)
cd /home/yueyu/scratch/GA
git clone https://github.com/rmdickson/juicer.git

cd /home/yueyu/scratch/GA/Hdeb2414_HIC_oldhifiasm.hic.hap1
ln -s /home/yueyu/scratch/GA/juicer/CPU/ scripts
cd scripts/common
wget https://hicfiles.tc4ga.com/public/juicer/juicer_tools.1.9.9_jcuda.0.8.jar
ln -s juicer_tools.1.9.9_jcuda.0.8.jar  juicer_tools.jar


cd /home/yueyu/scratch/GA/Hdeb2414_HIC_oldhifiasm.hic.hap2
ln -s /home/yueyu/scratch/GA/juicer/CPU/ scripts
#run juicier with out GPUS and more cores


#-----------
#   Hap 1
#-----------
cd /home/yueyu/scratch/GA/Hdeb2414_HIC_oldhifiasm.hic.hap1
ln -s /home/yueyu/scratch/GA/fastq 

nano Deb_Hap01_juicer.sh

#NEED TO BE IN THE CURRENT FOLDER TO RUN THE FOLLOWING CODE
-rwxr-x---. 1 yueyu yueyu   585 Jul 11 17:13 Deb_Hap01_juicer.sh
lrwxrwxrwx. 1 yueyu yueyu    28 Jul 11 17:10 fastq -> /home/yueyu/scratch/GA/fastq
drwxr-x---. 2 yueyu yueyu 25600 Jul 11 15:21 references
drwxr-x---. 2 yueyu yueyu 25600 Jul 11 15:21 restriction_sites
lrwxrwxrwx. 1 yueyu yueyu    34 Jul 11 17:25 scripts -> /home/yueyu/scratch/GA/juicer/CPU/

#----------------- Deb_Hap01_juicer.sh --------------

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=def-rieseber
#SBATCH --time=3-0
#SBATCH --ntasks-per-node=32
#SBATCH --mem=149G

cd /home/yueyu/scratch/GA/Hdeb2414_HIC_oldhifiasm.hic.hap1

module load StdEnv/2020 bwa/0.7.17 java/17.0.2 samtools/1.15.1

bash scripts/juicer.sh -D $PWD \
					   -g contigassembly \
					   -s DpnII \
					   -p restriction_sites/Hdeb2414_HIC_oldhifiasm.hic.hap1_DpnII.chrom.sizes \
					   -y restriction_sites/Hdeb2414_HIC_oldhifiasm.hic.hap1_DpnII.txt \
					   -z references/Hdeb2414_HIC_oldhifiasm.hic.hap1.fasta \
					   -t 32 \
					   -S early

#----------------- Deb_Hap01_juicer.sh (END) --------------

#Submitted batch job 31457074 (running 2024 July 11th 2:30pm)

#DONE


sbatch Deb_Hap01_juicer_2024Nov.sh
Submitted batch job 36925069
# submitted Nov 15th 9:02pm


#-----------
#   Hap 2
#-----------
cd /home/yueyu/scratch/GA/Hdeb2414_HIC_oldhifiasm.hic.hap2
ln -s /home/yueyu/scratch/GA/fastq 

nano Deb_Hap02_juicer.sh

#----------------- Deb_Hap02_juicer.sh --------------

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=def-rieseber
#SBATCH --time=3-0
#SBATCH --ntasks-per-node=32
#SBATCH --mem=149G

cd /home/yueyu/scratch/GA/Hdeb2414_HIC_oldhifiasm.hic.hap2

module load StdEnv/2020 bwa/0.7.17 java/17.0.2 samtools/1.15.1

bash scripts/juicer.sh -D $PWD \
					   -g contigassembly \
					   -s DpnII \
					   -p restriction_sites/Hdeb2414_HIC_oldhifiasm.hic.hap2_DpnII.chrom.sizes \
					   -y restriction_sites/Hdeb2414_HIC_oldhifiasm.hic.hap2_DpnII.txt \
					   -z references/Hdeb2414_HIC_oldhifiasm.hic.hap2.fasta \
					   -t 32 \
					   -S early

#----------------- Deb_Hap02_juicer.sh (END) --------------





 
#####################################
#
# Part 2: 3D-DNA
# 
#####################################

#-----------
#   Hap 1
#-----------
cd /home/yueyu/scratch/GA/Hdeb2414_HIC_oldhifiasm.hic.hap1
mkdir Hap1_3D_DNA

nano Deb_Hap01_3D_DNA.sh
#----------------- Deb_Hap01_3D_DNA.sh --------------
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=7-0
#SBATCH --cpus-per-task=15
#SBATCH --mem=100G

cd /home/yueyu/scratch/GA/Hdeb2414_HIC_oldhifiasm.hic.hap1

module load StdEnv/2020 python/3.10.2 java/17.0.2 lastz/1.04.03

export PATH="/Hap1_3D_DNA/3d-dna:$PATH"

source /home/yueyu/scratch/GA/3d-dna/3ddna/bin/activate

/home/yueyu/scratch/GA/3d-dna/run-asm-pipeline.sh -r 0 references/Hdeb2414_HIC_oldhifiasm.hic.hap1.fasta aligned/merged_nodups.txt

deactivate
#----------------- Deb_Hap01_3D_DNA.sh (END) --------------



#-----------
#   Hap 2
#-----------
cd /home/yueyu/scratch/GA/Hdeb2414_HIC_oldhifiasm.hic.hap2
mkdir Hap2_3D_DNA

nano Deb_Hap02_3D_DNA.sh

#----------------- Deb_Hap02_3D_DNA.sh --------------

#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=7-0
#SBATCH --cpus-per-task=15
#SBATCH --mem=100G

cd /home/yueyu/scratch/GA/Hdeb2414_HIC_oldhifiasm.hic.hap2

module load StdEnv/2020 python/3.10.2 java/17.0.2 lastz/1.04.03

export PATH="/Hap2_3D_DNA/3d-dna:$PATH"

source /home/yueyu/scratch/GA/3d-dna/3ddna/bin/activate

/home/yueyu/scratch/GA/3d-dna/run-asm-pipeline.sh -r 0 references/Hdeb2414_HIC_oldhifiasm.hic.hap2.fasta aligned/merged_nodups.txt

deactivate
#----------------- Deb_Hap02_3D_DNA.sh (END) --------------



 
#####################################
#
# Part 3: JuiceBox manual curation
# 
#####################################

#=====================
# Require file 1 .hic
# Require file 2 .assembly (the one ends with 0.assembly)
#=====================

# ---- A few things to check using JBAT (on large scale)
# 1. misjoins
# 2. translocations
# 3. inversions
# 4. fix chromosome boundaries

# ---- A few things to do using JBAT (on medium scale)
# 1. check if blocks align with chromosome number expected 
# 2. remove all scaffolds on the very beginning of genome
# 3. remove all scaffolds in between chromosomes
# 4. Define chromosome boundaries





#####################################
#
# Part 4: Review Assembly -> FASTA
# 
#####################################
#1. MIX_ASSEM <- args[1] #Reviewed assembly file
#2. PREFIX <- args[2]    #a random name to add to output
#3. MIX_FASTA <- args[3] #FASTA used to run Juicer (in the reference file)
#4. SAVE_DIR <- args[4]

#----- Hap 1 
MIX_ASSEM="Hdeb2414_HIC_oldhifiasm.hic.hap1/JBAT_reviewed_assembly/Hdeb2414_HIC_oldhifiasm.hic.hap1.0.review.assembly"
PREFIX="debilis_hap1"
MIX_FASTA="Hdeb2414_HIC_oldhifiasm.hic.hap1/references/Hdeb2414_HIC_oldhifiasm.hic.hap1.fasta"
SAVE_DIR="/home/yueyu/scratch/GA/Hdeb2414_HIC_oldhifiasm.hic.hap1/JBAT_reviewed_FASTA_Hap1"

screen -r run
module load StdEnv/2023
module load r/4.4.0
cd /home/yueyu/scratch/GA

Rscript ASM_TO_FASTA_Yue_Version_runfor_DEB_HAP1.R $MIX_ASSEM $PREFIX $MIX_FASTA $SAVE_DIR



#----- Hap 2 
MIX_ASSEM="Hdeb2414_HIC_oldhifiasm.hic.hap2/JBAT_reviewed_assembly/Hdeb2414_HIC_oldhifiasm.hic.hap2.0.review.assembly"
PREFIX="debilis_hap2"
MIX_FASTA="Hdeb2414_HIC_oldhifiasm.hic.hap2/references/Hdeb2414_HIC_oldhifiasm.hic.hap2.fasta"
SAVE_DIR="/home/yueyu/scratch/GA/Hdeb2414_HIC_oldhifiasm.hic.hap2/JBAT_reviewed_FASTA_Hap2"

screen -r run
module load StdEnv/2023
module load r/4.4.01
cd /home/yueyu/scratch/GA

Rscript ASM_TO_FASTA_Yue_Version_runfor_DEB_HAP2.R $MIX_ASSEM $PREFIX $MIX_FASTA $SAVE_DIR



#####################################
#
# Part 5: Reverse Compliment
#
# Part 6: Minimap + syri + plotSR
# 
# Part 7: Merged 1+2 -> Juicer again
# 
# Part 8: Reorder and rename CHR based on HA412
# 
# Part 9: Check assembly readiness: BUSCO
# 
# Part 10: Final Minimap
# 
#####################################
# Same code as in Code_2


#####################################
#
# Part 11: GENESPACE
# 
#####################################

# Liftoff to HA412


cd /AsteroidScratch/Annotation_YueYu/LIFTOFF/data

# -- Step 2: reference HA412 FASTA+GFF3

scp yueyu@narval.computecanada.ca:"/home/yueyu/scratch/HA412/Ha412HOv2.0-20181130.*" .
scp yueyu@narval.computecanada.ca:"/home/yueyu/scratch/HA412/HAN412_Final.sorted.*" .

mv Ha412HOv2.0-20181130.fasta HA412.fasta
mv Ha412HOv2.0-20181130.fasta.fai HA412.fasta.fai
mv HAN412_Final.sorted.gff3 HA412.gff3
mv HAN412_Final.sorted.pep HA412.pep

-rwxr-x--- 1 yue_yu 109M Mar 19 10:28 HA412.gff3
-rwxr-x--- 1 yue_yu  17M Mar 19 10:29 HA412.pep
-rwxr-x--- 1 yue_yu 3.1G Mar 19 10:31 HA412.fasta
-rwxr-x--- 1 yue_yu 986K Mar 19 10:31 HA412.fasta.fai


# -- Step 3.1: Query FASTA - Debilis
cd /AsteroidScratch/Annotation_YueYu/LIFTOFF/data

cp /AsteroidScratch/Annotation_YueYu/Yue_Final_DEB_FASTA_2025Feb07/debilis_hap?_toCHR_2024Dec13_RC.fasta .

mv debilis_hap1_toCHR_2024Dec13_RC.fasta debilis1.fasta
mv debilis_hap2_toCHR_2024Dec13_RC.fasta debilis2.fasta


# -- Step 3.2: Query FASTA - Praecox
cd /AsteroidScratch/Annotation_YueYu/LIFTOFF/data

scp yueyu@narval.computecanada.ca:"/home/yueyu/scratch/PRA_FINAL_FASTA/praecox?_corrected_20250407.fasta" .

mv praecox1_corrected_20250407.fasta praecox1.fasta
mv praecox2_corrected_20250407.fasta praecox2.fasta


# -- Step 4.1: LIFTOFF - Debilis - Hap 1  
cd /AsteroidScratch/Annotation_YueYu/LIFTOFF
source /build/Cellar/miniconda3/bin/activate liftoff

SPP_query=debilis1
SPP_ref=HA412

# Make folder first before run liftoff
# mkdir output

liftoff \
-g ./data/${SPP_ref}.gff3 -p 40 -o ./output/${SPP_query}_${SPP_ref}.gff3 \
-dir ${SPP_query}_${SPP_ref}_intermed -u ${SPP_query}_${SPP_ref}_unmapped_features.txt \
-infer_genes -copies -a 0.95 -s 0.95 -d 5.0 -flank 0.8 -polish \
./data/${SPP_query}.fasta ./data/${SPP_ref}.fasta


# -- Step 4.2: LIFTOFF - Debilis - Hap 2 
SPP_query=debilis2
SPP_ref=HA412

# -- Step 4.3: LIFTOFF - Praecox - Hap 1
SPP_query=praecox1
SPP_ref=HA412


# -- Step 4.4: LIFTOFF - Praecox - Hap 2 
SPP_query=praecox2
SPP_ref=HA412



# -- Step 5: RUN gff3 + FASTA --> protein.fa
cd /AsteroidScratch/Annotation_YueYu/LIFTOFF/output

SPP_query=debilis1
SPP_ref=HA412

/AsteroidScratch/Annotation_YueYu/software/gffread/gffread -g ../data/${SPP_query}.fasta \
 ${SPP_query}_${SPP_ref}.gff3_polished \
 -y ${SPP_query}_${SPP_ref}_proteins.pep

grep ">" debilis1_HA412_proteins.pep | wc -l
# 40,103


SPP_query=debilis2
SPP_ref=HA412

/AsteroidScratch/Annotation_YueYu/software/gffread/gffread -g ../data/${SPP_query}.fasta \
 ${SPP_query}_${SPP_ref}.gff3_polished \
 -y ${SPP_query}_${SPP_ref}_proteins.pep

grep ">" debilis2_HA412_proteins.pep | wc -l
# 40,075




SPP_query=praecox1
SPP_ref=HA412

/AsteroidScratch/Annotation_YueYu/software/gffread/gffread -g ../data/${SPP_query}.fasta \
 ${SPP_query}_${SPP_ref}.gff3_polished \
 -y ${SPP_query}_${SPP_ref}_proteins.pep

grep ">" praecox1_HA412_proteins.pep | wc -l
# 40,167



SPP_query=praecox2
SPP_ref=HA412

/AsteroidScratch/Annotation_YueYu/software/gffread/gffread -g ../data/${SPP_query}.fasta \
 ${SPP_query}_${SPP_ref}.gff3_polished \
 -y ${SPP_query}_${SPP_ref}_proteins.pep

grep ">" praecox2_HA412_proteins.pep | wc -l
# 40,155






# --  Step 6: Prepare BED and PEP

# ===================
#
#  HA412  -  44,544 genes
#
# ===================


# -- Step 2: PEP file
# check PEP file match gene names in BED (HA412)
export PATH=/AsteroidScratch/Annotation_YueYu/software:$PATH

seqkit grep -rvp "Ha412HOChr00" HA412.pep > HA412_forgenespace.pep
# seqkit Version: 2.10.0

# 44,544 PEP



# --  Step 7:  BED file
# convert GFF3 file -> BED (HA412)
cd /AsteroidScratch/Annotation_YueYu/GENESPACE

SP="HA412"
awk '$3 == "gene" {split($9, a, ";"); split(a[1], b, "="); print $1, $4, $5, b[2]}' OFS="\t" ${SP}.gff3 | grep -v "Ha412HOChr00" > ${SP}.bed

wc -l HA412.bed
# 46,223 GENES


# -- Subset BED based on PEP file
# -- Question: Why does PEP have less entry than GFF3 gene lines???

grep "^>" HA412_forgenespace.pep | sed 's/^>//' > HA412_overlap_genes.txt
wc -l HA412_overlap_genes.txt
# 44,544 overlap lines (correct)


awk 'NR==FNR {names[$1]; next} $4 in names' HA412_overlap_genes.txt HA412.bed > HA412_forgenespace.bed
# 44,544 (subset successful)

#rm HA412.bed HA412.gff3 HA412.pep HA412_overlap_genes.txt





# ===================
#
#   Debilis Hap 1  - 39,741 genes
#
# ===================

cp /AsteroidScratch/Annotation_YueYu/LIFTOFF/output/*HA412.gff3_polished .
cp /AsteroidScratch/Annotation_YueYu/LIFTOFF/output/*HA412_proteins.pep .

SP="debilis1_HA412"

# -------------------- Repeat code -----------------

# My BED file:
awk '$3 == "gene" {split($9, a, ";"); split(a[1], b, "="); split(b[2], c, ":"); print $1, $4, $5, c[2]}' OFS="\t" ${SP}.gff3_polished | grep -v "Ha412HOChr00" > ${SP}.bed
wc -l ${SP}.bed
# 42,783 debilis1_HA412.bed

# My PEP file:
sed 's/>mRNA:/>/g' ${SP}_proteins.pep > ${SP}_mid.pep
seqkit grep -rvp "Ha412HOChr00" ${SP}_mid.pep > ${SP}_forgenespace.pep

grep ">" ${SP}_forgenespace.pep | wc -l
# 39,741 debilis1_HA412


# Subset BED to PEP names
grep "^>" ${SP}_forgenespace.pep | sed 's/^>//' > ${SP}_GENE_NAMES.txt
wc -l ${SP}_GENE_NAMES.txt

awk 'NR==FNR {names[$1]; next} $4 in names' ${SP}_GENE_NAMES.txt ${SP}.bed > ${SP}_forgenespace.bed
wc -l ${SP}_forgenespace.bed


# remove original file
rm ${SP}.gff3_polished ${SP}_proteins.pep

# remove intermid files
rm ${SP}.bed ${SP}_mid.pep ${SP}_GENE_NAMES.txt

# -------------------- Repeat code (END) -----------------


# ===================
#
#   Debilis Hap 2  - 39,714 genes
#
# ===================
SP="debilis2_HA412"

# Repeat code above

# ===================
#
#   Praecox Hap 1  - 39,797 genes 
#
# ===================
SP="praecox1_HA412"

# Repeat code above

# ===================
#
#   Praecox Hap 2 - 39,796 genes
#
# ===================
SP="praecox2_HA412"

# Repeat code above


# -- Step 8 : Dir format (manual parsing the files)
cd /AsteroidScratch/Annotation_YueYu/GENESPACE

mkdir peptide
mkdir bed

mv *.bed ./bed
mv *.pep ./peptide

change all .pep end to .fa

# -- Note: Make sure the corresponding file name match in ./bed and ./peptide/

ls -thor ./bed
total 9.9M
-rw-rw-r-- 1 yue_yu 2.3M Mar 21 15:50 HA412_forgenespace.bed
-rw-rw-r-- 1 yue_yu 1.9M Mar 21 16:16 debilis1_HA412_forgenespace.bed
-rw-rw-r-- 1 yue_yu 1.9M Mar 21 16:29 debilis2_HA412_forgenespace.bed
-rw-rw-r-- 1 yue_yu 1.9M Mar 21 16:30 praecox1_HA412_forgenespace.bed
-rw-rw-r-- 1 yue_yu 1.9M Mar 21 16:31 praecox2_HA412_forgenespace.bed

ls -thor ./peptide
total 77M
-rw-rw-r-- 1 yue_yu 17M Mar 21 15:37 HA412_forgenespace.fa
-rw-rw-r-- 1 yue_yu 15M Mar 21 16:14 debilis1_HA412_forgenespace.fa
-rw-rw-r-- 1 yue_yu 15M Mar 21 16:28 debilis2_HA412_forgenespace.fa
-rw-rw-r-- 1 yue_yu 15M Mar 21 16:30 praecox1_HA412_forgenespace.fa
-rw-rw-r-- 1 yue_yu 15M Mar 21 16:31 praecox2_HA412_forgenespace.fa



# -- Step 9: Run GENESPACE

screen -S genespace
screen -r genespace

# -- activate environment for running GENESPACE
source /build/Cellar/miniconda3/bin/activate genespace

# -- Find location of MCScanX
which MCScanX
# /build/Cellar/miniconda3/envs/genespace/bin


nano genespace_run.R

R

# ---------------------- genespace_run.R ----------------------

###############################################
# -- change paths to those valid on your system
# genomeRepo 
wd <- "/AsteroidScratch/Annotation_YueYu/GENESPACE"
path2mcscanx <- "/build/Cellar/miniconda3/envs/genespace/bin"
###############################################

library(devtools)
library(GENESPACE)

#genomes2run <- c("HA412_forgenespace", "debilis1_HA412_forgenespace", "debilis2_HA412_forgenespace", "praecox1_HA412_forgenespace.pep", "praecox2_HA412_forgenespace.pep")

# -- initalize the run and QC the inputs (fast)
gpar <- init_genespace(
  wd = wd, nCores = 30,
  path2mcscanx = path2mcscanx)

# -- accomplish the run (more MEM required)
out <- run_genespace(gpar)

save.image("genespace_out.RData")

# RUNNING - 2025 April 7th - 2:13pm
# DONE - 2025 April 7th (ran for ~20 min)

# ---------------------- genespace_run.R (END) ----------------------

# If running in HPC
# Rscript genespace_run.R
# rm -r dotplots genespace_run.R orthofinder pangenes results riparian syntenicHits tmp


# 2025 Mrach 25th - 2pm
# ERROR 1: All protien.fa files for debilis and praecox have unexpected "." in the fasta file, can not be interpreted by orthofinder
# Possible cause: could be due to liftoff program when lifting genes from HA412, created gaps represented by . instead of -
# Solution: change all "." to "-" 

cd /AsteroidScratch/Annotation_YueYu/GENESPACE/peptide

for file in *.fa; do
  echo "Processing $file..."
  sed -i 's/\./-/g' "$file"
done


# -- Outputs from successful run

############################
GENESPACE run complete!  All results are stored in
/AsteroidScratch/Annotation_YueYu/GENESPACE in the following
subdirectories:
        syntenic block dotplots: /dotplots (...synHits.pdf)
        annotated blast files  : /syntenicHits
        annotated/combined bed : /results/combBed.txt
        syntenic block coords. : /results/blkCoords.txt
        syn. blk. by ref genome: /riparian/refPhasedBlkCoords.txt
        pan-genome annotations : /pangenes (...pangenes.txt.gz)
        riparian plots         : /riparian
        genespace param. list  : /results/gsParams.rda
############################


# -- Step 6 : check global synteny
cd riparian

scp *.pdf to local server



# -- Step 9 : CUSTOM - Visualize global synteny

# Resourses 1: https://htmlpreview.github.io/?https://github.com/jtlovell/tutorials/blob/main/genespaceGuide.html
# Resourses 2: https://htmlpreview.github.io/?https://github.com/jtlovell/tutorials/blob/main/riparianGuide.html

# -- Step 10: Plot only PRA and DEB

# - scale: genes
ripDat <- plot_riparian(
  gsParam = out, 
  refGenome = "debilis2_HA412_forgenespace",
  genomeIDs = c("debilis2_HA412_forgenespace", "debilis1_HA412_forgenespace", "praecox1_HA412_forgenespace", "praecox2_HA412_forgenespace"), 
  forceRecalcBlocks = FALSE)

pdf(file = "DEB_PRA_2025April07.pdf",width = 6, height = 4) 
#Higher resolution changes make here!
print(ripDat)
dev.off()


# - scale: physical distance
ripDat_pysical <- plot_riparian(
  gsParam = out, 
  refGenome = "debilis2_HA412_forgenespace",
  genomeIDs = c("debilis2_HA412_forgenespace", "debilis1_HA412_forgenespace", "praecox1_HA412_forgenespace", "praecox2_HA412_forgenespace"), 
  useOrder = FALSE,
  forceRecalcBlocks = FALSE)


pdf(file = "DEB_PRA_physical_dis_2025April07.pdf",width = 6, height = 4) 
#Higher resolution changes make here!
print(ripDat_pysical)
dev.off()


# -- Step 11: Highlight Inversions

ripDat_INV <- plot_riparian(
  gsParam = out, 
  refGenome = "debilis2_HA412_forgenespace",
  genomeIDs = c("debilis2_HA412_forgenespace", "debilis1_HA412_forgenespace", "praecox1_HA412_forgenespace", "praecox2_HA412_forgenespace"), 
  useOrder = FALSE,
  chrLabFontSize = 6,
  inversionColor = "deeppink")

pdf(file = "DEB_INVERTEDpra_INVERSION_HIGHLIGHT_2025April07.pdf",width = 6, height = 4) 
#Higher resolution changes make here!
print(ripDat_INV)
dev.off()




# -- Step 11: Extract all SV (INVERSION) between DEB and PRA

# E.4 in: https://htmlpreview.github.io/?https://github.com/jtlovell/tutorials/blob/main/genespaceGuide.html

cd /AsteroidScratch/Annotation_YueYu/GENESPACE/results

# -------- Subset, inversion with "-" orientation

# ----------
# Between two REF genome
# ----------

# REF == "DEB2" - USED AS REFERENCE GENOME FOR GBS DATA
# QUR == "PRA2" - USED AS REFERENCE GENOME FOR GBS DATA
awk -F',' '($1 == "debilis2_HA412_forgenespace" && $2 == "praecox2_HA412_forgenespace" && $20 == "-") {print $1, $2, $3, $4, $6, $7, $20, $21, $22}' $INPUT | sed 's/,/\t/g' > subset_DEB2_PRA2_syntenicBlock_coordinates.txt


# ----------
# Within H1 and H2 - DEB
# ----------
# REF == "DEB2"
# QUR == "DEB1"
awk -F',' '($1 == "debilis2_HA412_forgenespace" && $2 == "debilis1_HA412_forgenespace" && $20 == "-") {print $1, $2, $3, $4, $6, $7, $20, $21, $22}' $INPUT | sed 's/,/\t/g' > subset_DEB2_DEB1_syntenicBlock_coordinates.txt



# ----------
# Within H1 and H2 - PRA
# ----------
# REF == "PRA2"
# QUR == "PRA1"
awk -F',' '($1 == "praecox1_HA412_forgenespace" && $2 == "praecox2_HA412_forgenespace" && $20 == "-") {print $1, $2, $3, $4, $6, $7, $20, $21, $22}' $INPUT | sed 's/,/\t/g' > subset_PRA2_PRA1_syntenicBlock_coordinates.txt





#####################################
#
# Part 12: Circo plot
# 
#####################################


# ---- Overview: Things to plor in Circo plot -------

# a: minimap + Syri --> output to synteny ribbons (Syri - Yue)
# b: Gene density (Ji annotation) 1Mbp window 
# c: TE density (Ji anotation)
# d: Predicted centromere ( RepeatOBserver - Yue ) 

# ---------------------------------------------------

# -----------
#
# Color extracted for circos (18 options)
#
# -----------
#c1a6e5
#d9b456
#e0eac4
#d67536
#725535
#747d4a
#579192
#b6c9bc
#d5d58d
#c3bb81
#8ec2b1
#823635
#d9551e
#dc8526
#e5d935
#8a9a7a
#b5c499
#6fbe80




# -----------
#
# a: minimap + Syri
#
# -----------

# ======
#  DEB 
# ======

cd /home/yueyu/scratch/circo_plot/syri

nano minimap.align.code.sh
#--------- minimap.align.code.sh  --------- 
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=1-0
#SBATCH --cpus-per-task=1
#SBATCH --mem=300G

module load StdEnv/2023
module load minimap2/2.28

module load gcc/12.3
module load samtools/1.20

# IMPORTANT CHANGE!
cd /home/yueyu/scratch/circo_plot/syri

deb_hap1="/home/yueyu/scratch/Repeatobserver_2026Mar/DEB/hap1/debilis_hap1_toCHR_2024Dec13_RC.fasta"
deb_hap2="/home/yueyu/scratch/Repeatobserver_2026Mar/DEB/hap2/debilis_hap2_toCHR_2024Dec13_RC.fasta"

minimap2 -ax asm5 --eqx $deb_hap1 $deb_hap2 > deb_minimap_20260303.sam

samtools view -b deb_minimap_20260303.sam > deb_minimap_20260303.bam

# ---------------------------------------------------------------------
echo "Job finished with exit code $? at: `date`"
# ---------------------------------------------------------------------

#--------- minimap.align.code.sh (END) --------- 


cd /home/yueyu/scratch/circo_plot/syri
module load StdEnv/2020 gcc python/3.9 igraph
source ~/ENV/bin/activate

#CHANGE THIS!!
deb_hap1="/home/yueyu/scratch/Repeatobserver_2026Mar/DEB/hap1/debilis_hap1_toCHR_2024Dec13_RC.fasta"
deb_hap2="/home/yueyu/scratch/Repeatobserver_2026Mar/DEB/hap2/debilis_hap2_toCHR_2024Dec13_RC.fasta"

syri -c deb_minimap_20260303.bam -r $deb_hap1 -q $deb_hap2 -F B -k --prefix deb_20260303
# Completed 20206 March 04

# -- Summary print
cat deb_20260303syri.summary

# -- (OPTONAL) Subset column 11 for Annotation type 
# -- SYN: Syntenic region (1189)
awk -F'\t' '$11 == "SYN"' deb_20260303syri.out > deb_20260303syri_SYN.out
# -- INV: Inverted region (93)
awk -F'\t' '$11 == "INV"' deb_20260303syri.out > deb_20260303syri_INV.out
# -- TRANS: Translocated region (592) 
awk -F'\t' '$11 == "TRANS"' deb_20260303syri.out > deb_20260303syri_TRANS.out
# -- INVTR: Inverted translocated region (582) + 592 above = 1174 
awk -F'\t' '$11 == "INVTR"' deb_20260303syri.out > deb_20260303syri_INVTR.out
# -- DUP: Duplicated region (1894)
awk -F'\t' '$11 == "DUP"' deb_20260303syri.out > deb_20260303syri_DUP.out
# -- INVDP: Inverted duplicated region (1744) + 1894 = 3638 = 3424 + 214
awk -F'\t' '$11 == "INVDP"' deb_20260303syri.out > deb_20260303syri_INVDP.out


# Extract inversion > 500Kbp
awk -F'\t' '$11 == "INV"' deb_20260303syri.out > INVERSION_deb_20260303syri.out
# Move to excel -> calc START-END -> Extract Inversion above 500,000bp






# ======
#  PRA 
# ======

cd /home/yueyu/scratch/circo_plot_PRA/syri

nano minimap.align.code.pra.sh
#--------- minimap.align.code.pra.sh  --------- 
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=30:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=300G

module load StdEnv/2023
module load minimap2/2.28

module load gcc/12.3
module load samtools/1.20

# IMPORTANT CHANGE!
cd /home/yueyu/scratch/circo_plot_PRA/syri

pra_hap1="/home/yueyu/scratch/Repeatobserver_2026Mar/PRA/hap1/praecox1.fasta"
pra_hap2="/home/yueyu/scratch/Repeatobserver_2026Mar/PRA/hap2/praecox2.fasta"

minimap2 -ax asm5 --eqx $pra_hap1 $pra_hap2 > pra_minimap_20260303.sam

samtools view -b pra_minimap_20260303.sam > pra_minimap_20260303.bam

# ---------------------------------------------------------------------
echo "Job finished with exit code $? at: `date`"
# ---------------------------------------------------------------------

#--------- minimap.align.code.pra.sh (END) --------- 


tmux new-session -s syri_pra_narval3
tmux attach-session -t syri_pra_narval3

cd /home/yueyu/scratch/circo_plot_PRA/syri
module load StdEnv/2020 gcc python/3.9 igraph
source ~/ENV/bin/activate

#CHANGE THIS!!
pra_hap1="/home/yueyu/scratch/Repeatobserver_2026Mar/PRA/hap1/praecox1.fasta"
pra_hap2="/home/yueyu/scratch/Repeatobserver_2026Mar/PRA/hap2/praecox2.fasta"

syri -c pra_minimap_20260303.bam -r $pra_hap1 -q $pra_hap2 -F B -k --prefix pra_20260303

# -k: keep intermediate files
#  20206 March 06 completed



# -- Summary print
cat pra_20260303.summary

# -- Prep for Circos

# Extract inversion > 500Kbp
awk -F'\t' '$11 == "INV"' pra_20260303syri.out > INVERSION_pra_20260303syri.out
# Move to excel -> calc START-END -> Extract Inversion above 500,000bp




# -----------
#
# b: Gene density (1Mbp window)
#
# -----------

# Make genome.chrom.sizes file --> Sort
cd /AsteroidScratch/Annotation_YueYu/LIFTOFF/data

cut -f1,2 debilis1.fasta.fai | sort -k1,1 > /AsteroidScratch/Annotation_YueYu/COUNT_FOR_CIRCO/debilis1.genome.chrom.sizes
cut -f1,2 debilis2.fasta.fai | sort -k1,1 > /AsteroidScratch/Annotation_YueYu/COUNT_FOR_CIRCO/debilis2.genome.chrom.sizes
cut -f1,2 praecox1.fasta.fai | sort -k1,1 > /AsteroidScratch/Annotation_YueYu/COUNT_FOR_CIRCO/praecox1.genome.chrom.sizes
cut -f1,2 praecox2.fasta.fai | sort -k1,1 > /AsteroidScratch/Annotation_YueYu/COUNT_FOR_CIRCO/praecox2.genome.chrom.sizes


cd /AsteroidScratch/Annotation_YueYu/COUNT_FOR_CIRCO
cat debilis1.genome.chrom.sizes
DEB_Chr01	194309941
DEB_Chr02	170093763
DEB_Chr03	186452246
DEB_Chr04	214392016
DEB_Chr05	196277900
DEB_Chr06	198124136
DEB_Chr07	209806723
DEB_Chr08	192682052
DEB_Chr09	183801473
DEB_Chr10	191060906
DEB_Chr11	205358316
DEB_Chr12	243777270
DEB_Chr13	182728780
DEB_Chr14	171009448
DEB_Chr15	169149640
DEB_Chr16	217278156
DEB_Chr17	179148337

cat debilis2.genome.chrom.sizes
DEB_Chr01	191791176
DEB_Chr02	172919073
DEB_Chr03	188981451
DEB_Chr04	213863459
DEB_Chr05	195872972
DEB_Chr06	198360225
DEB_Chr07	208906435
DEB_Chr08	194776242
DEB_Chr09	186015485
DEB_Chr10	189972977
DEB_Chr11	205457193
DEB_Chr12	246124520
DEB_Chr13	181569588
DEB_Chr14	173326191
DEB_Chr15	168104269
DEB_Chr16	213642740
DEB_Chr17	182159180

cat praecox1.genome.chrom.sizes
PRA_chr01	170814815
PRA_chr02	165082315
PRA_chr03	185407962
PRA_chr04	186506442
PRA_chr05	174755163
PRA_chr06	173029008
PRA_chr07	158722642
PRA_chr08	178220724
PRA_chr09	168202118
PRA_chr10	185660960
PRA_chr11	179209817
PRA_chr12	209983383
PRA_chr13	198968881
PRA_chr14	154210255
PRA_chr15	160988758
PRA_chr16	208698999
PRA_chr17	151074730


cat praecox2.genome.chrom.sizes
PRA_chr01	180437650
PRA_chr02	175100439
PRA_chr03	187868931
PRA_chr04	202785240
PRA_chr05	187619502
PRA_chr06	204516750
PRA_chr07	165541277
PRA_chr08	173574252
PRA_chr09	166191611
PRA_chr10	187120831
PRA_chr11	184666794
PRA_chr12	234365546
PRA_chr13	214980597
PRA_chr14	154014317
PRA_chr15	157226489
PRA_chr16	196594418
PRA_chr17	170501492





# ======
# DEB H1 -- 3319, 1Mbp window -- COMPLETED
# ======
GENE="/AsteroidScratch/for_Eric/annotation_result/DEB_H1/braker.gff3"
CHROMSIZE="debilis1.genome.chrom.sizes"
OUT="gene_density_1Mbp_DEB_H1"

cd /AsteroidScratch/Annotation_YueYu/COUNT_FOR_CIRCO
grep -P "\tgene\t" $GENE | awk 'BEGIN{OFS="\t"}{print $1, $4, $5}' > gene.bed
# Change Chr01 --> "DEB_Chr01"
sed -i 's/^Chr/DEB_&/' gene.bed # ADD "DEB_" to all chromosomes
bedtools makewindows -g $CHROMSIZE -w 1000000 > windows.bed  
bedtools intersect -a windows.bed -b gene.bed -c > ${OUT}.bed
# The output (gene_density_1Mbp_DEB_H1.bed) will have columns: Chrom, Start, End, GeneCount. 


# ======
# DEB H2
# ======
GENE="/AsteroidScratch/for_Eric/annotation_result/DEB_H2/braker.gff3"
CHROMSIZE="debilis2.genome.chrom.sizes"
OUT="gene_density_1Mbp_DEB_H2"


# ======
# PRA H1
# ======
GENE="/AsteroidScratch/for_Eric/annotation_result/PRA_H1/praecox_H1_renamed_F.gff3" # Changed to this version Ji checked CHR order: 2026 April 2
CHROMSIZE="praecox1.genome.chrom.sizes"
OUT="gene_density_1Mbp_PRA_H1"

cd /AsteroidScratch/Annotation_YueYu/COUNT_FOR_CIRCO
grep -P "\tgene\t" $GENE | awk 'BEGIN{OFS="\t"}{print $1, $4, $5}' > gene.bed
# Change seq1 --> "PRA_chr01"
paste <(cut -f1 praecox1.genome.chrom.sizes) <(seq 1 17 | awk '{printf "Chr%d\n", $1}') > chr_map.txt
awk 'NR==FNR{map[$2]=$1; next} {if($1 in map) $1=map[$1]; print}' OFS="\t" chr_map.txt gene.bed > gene_chr_name.bed
bedtools makewindows -g $CHROMSIZE -w 1000000 > windows.bed  
bedtools intersect -a windows.bed -b gene_chr_name.bed -c > ${OUT}.bed

# ======
# PRA H2
# ======
GENE="/AsteroidScratch/for_Eric/annotation_result/PRA_H2/praecox_H2_renamed_F.gff3" 
CHROMSIZE="praecox2.genome.chrom.sizes"
OUT="gene_density_1Mbp_PRA_H2"

# Note: Gene density calculation was redone based on Ji update on 2026 April 02




# -----------
#
# c: TE density  -- COMPLETED
#
# -----------

# ======
# DEB H1 
# ======
REPEAT="/AsteroidScratch/for_Eric/annotation_result/DEB_H1/debilis_hap1_2024Dec13.fasta.out.gff"
CHROMSIZE="debilis1.genome.chrom.sizes"
OUT="repeat_density_1Mbp_DEB_H1"

cd /AsteroidScratch/Annotation_YueYu/COUNT_FOR_CIRCO
awk '$3=="dispersed_repeat" {print $1"\t"$4"\t"$5}' $REPEAT > dispersed_repeat.bed
bedtools makewindows -g $CHROMSIZE -w 1000000 > windows.bed  
bedtools intersect -a windows.bed -b dispersed_repeat.bed -c > ${OUT}.bed
# The output (gene_density_1Mbp_DEB_H1.bed) will have columns: Chrom, Start, End, GeneCount. 

# ======
# DEB H2
# ======
REPEAT="/AsteroidScratch/for_Eric/annotation_result/DEB_H2/debilis_hap2_2024Dec13.fasta.out.gff"
CHROMSIZE="debilis2.genome.chrom.sizes"
OUT="repeat_density_1Mbp_DEB_H2"


# ======
# PRA H1
# ======
REPEAT="/AsteroidScratch/for_Eric/annotation_result/PRA_H1/praecox_hap1_2025April07.fasta.out.gff"
CHROMSIZE="praecox1.genome.chrom.sizes"
OUT="repeat_density_1Mbp_PRA_H1"


# ======
# PRA H2
# ======
REPEAT="/AsteroidScratch/for_Eric/annotation_result/PRA_H2/praecox_hap2_2025April07.fasta.out.gff"
CHROMSIZE="praecox2.genome.chrom.sizes"
OUT="repeat_density_1Mbp_PRA_H2"




# -----------
#
# d: Centromere
#
# -----------

# ======
#  DEB 
# ======

# DEB - hap 1
cd /home/yueyu/scratch/Repeatobserver_2026Mar/DEB/hap1/output_chromosomes/debilis_H0-AT/Summary_output/histograms
le debilis_H0-AT_RepeatAbund_centromere_minrange.txt

# DEB - hap 2
cd /home/yueyu/scratch/Repeatobserver_2026Mar/DEB/hap2/output_chromosomes/debilis_H0-AT/Summary_output/histograms
le debilis_H0-AT_RepeatAbund_centromere_minrange.txt



# ======
#  PRA 
# ======

# PRA - hap 1
cd /home/yueyu/scratch/Repeatobserver_2026Mar/PRA/hap1/output_chromosomes/praecox_H0-AT/Summary_output/histograms
le praecox_H0-AT_RepeatAbund_centromere_minrange.txt

# PRA - hap 2
cd /home/yueyu/scratch/Repeatobserver_2026Mar/PRA/hap2/output_chromosomes/praecox_H0-AT/Summary_output/histograms
le praecox_H0-AT_RepeatAbund_centromere_minrange.txt











# -----------
#
# e: TB tools on local machine
#
# -----------

# -- Move all files to local machine

# IMPORTANT NOTE: change all chrom names for H2 genome files 

cd /Users/yueyu/Desktop/circo_2026March11

#================
# Chromosome size: 
#================
# Need to change all CHROM name in H2 genome to "DEB_Chr01_H2" for all following file
awk '{ $1=$1"_H2"; print }' OFS='\t' debilis2.genome.chrom.sizes > debilis2.H2.genome.chrom.sizes
awk '{ $1=$1"_H2"; print }' OFS='\t' praecox2.genome.chrom.sizes > praecox2.H2.genome.chrom.sizes

# Reverse Chr17_H2 -> Chr17_H1
tail -r debilis2.H2.genome.chrom.sizes > reverse_debilis2.H2.genome.chrom.sizes
tail -r praecox2.H2.genome.chrom.sizes > reverse_praecox2.H2.genome.chrom.sizes

cat debilis1.genome.chrom.sizes reverse_debilis2.H2.genome.chrom.sizes > deb_h1_h2_chrom_size.txt
cat praecox1.genome.chrom.sizes reverse_praecox2.H2.genome.chrom.sizes > pra_h1_h2_chrom_size.txt


#================
# gene + repeats: CHROM\tSTART\tEND  -> in "lines"
#================

scp yue_yu@sundance.zoology.ubc.ca:"/AsteroidScratch/Annotation_YueYu/COUNT_FOR_CIRCO/*" .

# H2 name change: gene
cd /Users/yueyu/Desktop/circo_2026March11/gene
awk '{ $1=$1"_H2"; print }' OFS='\t' gene_density_1Mbp_DEB_H2.bed > gene_density_1Mbp_DEB_H2_CHR_changed.bed
awk '{ $1=$1"_H2"; print }' OFS='\t' gene_density_1Mbp_PRA_H2.bed > gene_density_1Mbp_PRA_H2_CHR_changed.bed

# merge
cat gene_density_1Mbp_DEB_H1.bed gene_density_1Mbp_DEB_H2_CHR_changed.bed > deb_all_gene.txt
cat gene_density_1Mbp_PRA_H1.bed gene_density_1Mbp_PRA_H2_CHR_changed.bed > pra_all_gene.txt



# H2 name change: repeat
cd /Users/yueyu/Desktop/circo_2026March11/repeat
awk '{ $1=$1"_H2"; print }' OFS='\t' repeat_density_1Mbp_DEB_H2.bed> repeat_density_1Mbp_DEB_H2_CHR_changed.bed
awk '{ $1=$1"_H2"; print }' OFS='\t' repeat_density_1Mbp_PRA_H2.bed > repeat_density_1Mbp_PRA_H2_CHR_changed.bed

# merge
cat repeat_density_1Mbp_DEB_H1.bed repeat_density_1Mbp_DEB_H2_CHR_changed.bed > deb_all_repeat.txt
cat repeat_density_1Mbp_PRA_H1.bed repeat_density_1Mbp_PRA_H2_CHR_changed.bed > pra_all_repeat.txt






#================
# Centeromere: CHROM\tSTART\tEND\t"1"  -> in "Tiles"
#================

scp yueyu@narval.computecanada.ca:"/home/yueyu/scratch/Repeatobserver_2026Mar/PRA/hap1/output_chromosomes/praecox_H0-AT/Summary_output/histograms/praecox_H0-AT_RepeatAbund_centromere_minrange.txt" .
mv praecox_H0-AT_RepeatAbund_centromere_minrange.txt pra_h1_centromere.txt

scp yueyu@narval.computecanada.ca:"/home/yueyu/scratch/Repeatobserver_2026Mar/PRA/hap2/output_chromosomes/praecox_H0-AT/Summary_output/histograms/praecox_H0-AT_RepeatAbund_centromere_minrange.txt" .
mv praecox_H0-AT_RepeatAbund_centromere_minrange.txt pra_h2_centromere.txt

scp yueyu@narval.computecanada.ca:"/home/yueyu/scratch/Repeatobserver_2026Mar/DEB/hap1/output_chromosomes/debilis_H0-AT/Summary_output/histograms/debilis_H0-AT_RepeatAbund_centromere_minrange.txt" .
mv debilis_H0-AT_RepeatAbund_centromere_minrange.txt deb_h1_centromere.txt

scp yueyu@narval.computecanada.ca:"/home/yueyu/scratch/Repeatobserver_2026Mar/DEB/hap2/output_chromosomes/debilis_H0-AT/Summary_output/histograms/debilis_H0-AT_RepeatAbund_centromere_minrange.txt" .
mv debilis_H0-AT_RepeatAbund_centromere_minrange.txt deb_h2_centromere.txt


# Summarize range in R
R

a <- read.delim("deb_h1_centromere.txt", header = F)
a <- a[,c(2,3,4)]
a$V2 <- paste0("DEB_Chr", sprintf("%02d", as.integer(a$V2)))
write.table(a, file="deb_h1_centromere_ready_circo.txt", sep = "\t", col.names = F, row.names = F, quote=F)


a <- read.delim("deb_h2_centromere.txt", header = F)
a <- a[,c(2,3,4)]
a$V2 <- paste0("DEB_Chr", sprintf("%02d", as.integer(a$V2)))
write.table(a, file="deb_h2_centromere_ready_circo.txt", sep = "\t", col.names = F, row.names = F, quote=F)

a <- read.delim("pra_h1_centromere.txt", header = F)
a <- a[,c(2,3,4)]
a$V2 <- paste0("PRA_chr", sprintf("%02d", as.integer(a$V2)))
write.table(a, file="pra_h1_centromere_ready_circo.txt", sep = "\t", col.names = F, row.names = F, quote=F)

a <- read.delim("pra_h2_centromere.txt", header = F)
a <- a[,c(2,3,4)]
a$V2 <- paste0("PRA_chr", sprintf("%02d", as.integer(a$V2)))
write.table(a, file="pra_h2_centromere_ready_circo.txt", sep = "\t", col.names = F, row.names = F, quote=F)

# H2 name change
awk '{ $1=$1"_H2"; print }' OFS='\t' deb_h2_centromere_ready_circo.txt > deb_h2_centromere_ready_circo_CHR_changed.txt
awk '{ $1=$1"_H2"; print }' OFS='\t' pra_h2_centromere_ready_circo.txt > pra_h2_centromere_ready_circo_CHR_changed.txt

# merge
cat deb_h1_centromere_ready_circo.txt deb_h2_centromere_ready_circo_CHR_changed.txt > deb_all_centromere.tsv
cat pra_h1_centromere_ready_circo.txt pra_h2_centromere_ready_circo_CHR_changed.txt > pra_all_centromere.tsv


# Remove 8.3e+07 into numeric for TBtools to readin
awk '{printf "%s\t%d\t%d\n",$1,$2,$3}' deb_all_centromere.tsv | awk 'BEGIN{OFS="\t"}{print $0,1}' > deb_all_centromere_NUM.tsv
awk '{printf "%s\t%d\t%d\n",$1,$2,$3}' pra_all_centromere.tsv | awk 'BEGIN{OFS="\t"}{print $0,1}' > pra_all_centromere_NUM.tsv



#================
# syri  CHROM\tSTART\tEND\tCHROM_H2\tSTART\tEND  -> add in the beginning
#================

scp yueyu@narval.computecanada.ca:"/home/yueyu/scratch/circo_plot/syri/move_to_local/*.txt" .

cat deb*synOut.txt | grep "#" | awk 'BEGIN{OFS="\t"} {print $2,$3,$4,$6"_H2",$7,$8}' > deb_synOut_Circos.tsv
cat deb*invOut.txt | grep "#" | awk 'BEGIN{OFS="\t"} {print $2,$3,$4,$6"_H2",$7,$8}' > deb_invOut_Circos.tsv
cat deb*TLOut.txt | grep "#" | awk 'BEGIN{OFS="\t"} {print $2,$3,$4,$6"_H2",$7,$8}' > deb_TLOut_Circos.tsv
cat deb*invTLOut.txt | grep "#" | awk 'BEGIN{OFS="\t"} {print $2,$3,$4,$6"_H2",$7,$8}' > deb_invTLOut_Circos.tsv
cat deb*dupOut.txt | grep "#" | awk 'BEGIN{OFS="\t"} {print $2,$3,$4,$6"_H2",$7,$8}' > deb_dupOut_Circos.tsv
cat deb*invDupOut.txt | grep "#" | awk 'BEGIN{OFS="\t"} {print $2,$3,$4,$6"_H2",$7,$8}' > deb_invDupOut_Circos.tsv

# merge
cat deb*Circos.tsv > deb_all_Circos.tsv

# OPTIONAL: add color to each type of SV
# OPTIONAL: subset syntenic region to >5Mbp
# OPTIONAL: only plot SYN + add color to each chromosome + highlight candidtae INV





scp yueyu@narval.computecanada.ca:"/home/yueyu/scratch/circo_plot_PRA/syri/move_to_local/*.txt" .

cat pra*synOut.txt | grep "#" | awk 'BEGIN{OFS="\t"} {print $2,$3,$4,$6"_H2",$7,$8}' > pra_synOut_Circos.tsv
cat pra*invOut.txt | grep "#" | awk 'BEGIN{OFS="\t"} {print $2,$3,$4,$6"_H2",$7,$8}' > pra_invOut_Circos.tsv
cat pra*TLOut.txt | grep "#" | awk 'BEGIN{OFS="\t"} {print $2,$3,$4,$6"_H2",$7,$8}' > pra_TLOut_Circos.tsv
cat pra*invTLOut.txt | grep "#" | awk 'BEGIN{OFS="\t"} {print $2,$3,$4,$6"_H2",$7,$8}' > pra_invTLOut_Circos.tsv
cat pra*dupOut.txt | grep "#" | awk 'BEGIN{OFS="\t"} {print $2,$3,$4,$6"_H2",$7,$8}' > pra_dupOut_Circos.tsv
cat pra*invDupOut.txt | grep "#" | awk 'BEGIN{OFS="\t"} {print $2,$3,$4,$6"_H2",$7,$8}' > pra_invDupOut_Circos.tsv
# merge
cat pra*Circos.tsv > pra_all_Circos.tsv


# COMPLETED

# END