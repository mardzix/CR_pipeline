#!/usr/bin/python3

#SBATCH-n 16
#SBATCH -t 0-10:00
#SBATCH -A snic2017-7-345
#SBATCH -o slurm_%j.out \
#SBATCH -e slurm_%j.err \

import sys
import os
import fnmatch
import argparse
from datetime import datetime

parser = argparse.ArgumentParser("Cut and run pipeline")

parser.add_argument("-i","--input",dest="input",type=str,help="Path to input file folder. The script will look into subfolders and retieve all fastq files")
parser.add_argument("-p","--pattern",dest="pattern",type=str,help="Regexp pattern for which the script should look for")
parser.add_argument("-o","--output",dest="out",type=str,help="Output folder")
parser.add_argument("-n","--name",dest="name",type=str,help="Name of the experiment",default=False)
args = parser.parse_args()



FASTQ_PATH      = args.input    #"/proj/uppstore2017150/private/marek/G.CasteloBranco_19_05-P13556/raw/P13556_1001/02-FASTQ/190625_A00621_0091_BHCKNLDRXX/"
PATTERN         = args.pattern   #"P13556*.fastq.gz"
OUTDIR          = args.out       #"/proj/uppstore2017150/private/marek/G.CasteloBranco_19_05-P13556/CR_pipeline_MB/"
EXPERIMENT_NAME = args.name

CR_PATH          = "~/bin/CR_pipeline"

TRIMMOMATIC      = "trimmomatic/0.36"
SAMTOOLS         = "samtools/1.9"
BOWTIE2          = "bowtie2/2.2.9"
PICARD           = "picard/2.10.3"
PICARD_JAR_PATH  = "/sw/apps/bioinfo/picard/2.10.3/milou/picard.jar" 
DEEPTOOLS        = "deepTools/3.1.0"
MACS             = "MACS/2.1.2"
BEDOPS           = "BEDOPS/2.4.28"
MEME             = "MEMESuite/5.0.1"
BEDTOOLS         = "BEDTools/2.27.1"
OPENMPI          = "openmpi/3.1.1"

NTHREADS         = 16
BOWTIE2_INDEX    =  "/proj/uppstore2017150/private/marek/index/bowtie2/mm10_iGenomes/Mus_musculus/Ensembl/GRCm38/Sequence/Bowtie2Index/mm10"
GENOME_FA        = "/proj/uppstore2017150/private/marek/index/bowtie2/mm10_iGenomes/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta/genome.fa"
#BOWTIE2_SPIKEIN  =  "/proj/uppstore2017150/private/marek/index/bowtie2/cerevisiae/s_cerevisiae"
BOWTIE2_SPIKEIN  =  "/proj/uppstore2017150/private/marek/index/bowtie2/ecoli/e_coli"
ADAPTER_PATH     = "/home/marek/bin/cutruntools/adapters"
SIZE             = 120
GENOME           = "mm10"


exec(open('/home/marek/bin/module/init/python.py').read())
module('load','bioinfo-tools')
module('load',TRIMMOMATIC)
module('load',SAMTOOLS)
module('load',BOWTIE2)
module('load',PICARD)
module('load',DEEPTOOLS)
module('load',MACS)
module('load',BEDOPS)
module('load',BEDTOOLS)
module('load',MEME)

# MEME gave error when loading with pythonic module()
os.system("module load " + MEME)


def search_for_fastq(path,pattern):
  result = []
  for r,d,f in os.walk(path):
    for name in f:
      if fnmatch.fnmatch(name, pattern):
        result.append(os.path.join(r, name))
  return (result)
  

def main(OUTDIR=OUTDIR,PATTERN=PATTERN,FASTQ_PATH=FASTQ_PATH,EXPERIMENT_NAME=EXPERIMENT_NAME):
  startTime = datetime.now()
    
  FASTQ_FILES=search_for_fastq(FASTQ_PATH,PATTERN)
  sys.stderr.write("*** {0} *** MB_pipeline: Files to be analyzed \n{1}\n{2}\n".format(datetime.now(),FASTQ_FILES[0],FASTQ_FILES[1]))
    
  if len(FASTQ_FILES) != 2:
    sys.exit("*** MB_pipeline: More than 2 fastq files present in root folder")
  
  if not EXPERIMENT_NAME:
    EXPERIMENT_NAME = "_".join(os.path.basename(FASTQ_FILES[1]).split("_")[0:2])
  
  if not os.path.exists(OUTDIR):
    os.mkdir("{}".format(OUTDIR))
  
  OUTDIR = OUTDIR + EXPERIMENT_NAME
  
  if not os.path.exists(OUTDIR):
    os.mkdir("{}/".format(OUTDIR))
  
  
  
  ##############
  # TRRIMMING  #
  ##############
  
  if not os.path.exists(OUTDIR + "/trimmed/"):
    os.mkdir("{}/trimmed".format(OUTDIR))
    
  TRIM_OUTDIR     = OUTDIR + "/trimmed/"
  TRIM_FILE1       = TRIM_OUTDIR + os.path.basename(FASTQ_FILES[0]).replace(".fastq.gz","_trimmed.fastq.gz")
  TRIM_UNPAIRED1   = TRIM_OUTDIR + os.path.basename(FASTQ_FILES[0]).replace(".fastq.gz","_unpaired.fastq.gz")
  TRIM_FILE2       = TRIM_OUTDIR + os.path.basename(FASTQ_FILES[1]).replace(".fastq.gz","_trimmed.fastq.gz")
  TRIM_UNPAIRED2   = TRIM_OUTDIR + os.path.basename(FASTQ_FILES[1]).replace(".fastq.gz","_unpaired.fastq.gz")
  
  TRIM_CMD = "trimmomatic PE -threads {7} \
                             -phred33 \
                             {0} \
                             {1} \
                             {2} \
                             {3} \
                             {4} \
                             {5} \
                             ILLUMINACLIP:{6}/Truseq3.PE.fa:2:15:4:4:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:25 2>&1".format(FASTQ_FILES[0],
                                                                                                                                       FASTQ_FILES[1],
                                                                                                                                      TRIM_FILE1,
                                                                                                                                      TRIM_UNPAIRED1,
                                                                                                                                      TRIM_FILE2,
                                                                                                                                      TRIM_UNPAIRED2,
                                                                                                                                      ADAPTER_PATH,
                                                                                                                                      NTHREADS)
  sys.stderr.write("*** {0} ***MB_pipeline: Performing trimming using {1}\n".format(datetime.now(),TRIMMOMATIC))
  os.system(TRIM_CMD)
  
  #############
  # Alignment #
  #############
  
  if not os.path.exists("{}/bowtie2_align/".format(OUTDIR)):
    os.mkdir("{}/bowtie2_align/".format(OUTDIR))
  
  BOWTIE_CMD = "bowtie2 -p {5} \
                         --dovetail \
                        --phred33 \
                        -x {0} \
                        -1 {1} \
                        -2 {2} \
                        | samtools view -bS - > {3}/bowtie2_align/{4}.bam".format(BOWTIE2_INDEX,
                                                                                  TRIM_FILE1,
                                                                                  TRIM_FILE2,
                                                                                  OUTDIR,
                                                                                  EXPERIMENT_NAME,
                                                                                  NTHREADS,
                                                                                  TRIM_OUTDIR)
  
  
  BOWTIE_SPIKEIN     = "bowtie2 -p {5} \
                        --dovetail \
                        --phred33 \
                        -x {0} \
                        -1 {1} \
                        -2 {2} \
                        | samtools view  -f 2 -bS - > {3}/bowtie2_align/{4}_spikein.bam".format(BOWTIE2_SPIKEIN,
                                                                                                TRIM_FILE1,
                                                                                                TRIM_FILE2,
                                                                                                OUTDIR,
                                                                                                EXPERIMENT_NAME,
                                                                                                NTHREADS,
                                                                                                TRIM_OUTDIR)
  SPIKEIN_COUNT        = "samtools view -q 1 -f 2 {0}/bowtie2_align/{1}_spikein.bam | wc -l > {0}/bowtie2_align/{1}_spikein_norm.txt".format(OUTDIR,EXPERIMENT_NAME)
  
  sys.stderr.write("*** {0} *** MB_pipeline: Mapping trimmed reads with {1}\n".format(datetime.now(),BOWTIE2))
  os.system(BOWTIE_CMD)
  sys.stderr.write("*** {0} *** MB_pipeline: Mapping trimmed reads to the spike-in genome with {1}\n".format(datetime.now(),BOWTIE2))
  os.system(BOWTIE_SPIKEIN)
  sys.stderr.write("*** {0} *** MB_pipeline: Calculating spike-in normalisation factor {1}\n".format(datetime.now(),BOWTIE2))
  os.system(SPIKEIN_COUNT)
  
  #####################################
  # PICARD MANIPULATIONS OF BAM FILES #
  #####################################
  
  if not os.path.exists("{}/bowtie2_align/sorted/".format(OUTDIR)):
    os.mkdir("{}/bowtie2_align/sorted".format(OUTDIR))
  
  PICARD_SORT = "java -jar {0} SortSam INPUT={1}/bowtie2_align/{2}.bam \
                                       OUTPUT={1}/bowtie2_align/sorted/{2}_sorted.bam \
                                       SORT_ORDER=coordinate".format(PICARD_JAR_PATH,
                                                                     OUTDIR,
                                                                     EXPERIMENT_NAME)
  
  # Picard mark is not run anymore (os.system(PICARD_MARK) to run it)
  PICARD_MARK = "java -jar {0} MarkDuplicates INPUT={1}/bowtie2_align/sorted/{2}_sorted.bam \
                                              OUTPUT={1}/bowtie2_align/sorted/{2}_marked-dup.bam \
                                              METRICS_FILE={1}/bowtie2_align/sorted/{2}_marked-dup_metrics.txt".format(PICARD_JAR_PATH,
                                                                                                               OUTDIR,
                                                                                                               EXPERIMENT_NAME)
  
  PICARD_REMOVE_DUP = "java -jar {0} MarkDuplicates INPUT={1}/bowtie2_align/sorted/{2}_sorted.bam\
                                                    OUTPUT={1}/bowtie2_align/sorted/{2}_dedup.bam \
                                                    METRICS_FILE={1}/bowtie2_align/sorted/{2}_dedup_metrics.txt\
                                                    REMOVE_DUPLICATES=true".format(PICARD_JAR_PATH,
                                                                                   OUTDIR,
                                                                                   EXPERIMENT_NAME)
  
  sys.stderr.write("*** {0} *** MB_pipeline: Sorting bam files {1} with {2}\n".format(datetime.now(),EXPERIMENT_NAME,PICARD))
  os.system(PICARD_SORT)


  sys.stderr.write("*** {0} *** MB_pipeline: Removing duplicates {1} with {2}\n".format(datetime.now(),EXPERIMENT_NAME,PICARD))
  os.system(PICARD_REMOVE_DUP)
  
  #################################### Filter <120bp reads (sub-nucleosome fraction)
  
  SAMTOOLS_120_FILTER_CMD          = "samtools view -h {0}/bowtie2_align/sorted/{1}_sorted.bam | awk -f {3}/bin/filter_below.awk | samtools view -bS - > {0}/bowtie2_align/sorted/{1}_sorted_120bp.bam".format(OUTDIR, EXPERIMENT_NAME, SIZE,CR_PATH)
  SAMTOOLS_120_FILTER_CMD_MARKED   = "samtools view -h {0}/bowtie2_align/sorted/{1}_marked-dup.bam | awk -f {3}/bin/filter_below.awk | samtools view -bS - > {0}/bowtie2_align/sorted/{1}_marked-dup_120bp.bam".format(OUTDIR, EXPERIMENT_NAME, SIZE,CR_PATH)
  SAMTOOLS_120_FILTER_CMD_DEDUP    = "samtools view -h {0}/bowtie2_align/sorted/{1}_dedup.bam  | awk -f {3}/bin/filter_below.awk | samtools view -bS - > {0}/bowtie2_align/sorted/{1}_dedup_120bp.bam".format(OUTDIR, EXPERIMENT_NAME, SIZE, CR_PATH)
  
  sys.stderr.write("*** {0} *** MB_pipeline: Extracting reads < 120 bp from {1} with {2} \n".format(datetime.now(),EXPERIMENT_NAME, SAMTOOLS))

  os.system(SAMTOOLS_120_FILTER_CMD)
  os.system(SAMTOOLS_120_FILTER_CMD_DEDUP)
  
  
  ###### Index all bam files
  
  INDEX1  =  "samtools index {0}/bowtie2_align/sorted/{1}_sorted.bam".format(OUTDIR,EXPERIMENT_NAME)
  INDEX2  =  "samtools index {0}/bowtie2_align/sorted/{1}_marked-dup.bam".format(OUTDIR,EXPERIMENT_NAME)
  INDEX3  =  "samtools index {0}/bowtie2_align/sorted/{1}_dedup.bam".format(OUTDIR,EXPERIMENT_NAME)
  INDEX4  =  "samtools index {0}/bowtie2_align/sorted/{1}_sorted_120bp.bam".format(OUTDIR,EXPERIMENT_NAME)
  INDEX5  =  "samtools index {0}/bowtie2_align/sorted/{1}_marked-dup_120bp.bam".format(OUTDIR,EXPERIMENT_NAME)
  INDEX6  =  "samtools index {0}/bowtie2_align/sorted/{1}_dedup_120bp.bam".format(OUTDIR,EXPERIMENT_NAME)

  sys.stderr.write("*** {0} *** MB_pipeline: Indexing bam files from {1} with {2} \n".format(datetime.now(),EXPERIMENT_NAME, SAMTOOLS))  
  os.system(INDEX1)
  os.system(INDEX3)
  os.system(INDEX4)
  os.system(INDEX6)  
  
  ################################
  # Create bw for genome browser #
  ################################
  
  
  if not os.path.exists("{}/bigwig/".format(OUTDIR)):
    os.mkdir("{}/bigwig/".format(OUTDIR))
  
  SPIKEIN_COUNT = int(open("{0}/bowtie2_align/{1}_spikein_norm.txt".format(OUTDIR,EXPERIMENT_NAME)).readline().strip())
  if SPIKEIN_COUNT == 0:
    SPIKEIN_COUNT = 1
  
  SCALE_FACTOR = 1000.0/SPIKEIN_COUNT
#  SCALE_FACTOR = 100
  
  BIGWIG_CMD1 = "bamCoverage --scaleFactor {2} --numberOfProcessors {3} --ignoreDuplicates --minMappingQuality 5 -b {0}/bowtie2_align/sorted/{1}_sorted.bam -o {0}/bigwig/{1}.bw --binSize 10 --centerReads --smoothLength 50" .format(OUTDIR,EXPERIMENT_NAME,SCALE_FACTOR,NTHREADS)
  BIGWIG_CMD2 = "bamCoverage --scaleFactor {2} --numberOfProcessors {3} --ignoreDuplicates --minMappingQuality 5 -b {0}/bowtie2_align/sorted/{1}_marked-dup.bam -o {0}/bigwig/{1}_marked-dup.bw --binSize 10 --centerReads --smoothLength 50" .format(OUTDIR,EXPERIMENT_NAME,SCALE_FACTOR,NTHREADS)
  BIGWIG_CMD3 = "bamCoverage --scaleFactor {2} --numberOfProcessors {3} --ignoreDuplicates --minMappingQuality 5 -b {0}/bowtie2_align/sorted/{1}_dedup.bam -o {0}/bigwig/{1}_dedup.bw --binSize 10 --centerReads --smoothLength 50" .format(OUTDIR,EXPERIMENT_NAME,SCALE_FACTOR,NTHREADS)
  BIGWIG_CMD4 = "bamCoverage --scaleFactor {2} --numberOfProcessors {3} --ignoreDuplicates --minMappingQuality 5 -b {0}/bowtie2_align/sorted/{1}_sorted_120bp.bam -o {0}/bigwig/{1}_sorted_120bp.bw --binSize 10 --centerReads --smoothLength 50" .format(OUTDIR,EXPERIMENT_NAME,SCALE_FACTOR,NTHREADS)
  BIGWIG_CMD5 = "bamCoverage --scaleFactor {2} --numberOfProcessors {3} --ignoreDuplicates --minMappingQuality 5 -b {0}/bowtie2_align/sorted/{1}_marked-dup_120bp.bam -o {0}/bigwig/{1}_marked-dup_120bp.bw --binSize 10 --centerReads --smoothLength 50" .format(OUTDIR,EXPERIMENT_NAME,SCALE_FACTOR,NTHREADS)
  BIGWIG_CMD6 = "bamCoverage --scaleFactor {2} --numberOfProcessors {3} --ignoreDuplicates --minMappingQuality 5 -b {0}/bowtie2_align/sorted/{1}_dedup_120bp.bam -o {0}/bigwig/{1}_dedup_120bp.bw --binSize 10 --centerReads --smoothLength 50" .format(OUTDIR,EXPERIMENT_NAME,SCALE_FACTOR,NTHREADS)
  
  BIGWIG_CMD_RPKM = "bamCoverage --normalizeUsing RPKM --numberOfProcessors {3} --ignoreDuplicates --minMappingQuality 5 -b {0}/bowtie2_align/sorted/{1}_dedup.bam -o {0}/bigwig/{1}_dedup_RPKM.bw --binSize 10 --centerReads --smoothLength 50" .format(OUTDIR,EXPERIMENT_NAME,SCALE_FACTOR,NTHREADS)
  
  sys.stderr.write("*** {0} *** MB_pipeline: Generating bigwig files {1} using bamCoverage from {2} \n".format(datetime.now(),EXPERIMENT_NAME, DEEPTOOLS))    
  
  os.system(BIGWIG_CMD1)  
  os.system(BIGWIG_CMD3)
  os.system(BIGWIG_CMD4)
  os.system(BIGWIG_CMD6)
  os.system(BIGWIG_CMD_RPKM)
  
  if not os.path.exists("{}/macs/".format(OUTDIR)):
    os.mkdir("{}/macs/".format(OUTDIR))
    
  if not os.path.exists("{}/macs.nodup/".format(OUTDIR)):
    os.mkdir("{}/macs.nodup/".format(OUTDIR))
  
  
  ##########################
  # Peak calling with MACS #
  ##########################
  
  MACS_CMD = "macs2 callpeak  -t {0}/bowtie2_align/sorted/{1}_sorted.bam \
                              -g {2} \
                              -f BAMPE \
                              -n {1} \
                              --outdir {0}/macs \
                              -q 0.01 \
                              -B \
                              --SPMR \
                              --keep-dup all 2>&1".format(OUTDIR,EXPERIMENT_NAME,GENOME.rstrip("0123456789"))
  
  MACS_CMD_BROAD = "macs2 callpeak  -t {0}/bowtie2_align/sorted/{1}_sorted.bam \
                                    -g {2} \
                                    -f BAMPE \
                                    -n {1} \
                                    --outdir {0}/macs_broad \
                                    -q 0.01 \
                                    -B \
                                    --SPMR \
                                    --broad \
                                    --keep-dup all 2>&1".format(OUTDIR,EXPERIMENT_NAME,GENOME.rstrip("0123456789"))
  
  MACS_CMD_DEDUP = "macs2 callpeak  -t {0}/bowtie2_align/sorted/{1}_dedup.bam \
                                    -g {2} \
                                    -f BAMPE \
                                    -n {1} \
                                    --outdir {0}/macs.nodup \
                                    -q 0.01 \
                                    -B \
                                    --SPMR \
                                    --keep-dup all 2>&1".format(OUTDIR,EXPERIMENT_NAME,GENOME.rstrip("0123456789"))
  
  MACS_CMD_DEDUP_BROAD = "macs2 callpeak  -t {0}/bowtie2_align/sorted/{1}_dedup.bam \
                                          -g {2} \
                                          -f BAMPE \
                                          -n {1} \
                                          --outdir {0}/macs.nodup_broad \
                                          -q 0.01 \
                                          -B \
                                          --SPMR \
                                          --broad \
                                          --keep-dup all 2>&1".format(OUTDIR,EXPERIMENT_NAME,GENOME.rstrip("0123456789"))
                                          
  MACS_CMD_DEDUP_120 = "macs2 callpeak  -t {0}/bowtie2_align/sorted/{1}_dedup_120bp.bam \
                                        -g {2} \
                                        -f BAMPE \
                                        -n {1} \
                                        --outdir {0}/macs.nodup.120 \
                                        -q 0.01 \
                                        -B \
                                        --SPMR \
                                        --keep-dup all 2>&1".format(OUTDIR,EXPERIMENT_NAME,GENOME.rstrip("0123456789"))
  
  sys.stderr.write("*** {0} *** MB_pipeline: Peak calling for {1} using {2} \n".format(datetime.now(),EXPERIMENT_NAME,MACS))    
  os.system(MACS_CMD)
  os.system(MACS_CMD_BROAD)
  os.system(MACS_CMD_DEDUP)
  os.system(MACS_CMD_DEDUP_BROAD)
  os.system(MACS_CMD_DEDUP_120)

  # Plot heatmap of reads around the called peaks
  sys.stderr.write("*** {0} *** MB_pipeline: Plotting heatmap around peaks for {1} \n".format(datetime.now(),EXPERIMENT_NAME,MACS))    
  
  HEATMAP_PLOT = "computeMatrix reference-point -S {0}/bigwig/{1}_dedup.bw -R {0}/macs.nodup/{1}_summits.bed -a 1000 -b 1000 -out {0}/macs.nodup/{1}_matrix -bs 25 -p 1 --missingDataAsZero; \
  plotHeatmap -m {0}/macs.nodup/{1}_matrix --colorList blue,yellow,red --heatmapHeight 25 --heatmapWidth 5 --samplesLabel {1} -out {0}/macs.nodup/{1}_matrix.png --whatToShow 'heatmap and colorbar'  --sortUsing max  --refPointLabel '' ".format(OUTDIR,EXPERIMENT_NAME)
  os.system(HEATMAP_PLOT)

  #############################
  # MEME to search for motifs # 
  #############################
  
  MACS_INDIR = "macs.nodup"
  MEME_OUTDIR = "MEME.nodup"

  if not os.path.exists("{0}/{1}/".format(OUTDIR,MEME_OUTDIR)):
    os.mkdir("{0}/{1}/".format(OUTDIR,MEME_OUTDIR))
  


  sys.stderr.write("*** {0} *** MB_pipeline: Preparing fasta files for motif search {1} using {2} \n".format(datetime.now(),EXPERIMENT_NAME,MEME))

  SUMMITS_PADDED_CMD = "bedops --range 150 -u {0}/{2}/{1}_summits.bed  > {0}/{3}/{1}_summits_padded.bed".format(OUTDIR,EXPERIMENT_NAME,MACS_INDIR,MEME_OUTDIR)
  os.system(SUMMITS_PADDED_CMD)

  SUMMITS_PADDED_GETFASTA_CMD = "bedtools getfasta -fi {0} -bed {1}/{3}/{4}_summits_padded.bed -fo {1}/{3}/{4}_summits_padded.fa".format(GENOME_FA,OUTDIR,MACS_INDIR,MEME_OUTDIR,EXPERIMENT_NAME)
  os.system(SUMMITS_PADDED_GETFASTA_CMD)

  # Filte N and lower case atcg contatining fasta entries
  FILTER_FA_CMD = "python {4}/bin/filter.py {0}/{3}/{1}_summits_padded.fa 300 > {0}/{3}/{1}_summits_repeats.bed".format(OUTDIR,EXPERIMENT_NAME,MACS_INDIR,MEME_OUTDIR,CR_PATH)
  os.system(FILTER_FA_CMD)

  SUMMITS_FILTER_CMD = "bedops -n 1 {0}/{2}/{1}_summits.bed {0}/{3}/{1}_summits_repeats.bed | sort-bed - > {0}/{3}/{1}_summits_filtered.bed".format(OUTDIR,EXPERIMENT_NAME,MACS_INDIR,MEME_OUTDIR)
  os.system(SUMMITS_FILTER_CMD)

  # Filter blacklisted regions
  FINAL_SUMMIT_FILTER_CMD = "cat {0}/{3}/{1}_summits_filtered.bed | grep -v -e 'chrM' | sort-bed - | bedops -n 1 - {5}/bin/blacklist/{4}.blacklist.bed > {0}/{3}/{1}_summits_filtered_final.bed".format(OUTDIR,EXPERIMENT_NAME,MACS_INDIR,MEME_OUTDIR,GENOME,CR_PATH)
  os.system(FINAL_SUMMIT_FILTER_CMD)

  # Get 1000 random peaks from top 5000 peaks
  RANDOMIZE = "sort -k5gr {0}{3}/{1}_summits_filtered.bed | head -5000 | shuf | head -1000 | sort-bed -  > {0}/{3}/{1}_1000_random_summits.bed".format(OUTDIR,EXPERIMENT_NAME,MACS_INDIR,MEME_OUTDIR) 
  os.system(RANDOMIZE)
  
  # Get fasta for random 1000 / top 5000
  PICK_1000 = "sort -k5gr {0}/{3}/{1}_summits_filtered.bed | head -1000 | sort-bed -  > {0}/{3}/{1}_1000_random_summits.bed; \
               bedops --range 150 -u {0}/{3}/{1}_1000_random_summits.bed > {0}/{3}/{1}_1000_random_summits_padded.bed; \
               bedtools getfasta -fi {4} -bed {0}/{3}/{1}_1000_random_summits_padded.bed -fo {0}/{3}/{1}_1000_random_summits_padded.fa".format(OUTDIR,EXPERIMENT_NAME,MACS_INDIR,MEME_OUTDIR,GENOME_FA)
  os.system(PICK_1000)
  
  PICK_5000 = "sort -k5gr {0}/{3}/{1}_summits_filtered.bed | head -5000 | sort-bed -  > {0}/{3}/{1}_5000_random_summits.bed; \
               bedops --range 150 -u {0}/{3}/{1}_5000_random_summits.bed > {0}/{3}/{1}_5000_random_summits_padded.bed; \
               bedtools getfasta -fi {4} -bed {0}/{3}/{1}_5000_random_summits_padded.bed -fo {0}/{3}/{1}_5000_random_summits_padded.fa".format(OUTDIR,EXPERIMENT_NAME,MACS_INDIR,MEME_OUTDIR,GENOME_FA)
  os.system(PICK_5000)
  
  PICK_ALL = "bedops --range 150 -u {0}/{3}/{1}_summits_filtered.bed | sort-bed - > {0}/{3}/{1}_all_summits_padded.bed; \
            bedtools getfasta -fi {4} -bed {0}/{3}/{1}_all_summits_padded.bed -fo {0}/{3}/{1}_all_summits_padded.fa".format(OUTDIR,EXPERIMENT_NAME,MACS_INDIR,MEME_OUTDIR,GENOME_FA)
  os.system(PICK_ALL)

  
  sys.stderr.write("*** {0} *** MB_pipeline: Running motif search on {1} using {2} \n".format(datetime.now(),EXPERIMENT_NAME,MEME))    

  MEME_RUN_1000 = "meme-chip -oc {0}/{3}/MEME_1000/ -dreme-m 10 -meme-nmotifs 10 {0}/{3}/{1}_1000_random_summits_padded.fa".format(OUTDIR,EXPERIMENT_NAME,MACS_INDIR,MEME_OUTDIR)
  os.system(MEME_RUN_1000 )

  MEME_RUN_5000 = "meme-chip -oc {0}/{3}/MEME_5000/ -dreme-m 10 -meme-nmotifs 10 {0}/{3}/{1}_5000_random_summits_padded.fa".format(OUTDIR,EXPERIMENT_NAME,MACS_INDIR,MEME_OUTDIR)
  os.system(MEME_RUN_5000)

  MEME_RUN_ALL = "meme-chip -oc {0}/{3}/MEME_ALL/ -dreme-m 10 -meme-nmotifs 10 {0}/{3}/{1}_all_summits_padded.fa".format(OUTDIR,EXPERIMENT_NAME,MACS_INDIR,MEME_OUTDIR)
  os.system(MEME_RUN_ALL)




  ########################################################################### 
  # MEME to search for motifs in < 120bp fragments # TODO -> SIMPLIFY THIS? #
  ###########################################################################
  
  # 120 bp

  MACS_INDIR = "macs.nodup.120"
  MEME_OUTDIR = "MEME.nodup.120"

  if not os.path.exists("{0}/{1}/".format(OUTDIR,MEME_OUTDIR)):
    os.mkdir("{0}/{1}/".format(OUTDIR,MEME_OUTDIR))
  


  sys.stderr.write("*** {0} *** MB_pipeline: Preparing fasta files for motif search {1} using {2} \n".format(datetime.now(),EXPERIMENT_NAME,MEME))

  SUMMITS_PADDED_CMD = "bedops --range 150 -u {0}/{2}/{1}_summits.bed  > {0}/{3}/{1}_summits_padded.bed".format(OUTDIR,EXPERIMENT_NAME,MACS_INDIR,MEME_OUTDIR)
  os.system(SUMMITS_PADDED_CMD)

  SUMMITS_PADDED_GETFASTA_CMD = "bedtools getfasta -fi {0} -bed {1}/{3}/{4}_summits_padded.bed -fo {1}/{3}/{4}_summits_padded.fa".format(GENOME_FA,OUTDIR,MACS_INDIR,MEME_OUTDIR,EXPERIMENT_NAME)
  os.system(SUMMITS_PADDED_GETFASTA_CMD)

  # Filte N and lower case atcg contatining fasta entries
  FILTER_FA_CMD = "python {4}/bin/filter.py {0}/{3}/{1}_summits_padded.fa 300 > {0}/{3}/{1}_summits_repeats.bed".format(OUTDIR,EXPERIMENT_NAME,MACS_INDIR,MEME_OUTDIR,CR_PATH)
  os.system(FILTER_FA_CMD)

  SUMMITS_FILTER_CMD = "bedops -n 1 {0}/{2}/{1}_summits.bed {0}/{3}/{1}_summits_repeats.bed | sort-bed - > {0}/{3}/{1}_summits_filtered.bed".format(OUTDIR,EXPERIMENT_NAME,MACS_INDIR,MEME_OUTDIR)
  os.system(SUMMITS_FILTER_CMD)

  # Filter blacklisted regions
  FINAL_SUMMIT_FILTER_CMD = "cat {0}/{3}/{1}_summits_filtered.bed | grep -v -e 'chrM' | sort-bed - | bedops -n 1 - {5}/bin/blacklist/{4}.blacklist.bed > {0}/{3}/{1}_summits_filtered_final.bed".format(OUTDIR,EXPERIMENT_NAME,MACS_INDIR,MEME_OUTDIR,GENOME,CR_PATH)
  os.system(FINAL_SUMMIT_FILTER_CMD)

  # Get 1000 random peaks from top 5000 peaks
  RANDOMIZE = "sort -k5gr {0}{3}/{1}_summits_filtered.bed | head -5000 | shuf | head -1000 | sort-bed -  > {0}/{3}/{1}_1000_random_summits.bed".format(OUTDIR,EXPERIMENT_NAME,MACS_INDIR,MEME_OUTDIR) 
  os.system(RANDOMIZE)
  
  # Get fasta for random 1000 / top 5000
  PICK_1000 = "sort -k5gr {0}/{3}/{1}_summits_filtered.bed | head -1000 | sort-bed -  > {0}/{3}/{1}_1000_random_summits.bed; \
               bedops --range 150 -u {0}/{3}/{1}_1000_random_summits.bed > {0}/{3}/{1}_1000_random_summits_padded.bed; \
               bedtools getfasta -fi {4} -bed {0}/{3}/{1}_1000_random_summits_padded.bed -fo {0}/{3}/{1}_1000_random_summits_padded.fa".format(OUTDIR,EXPERIMENT_NAME,MACS_INDIR,MEME_OUTDIR,GENOME_FA)
  os.system(PICK_1000)
  
  PICK_5000 = "sort -k5gr {0}/{3}/{1}_summits_filtered.bed | head -5000 | sort-bed -  > {0}/{3}/{1}_5000_random_summits.bed; \
               bedops --range 150 -u {0}/{3}/{1}_5000_random_summits.bed > {0}/{3}/{1}_5000_random_summits_padded.bed; \
               bedtools getfasta -fi {4} -bed {0}/{3}/{1}_5000_random_summits_padded.bed -fo {0}/{3}/{1}_5000_random_summits_padded.fa".format(OUTDIR,EXPERIMENT_NAME,MACS_INDIR,MEME_OUTDIR,GENOME_FA)
  os.system(PICK_5000)
  
  PICK_ALL = "bedops --range 150 -u {0}/{3}/{1}_summits_filtered.bed | sort-bed - > {0}/{3}/{1}_all_summits_padded.bed; \
            bedtools getfasta -fi {4} -bed {0}/{3}/{1}_all_summits_padded.bed -fo {0}/{3}/{1}_all_summits_padded.fa".format(OUTDIR,EXPERIMENT_NAME,MACS_INDIR,MEME_OUTDIR,GENOME_FA)
  os.system(PICK_ALL)

  
  sys.stderr.write("*** {0} *** MB_pipeline: Running motif search on {1} using {2} \n".format(datetime.now(),EXPERIMENT_NAME,MEME))    

  MEME_RUN_1000 = "meme-chip -oc {0}/{3}/MEME_1000/ -dreme-m 10 -meme-nmotifs 10 {0}/{3}/{1}_1000_random_summits_padded.fa".format(OUTDIR,EXPERIMENT_NAME,MACS_INDIR,MEME_OUTDIR)
  os.system(MEME_RUN_1000 )

  MEME_RUN_5000 = "meme-chip -oc {0}/{3}/MEME_5000/ -dreme-m 10 -meme-nmotifs 10 {0}/{3}/{1}_5000_random_summits_padded.fa".format(OUTDIR,EXPERIMENT_NAME,MACS_INDIR,MEME_OUTDIR)
  os.system(MEME_RUN_5000)

  MEME_RUN_ALL = "meme-chip -oc {0}/{3}/MEME_ALL/ -dreme-m 10 -meme-nmotifs 10 {0}/{3}/{1}_all_summits_padded.fa".format(OUTDIR,EXPERIMENT_NAME,MACS_INDIR,MEME_OUTDIR)
  os.system(MEME_RUN_ALL)



  #########################
  # Footprinting analysis #
  #########################
  # TODO
  
  sys.stderr.write("*** Pipeline runtime:\n{}\n".format(datetime.now() - startTime))  
  
  
  
  
  
  


main()
