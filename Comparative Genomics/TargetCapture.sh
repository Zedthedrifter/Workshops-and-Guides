#!/bin/bash
#SBATCH --job-name=my_array_job
#SBATCH --array=4-5       # Run tasks 1 through 10
#SBATCH --cpus-per-task=4     # Request 4 CPUs per task
#SBATCH --mem=8G            # we don't need much memory for this

function setup_conda {

conda create  --name target_capture
conda install --name target_capture -y fastp
conda install --name target_capture -y bwa
conda install --name target_capture -y samtools
conda install --name target_capture -y bcftools
conda install --name target_capture -y iqtree
conda install --name target_capture -y seqkit
conda install --name target_capture -y blast
}

function download_ena {

OUTDIR=$1

#wget -nc -P $OUTDIR ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR776/ERR776786/ERR776786_1.fastq.gz
#wget -nc -P $OUTDIR ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR776/ERR776786/ERR776786_2.fastq.gz 
#wget -nc -P $OUTDIR ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR776/ERR776790/ERR776790_1.fastq.gz
#wget -nc -P $OUTDIR ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR776/ERR776790/ERR776790_2.fastq.gz
#wget -nc -P $OUTDIR ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR776/ERR776780/ERR776780_1.fastq.gz
#wget -nc -P $OUTDIR ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR776/ERR776780/ERR776780_2.fastq.gz
wget -nc -P $OUTDIR ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR776/ERR776792/ERR776792_1.fastq.gz
wget -nc -P $OUTDIR ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR776/ERR776822/ERR776822_1.fastq.gz
wget -nc -P $OUTDIR ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR776/ERR776822/ERR776822_2.fastq.gz
wget -nc -P $OUTDIR ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR776/ERR776792/ERR776792_2.fastq.gz

}

#CLEANING UP READS AND QC REPORT

function fastp_qc {

#SKIP TRIM ADAPTOR
#annotion 
INDIR=$1
OUTDIR=$2
SAMPLE=$3

fastp -i ${INDIR}/${SAMPLE}_1.fastq.gz -I ${INDIR}/${SAMPLE}_2.fastq.gz \
      -o ${OUTDIR}/${SAMPLE}_R1.fq.gz   -O ${OUTDIR}/${SAMPLE}_R2.fq.gz \
      --unpaired1 ${OUTDIR}/${SAMPLE}_unpaired1.fq.gz \
      --unpaired2 ${OUTDIR}/${SAMPLE}_unpaired2.fq.gz \
      -h ${OUTDIR}/${SAMPLE}.html \
      --json ${OUTDIR}/${SAMPLE}.json \
      --detect_adapter_for_pe \
      --qualified_quality_phred 20 \
      --unqualified_percent_limit 40\
      -c \
      --cut_front \
      --cut_front_mean_quality 20 \
      --cut_tail \
      --cut_front_mean_quality 20 \
      --cut_right \
      --cut_right_window_size 4 \
      --cut_right_mean_quality 20 \
      --length_required 36 \
      --trim_poly_g \
      --thread 8 \
      > ${OUTDIR}/${SAMPLE}.fastp.log

rm -f ${OUTDIR}/${SAMPLE}_unpaired*.fq.gz ${OUTDIR}/${SAMPLE}.json 
# Enable 5' trimming (LEADING equivalent)
# Quality threshold for 5' trim
# Enable 3' trimming (TRAILING equivalent)
# Quality threshold for 3' trim  
# Enable sliding window trim (SLIDINGWINDOW equivalent)
# Window size = 4
# Mean quality threshold = 15
# MINLEN equivalent
#--trim_poly_g                    force polyG tail trimming, by default trimming is automatically enabled for Illumina NextSeq/NovaSeq data
}

#MAPPING
#BOWTIE2 MAPPING PIPELINE
function bowtie2_map {

INDIR=$1
IDX=$2
OUTDIR=$3
SAMPLE=$4

minQ=10

echo 'RUNNING BOWTIE2 MAPPING ON' $SAMPLE
bowtie2-align-s --wrapper basic-0 \
                -x $IDX \
                -1 ${INDIR}/${SAMPLE}_R1.fq.gz \
                -2 ${INDIR}/${SAMPLE}_R2.fq.gz \
                -p 20 \
                -N 1 \
                -L 20 \
                --threads 8 \
                --phred33 \
                --very-sensitive-local \
                --no-discordant \
                --no-mixed \
                --no-unal \
                --time \
                --rg-id ${SAMPLE} \
                --rg SM:${SAMPLE} \
                --rg PL:'ILLUMINA' |\
                samtools view -Sbh -F 4 -@ 8 -o $OUTDIR/${SAMPLE}.all.bam #all mapped reads
#number of all mapped reads
#samtools view -F 4 -c $OUTDIR/${prefix}_${SLURM_ARRAY_TASK_ID}.all.bam > $log

#only paired mapped reads with q>minQ
samtools view $OUTDIR/${SAMPLE}.all.bam -Sbh -F 4 -f 3 -q $minQ -@ 8 |samtools sort  -@ 8 -o $OUTDIR/${SAMPLE}.sorted.bam

echo 'RUNNING SAMTOOLS INDEXING ON' ${SAMPLE}
samtools index $OUTDIR/${SAMPLE}.sorted.bam
}

#BWA MAPPING PIPELINE
function bwa_map {

INDIR=$1
REF=$2
OUTDIR=$3
SAMPLE=$4

minQ=10

echo 'RUNNING BWA MAPPING ON' ${SAMPLE}
bwa mem $REF \
        ${INDIR}/${SAMPLE}_R1.fq.gz \
        ${INDIR}/${SAMPLE}_R2.fq.gz \
        -t 8 \
        -k 20 |\
        samtools view -Sbh -F 4 -@ 8 -o $OUTDIR/${SAMPLE}.all.bam

#only paired mapped reads with q> minQ
#add heading manually to bwa outputs
  samtools view $OUTDIR/${SAMPLE}.all.bam -Sbh -F 4 -f 3 -q $minQ -@ 8 \
 |samtools addreplacerg \
         -r "ID:${SAMPLE}" \
         -r "SM:${SAMPLE}" \
         -r "LB:${SAMPLE}" \
         -r "PL:ILLUMINA" \
         -O bam \
         -o - \
         - \
 | samtools sort  -@ 8 -o $OUTDIR/${SAMPLE}.sorted.bam

#INDEX
echo 'RUNNING SAMTOOLS INDEXING ON' ${SAMPLE}
samtools index $OUTDIR/${SAMPLE}.sorted.bam
#rm $OUTDIR/${SAMPLE}.all.bam
}

#SNP calling : 
function snp_calling {

INDIR=$1
OUTDIR=$2
REF=$3
SAMPLE=$4

#FILTERING
minQ=20
minDP=10
MIN_FMT_DP=5
MIN_MINOR_ALLELE_FREQ=0.01
MIN_MAPPING_QUALITY=30


echo 'calling SNPs'
#ls $INDIR/${prefix}*sorted.bam > bamlist.txt #need to change back to sorted.bam
#can ues a --bam-list $bamlist optio if you're doing all samples at the same time
bcftools mpileup -Ou -f $REF $INDIR/${SAMPLE}.sorted.bam --threads 8 \
                 --annotate INFO/AD,FORMAT/DP,FORMAT/AD \
 |   bcftools call -Ou -mv \
 |   bcftools view -i '1==1' 2>/dev/null \
 |   bcftools filter -e "QUAL<${minQ} || INFO/DP<${minDP}" \
 |   bcftools +fill-tags -- -t AC,AN,AF,MAF,HWE  \
 |   bcftools filter   -e "INFO/MAF[0] < ${MIN_MINOR_ALLELE_FREQ}"  \
 |   bcftools filter   -e "INFO/MQ < ${MIN_MAPPING_QUALITY}"  \
 |   bcftools filter   -e 'INFO/VDB < 0.1'  -Oz -o $OUTDIR/${SAMPLE}.flt3.vcf.gz

echo "Indexing filtered VCF..."
    bcftools index $OUTDIR/${SAMPLE}.flt3.vcf.gz -f
}

#MAKE CONSENSUS FROM VCF
#THIS ONE NEEDS MORE MODIFICATION. WE CAN WORRY ABOUT IT LATER
function extract_consensus { 

BAMINDIR=$1
VCFINDIR=$2
OUTDIR=$3
REF=$4
SAMPLE=$5

echo "EXTRACT CONSENSUS SEQUENCE FROM $REF REF AND $SAMPLE"

MIN_COVERAGE=2
BAM_FILE=$BAMINDIR/${SAMPLE}.sorted.bam
#FOR LOOP OVER ALL THE GENES IN THE BAIT SET

mkdir ${OUTDIR}/${SAMPLE}
for HEADING in $(grep '>' $REF)

do

  CHRM=${HEADING/>/};echo "Processing $CHRM ..."
  mkdir ${OUTDIR}/${CHRM}_${SAMPLE}_tmp #a temporary directory so the files won't mix and interfere with each other
  #MARK LOW COVERAGE AREA
  samtools depth -a -r $CHRM $BAM_FILE | awk -v min="$MIN_COVERAGE" '$3 < min {print $1"\t"$2-1"\t"$2}'  > ${OUTDIR}/${CHRM}_${SAMPLE}_tmp/low_coverage.bed
   
  # Generate consensus for this sample
  #mask the regions with low coverage and get consensus for one sample
  samtools faidx $REF $CHRM | \
  bcftools consensus -s ${SAMPLE} $VCFINDIR/${SAMPLE}.flt3.vcf.gz \
                     -I \
                     --mask ${OUTDIR}/${CHRM}_${SAMPLE}_tmp/low_coverage.bed \
                     -o ${OUTDIR}/${CHRM}_${SAMPLE}_tmp/${SAMPLE}.fasta
          
  # Add sample name to header
  sed -i "1s/^>.*/>${SAMPLE}/" ${OUTDIR}/${CHRM}_${SAMPLE}_tmp/${SAMPLE}.fasta
  #collect the results
  cp ${OUTDIR}/${CHRM}_${SAMPLE}_tmp/${SAMPLE}.fasta ${OUTDIR}/${SAMPLE}/${CHRM}.fasta
  #remove temporary files
  rm -rf ${OUTDIR}/${CHRM}_${SAMPLE}_tmp
done
}

#COMPILE GENE AND ALIGN
function easy353_mafft {

INDIR=$1
OUTDIR=$2
GENE=$3

echo "compiling and aligning ${GENE}"
cat $INDIR/*/${GENE}.fasta > ${OUTDIR}/tmp_${GENE}.fasta
mafft --maxiterate 10000 $OUTDIR/tmp_${GENE}.fasta > $OUTDIR/${GENE}.fasta
rm -f $OUTDIR/tmp_${GENE}.fasta
#remove any space that might cause trouble later
sed -i 's/ //g' $OUTDIR/${GENE}.fasta
}

#MAKE PHYLOGENY
function iqtree_per_gene {

INDIR=$1
OUTDIR=$2
GENE=$3

mkdir $OUTDIR/$GENE
cp $INDIR/${GENE}.fasta $OUTDIR/$GENE
iqtree -s $OUTDIR/$GENE/${GENE}.fasta -bb 1000 -redo -safe
mv $OUTDIR/$GENE/*treefile $OUTDIR
rm -rf $OUTDIR/$GENE


}
#=====================================================================================================

function main {

WORKDIR=$SCRATCH/Target_capture #your working directory for this exercise  
METADATA=$HOME/projects/rbge/atam/Target_Capture_tutorial/metadata.csv
BAITS=$WORKDIR/CK_GIT/Ref.fna
REF=$WORKDIR/References #YOU CAN BUILD A REFERENCE DIR IN PROJECT. WE WILL TALK ABOUT THAT LATER


DATA=$WORKDIR/00_DATA
RESULT1=$WORKDIR/01_CLEANED
RESULT2=$WORKDIR/02_MAPPING
RESULT3=$WORKDIR/03_VCF
RESULT4=$WORKDIR/04_EXTRACTS
RESULT5=$WORKDIR/05_ALIGNMENTS
RESULT6=$WORKDIR/06_TREES
#setup 
function setup_dir {

mkdir $WORKDIR/CK_GIT -p
mkdir $REF -p
mkdir $DATA -p
mkdir $RESULT1 -p
mkdir $RESULT2 -p
mkdir $RESULT3 -p
mkdir $RESULT4 -p
mkdir $RESULT5 -p
mkdir $RESULT6 -p
}

setup_dir
#===========================================================================================
#EXECUTE 

#setup_conda
#download_ena $DATA
#GET REFERENCE 
#git clone https://github.com/ckidner/Targeted_enrichment.git $WORKDIR/CK_GIT


#ANALYSIS

#SAMPLE ARRAY: 
SLURM_ARRAY_TASK_ID=3
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${METADATA} | cut -d ',' -f1)
echo $SAMPLE

#STEP 1: CLEANING UP READS AND QC REPORT; SAMPLE ARRAY
#fastp_qc $DATA $RESULT1 $SAMPLE

#STEP 2: MAPPING: SAMPLE ARRAY
#make index
#bwa index $BAITS
#bwa_map $RESULT1 $BAITS $RESULT2 $SAMPLE

#STEP 3: SNP CALLING: SAMPLE ARRAY
snp_calling $RESULT2 $RESULT3 $BAITS $SAMPLE

#STEP 4: CONSENSUS; SAMPLE ARRAY
#extract_consensus  $RESULT2 $RESULT3 $RESULT4 $BAITS $SAMPLE

#===========================================================================================
#GENE ARRAY
GENE=$(grep '>' $BAITS|sed -n "${SLURM_ARRAY_TASK_ID}p"|cut -d '>' -f 2) #get the Nth gene in the bait set (fasta), remove > from heading
echo $GENE
#GENE=comp41081_c0_seq1

#STEP 5: ALIGNMENT; GENE ARRAY
#easy353_mafft $RESULT4 $RESULT5 $GENE

#STEP 6: PHYLOGENY: GENE ARRAY
#iqtree_per_gene $RESULT5 $RESULT6 $GENE
}

main

