#!/bin/bash

#USER INPUTS
METADATA=$2
#GET SAMPLE
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${METADATA} | cut -d ',' -f1)
echo $SAMPLE

##########################################################################################################
function captus_clean {

INDIR=$1
OUTDIR=$2
SAMPLE=$3

#ANALYSIS
ls $INDIR/${SAMPLE}_R1.fq.gz
captus clean --reads $INDIR/${SAMPLE}_R1.fq.gz $INDIR/${SAMPLE}_R2.fq.gz \
             --out $OUTDIR \
             --overwrite
}

#ASSEMBLY
function captus_assembly {

INDIR=$1
OUTDIR=$2
SAMPLE=$3

#ANALYSIS
ls $INDIR/${SAMPLE}_R1.fq.gz
captus assemble --reads $INDIR/${SAMPLE}_R1.fq.gz $INDIR/${SAMPLE}_R2.fq.gz \
                --out $OUTDIR \
                --k_list 55,65 \
                --min_count 3 \
                --prune_level 2 \
                --tmp_dir $OUTDIR/${SAMPLE} \
                --overwrite

rm -rf $OUTDIR/${SAMPLE}
#it's important to keep temporary folders independent of each other when you run things in parallel
}

#EXTRACTION 
function captus_extract {

INDIR=$1
OUTDIR=$2
BAITS=$3
SAMPLE=$4

cp $INDIR/${SAMPLE}__captus-asm/01_assembly/assembly.fasta $INDIR/${SAMPLE}__captus-asm/01_assembly/${SAMPLE}.fasta
captus extract -a $OUTDIR/tmp_${SAMPLE} \
               --fastas $INDIR/${SAMPLE}__captus-asm/01_assembly/${SAMPLE}.fasta \
               --out $OUTDIR/ \
               --nuc_refs $BAITS \
               --ptd_refs SeedPlantsPTD \
               --mit_refs SeedPlantsMIT \
               --nuc_min_identity 50 \
               --overwrite
               
rm -rf $OUTDIR/tmp_${SAMPLE} $INDIR/${SAMPLE}__captus-asm/01_assembly/${SAMPLE}.fasta
}

#ALIGNMENT
function captus_align {

INDIR=$1
OUTDIR=$2

captus align --captus_extractions_dir $INDIR \
             --out $OUTDIR \
             --min_samples 4 \
             --max_paralogs 5 \
             --align_method mafft_auto \
             --overwrite

#for mafft configuration error:
#conda remove mafft --force
#conda install -c bioconda mafft
}

function IQTree {

INDIR=$1
OUTDIR=$2

FASTA=$(ls $INDIR|sed -n "${SLURM_ARRAY_TASK_ID}p"|cut -d '>' -f 2) #get the Nth gene in the bait set (fasta), remove > from heading
echo $FASTA
GENE=${FASTA/.fna/}

mkdir $OUTDIR/tmp_${GENE}
cp $INDIR/$FASTA $OUTDIR/tmp_${GENE}
iqtree -s $OUTDIR/tmp_${GENE}/$FASTA -bb 1000 -redo -safe
mv $OUTDIR/tmp_${GENE}/*treefile $OUTDIR/${GENE}.treefile
rm -rf $OUTDIR/tmp_${GENE}
}
###############################################################################
# Check command-line argument and call the function
case "$1" in
    captus_clean)
        captus_clean $3 $4 $SAMPLE # 
        ;;
    captus_assembly)
        captus_assembly $3 $4 $SAMPLE # 
        ;;
    captus_extract)
        captus_extract $3 $4 $5 $SAMPLE # 
        ;;
    captus_align)
        captus_align $3 $4 # 
        ;;
    IQTree)
        IQTree $3 $4 # 
        ;;
    *)
        echo "Usage: $0 {captus_clean|captus_assembly|captus_extract|IQTree} [args]"
        exit 1
        ;;
esac
