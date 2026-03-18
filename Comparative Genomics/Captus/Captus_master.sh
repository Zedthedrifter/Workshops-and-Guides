#!/bin/bash

function main {

WORKDIR=$SCRATCH/Captus_tutorial #your working directory for this exercise  
METADATA=/mnt/shared/projects/rbge/zedchen/workshop/Captus/metadata_captus.csv
SCRIPTS=/mnt/shared/projects/rbge/zedchen/workshop/Captus
BAITS=$SCRIPTS/Rhodo_baits.fasta


DATA=$SCRATCH/Captus_tutorial/00_raw_reads
RESULT1=$WORKDIR/01_CLEANED
RESULT2=$WORKDIR/02_ASSEMBLY
RESULT3=$WORKDIR/03_EXTRACT
RESULT4=$WORKDIR/04_ALIGN
RESULT5=$WORKDIR/05_IQTREE
RESULT=$WORKDIR/0
#setup 
function setup_dir {

mkdir $RESULT1 -p
mkdir $RESULT2 -p
mkdir $RESULT3 -p
mkdir $RESULT4 -p
mkdir $RESULT5 -p
}



#=========WORKFLOW============================

setup_dir

#CONFIG
$HOME/apps/conda/bin/conda init
source $HOME/apps/conda/bin/activate captus

#STEP 0: cleaning ARRAY
function captus_clean {

JOB1=$(sbatch --array=1-4 --cpus-per-task=4 --mem=4G --job-name="CLEANING" $SCRIPTS/Functions_captus.sh captus_clean $METADATA $DATA $RESULT1  | awk '{print $4}')
echo "Job $JOB1 submitted"
}

#STEP 1: assembly ARRAY
function captus_assemble {

JOB2=$(sbatch --array=1-4 --cpus-per-task=4 --mem=4G --job-name="ASSEMBLY" $SCRIPTS/Functions_captus.sh captus_assembly $METADATA $RESULT1 $RESULT2  | awk '{print $4}')
echo "Job $JOB2 submitted"
}

#STEP 2: EXTRACTION ARRAY
function captus_extract {

JOB2=$(sbatch --array=2,3 \
              --cpus-per-task=4 \
              --mem=16G \
              --job-name="EXTRACT" \
              $SCRIPTS/Functions_captus.sh captus_extract $METADATA $RESULT2 $RESULT3 $BAITS  \
              | awk '{print $4}')

echo "Job $JOB2 submitted"
}

#STEP 3: ALIGNMENT, submit as SINGLE!!!! automatically parallelized by captus
function captus_align {

JOB2=$(sbatch --cpus-per-task=4 \
              --mem=16G \
              --job-name="EXTRACT" \
              $SCRIPTS/Functions_captus.sh captus_align $METADATA $RESULT3 $RESULT4  \
              | awk '{print $4}')

echo "Job $JOB2 submitted"
}

#STEP 4: MAKE TREE
function IQTree {

INDIR=$1
OUTDIR=$2
ARRAY=$3

mkdir $OUTDIR
JOB2=$(sbatch --array=${ARRAY} \
              --cpus-per-task=2 \
              --mem=2G \
              --job-name="iqtree" \
              $SCRIPTS/Functions_captus.sh IQTree $METADATA $INDIR $OUTDIR  \
              | awk '{print $4}')

echo "Job $JOB2 submitted"
}


#=========EXECUTION============================
#captus_clean
#captus_assemble
#captus_extract
#captus_align
IQTree $RESULT4/03_trimmed/06_informed/01_coding_NUC/02_NT $RESULT5/coding_NUC '1-152'
IQTree $RESULT4/03_trimmed/06_informed/02_coding_PTD/02_NT $RESULT5/coding_PTD '1-52'
IQTree $RESULT4/03_trimmed/06_informed/03_coding_MIT/02_NT $RESULT5/coding_MIT '1-34'
}

main

