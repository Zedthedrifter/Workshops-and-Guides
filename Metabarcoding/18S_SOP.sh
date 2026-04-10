#!/bin/bash

function main {

#conda create -n metagenomics -c bioconda -c conda-forge vsearch  -y # or install vsearch into one of your existent env: conda install vsearch -y

#UNIVERSAL VARIABLES
EMAIL=YOUR_EMAIL@rbge.org.uk
WORKDIR=$SCRATCH/Lichen_18S_metabarcoding
SCRIPTS=$YOUR_SCRIPT_DIR_ON_PROJECT
METADATA=$SCRIPTS/Eukaryome_samples.csv #COPIED FROM MY WORKSHOP DIR

#SHARED WITH ZED
DATA=/mnt/shared/projects/rbge/zedchen/workshop/18S_metabarcoding/18S_PacBio
PR2=/mnt/shared/projects/rbge/zedchen/workshop/18S_metabarcoding/REF/pr2_version_5.1.1_SSU_UTAX.fasta.gz #or if you downloaded it to your own directory $WORKDIR/REF

#CONSTANTS
sratools=$HOME/apps/manual/sratoolkit.3.2.1-ubuntu64/bin/
pip=$HOME/apps/env/easy353/bin/pip
STAR=$HOME/apps/RNA_polish/bin/STAR/bin/Linux_x86_64/STAR

#subdirectories
RESULT2=$WORKDIR/02_REORIENT #shared with Master_18S
RESULT3=$WORKDIR/03_DEREP #shared with Master_18S
RESULT4=$WORKDIR/04_UNOISE
RESULT5=$WORKDIR/05_NOCHIM
RESULT6=$WORKDIR/06_AVS_SUM

#######################
function setup_dir {

mkdir -p $WORKDIR/REF
mkdir $RESULT1
mkdir $RESULT2
mkdir $RESULT3
mkdir $RESULT4
mkdir $RESULT5
mkdir $RESULT6
}
setup_dir

#STEP_0 DOWNLOAD
#wget -nc -P $WORKDIR/REF https://github.com/pr2database/pr2database/releases/download/v5.1.1/pr2_version_5.1.1_SSU_UTAX.fasta.gz

#CONFIG
$HOME/apps/conda/bin/conda init

#might not need this step 
function simple_QC {

source $HOME/apps/conda/bin/activate metagenomics

JOB=$(sbatch --array=2-67 \
              --cpus-per-task=2 \
              --mem=1G \
              --job-name="CLEANING" \
              $SCRIPTS/meta_functions.sh simple_QC $METADATA $DATA $RESULT1  \
              | awk '{print $4}')
echo "Job $JOB submitted"
}

function v_reorient {

source $HOME/apps/conda/bin/activate metagenomics

JOB=$(sbatch --array=2-67 \
              --cpus-per-task=2 \
              --mem=2G \
              --job-name="REORIENT" \
              $SCRIPTS/meta_functions.sh v_reorient $METADATA $DATA $RESULT2 $PR2 \
              | awk '{print $4}')
echo "Job $JOB submitted"
}

function read_count_REORIENT {


for SAMPLE in $(cat ${METADATA} | cut -d ',' -f4)
do
grep -c '>' $RESULT2/${SAMPLE}_18S.fa
done
}

function derep {

source $HOME/apps/conda/bin/activate metagenomics

JOB=$(sbatch --array=2-67 \
              --cpus-per-task=2 \
              --mem=200M \
              --job-name="derep" \
              $SCRIPTS/meta_functions.sh derep $METADATA $RESULT2 $RESULT3 \
              | awk '{print $4}')
echo "Job $JOB submitted"
}

#CONFIRM READ COUNT DIDN'T CHANGE
function read_count_WITH_SIZE {

INDIR=$1

echo $INDIR
for SAMPLE in $(tail -n +2 ${METADATA} | cut -d ',' -f4)
do
grep '>' $INDIR/${SAMPLE}*.fa|rev|cut -d ';' -f 1|rev|cut -d '=' -f 2|awk '{sum+=$1} END {print sum}'
done
}

function unoise {

source $HOME/apps/conda/bin/activate metagenomics

minsize=1

JOB=$(sbatch --array=2-53,55-67 \
              --cpus-per-task=2 \
              --mem=12G \
              --mail-user=$EMAIL \
              --mail-type=END,FAIL \
              --job-name="unnoise" \
              $SCRIPTS/meta_functions.sh unoise $METADATA $RESULT3 $RESULT4 $minsize \
              | awk '{print $4}')
echo "Job $JOB submitted"

#these two takes more mem
JOB=$(sbatch --array=54 \
              --cpus-per-task=2 \
              --mem=12G \
              --mail-user=$EMAIL \
              --mail-type=END,FAIL \
              --job-name="unnoise" \
              $SCRIPTS/meta_functions.sh unoise $METADATA $RESULT3 $RESULT4 $minsize \
              | awk '{print $4}')
echo "Job $JOB submitted"

}

function chimera {

source $HOME/apps/conda/bin/activate metagenomics

minsize=1

JOB=$(sbatch --array=2-67 \
              --cpus-per-task=2 \
              --mem=2G \
              --job-name="chimera" \
              --mail-user=$EMAIL \
              --mail-type=END,FAIL \
              $SCRIPTS/meta_functions.sh chimera $METADATA $RESULT4 $RESULT5 $minsize \
              | awk '{print $4}')
echo "Job $JOB submitted"

}

#COMPILE ALL THE READS AND UNOISE AND CLASSIFY
function com_UNOISE {  #ONLY ON COMPILED

source $HOME/apps/conda/bin/activate metagenomics

minsize=1 #minsize=4
#68 is all the samples compiled
JOB=$(sbatch  --cpus-per-task=4 \
              --mem=10G \
              --job-name="com_UNOISE" \
              --mail-user=$EMAIL \
              --mail-type=END,FAIL \
              $SCRIPTS/meta_functions.sh com_UNOISE $METADATA $RESULT5 $RESULT6 $minsize \
              | awk '{print $4}')
echo "Job $JOB submitted"

}


function taxonomy { #ONLY ON COMPILED

source $HOME/apps/conda/bin/activate metagenomics

JOB=$(sbatch  --cpus-per-task=4 \
              --mem=2G \
              --job-name="taxonomy" \
              $SCRIPTS/meta_functions.sh taxonomy $METADATA $RESULT6/asv_sequences.fa $PR2 \
              | awk '{print $4}')
echo "Job $JOB submitted"

}

###
function asv_counts {

source $HOME/apps/conda/bin/activate metagenomics

minsize=1

#USE THE PRECOMPILED/UNOISED-TOGETHER FILE
JOB=$(sbatch --array=2-67 \
              --cpus-per-task=2 \
              --mem=16G \
              --job-name="asv_counts" \
              $SCRIPTS/meta_functions.sh asv_counts $METADATA $RESULT5 $RESULT6 $minsize \
              | awk '{print $4}')

echo "Job $JOB submitted"
}

#CONFIRM READ COUNT DIDN'T CHANGE
function read_count_ASV_COUNT {

INDIR=$1

echo $INDIR
for SAMPLE in $(tail -n +2 ${METADATA} | cut -d ',' -f4)
do
cat $INDIR/${SAMPLE}_asv_counts.txt|cut -f 2|awk '{sum+=$1} END {print sum}'
done
}

function compile_table { #ONLY ON COMPILED

#68 is all the samples compiled
JOB=$(sbatch  --cpus-per-task=1 \
              --mem=2K \
              --job-name="compile_table" \
              $SCRIPTS/meta_functions.sh compile_table $METADATA $RESULT6 $RESULT6/asv_sequences.csv \
              | awk '{print $4}')
echo "Job $JOB submitted"

}


#EXECUTE =============================================================
#v_reorient #> 99.9% sequences oriented. remove the wrongly oriented sequences
#read_count_REORIENT 

#derep
#read_count_WITH_SIZE $RESULT3

#unoise  #use minsize =1 SHOULDN'T REMOVE ANYTHING
#read_count_WITH_SIZE $RESULT4

#chimera 
#read_count_WITH_SIZE $RESULT5

#==COMPILE, CLASSIFY, INDEX==============
#com_UNOISE
#taxonomy
#NOW MAKEK asv_sequences.csv TO SAVE YOURSELF SOME TIME (OR JUST USE MINE, SHOULD BE THE SAME)

#================
#asv_counts
#read_count_ASV_COUNT $RESULT6
#compile_table #grep -wFf to further filter down samples
$SCRIPTS/meta_functions.sh read_count_chloro $METADATA $RESULT6/asv_counts.chloro.csv
}

main