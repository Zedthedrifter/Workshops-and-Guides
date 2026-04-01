#!/bin/bash

function main {

#conda create -n captus -c bioconda -c conda-forge captus iqtree -y# T

#UNIVERSAL VARIABLES
WORKDIR=$SCRATCH/Lichen_18S_metabarcoding
SCRIPTS=/mnt/shared/projects/rbge/zedchen/Lichens/18S_metabarcode
METADATA=$SCRIPTS/Eukaryome_samples.csv
PR2=$WORKDIR/REF/pr2_version_5.1.1_SSU_UTAX.fasta.gz

#CONSTANTS
sratools=$HOME/apps/manual/sratoolkit.3.2.1-ubuntu64/bin/
pip=$HOME/apps/env/easy353/bin/pip
STAR=$HOME/apps/RNA_polish/bin/STAR/bin/Linux_x86_64/STAR

#subdirectories
DATA=$SCRATCH/Lichen_18S_metabarcoding/SSU_fastq_links
RESULT1=$WORKDIR/01_clean_reads
RESULT2=$WORKDIR/02_REORIENT
RESULT3=$WORKDIR/03_DEREP
RESULT4=$WORKDIR/04_UNOISE
RESULT5=$WORKDIR/05_NOCHIM
RESULT6=$WORKDIR/06_AVS_SUM
#COMPILE FROM DEREP STAGE AND UNOISE WITH COMPILED READS
RESULT7=$WORKDIR/07_CO_UNOISE
RESULT8=$WORKDIR/08_AVS_ALT
#######################
function setup_dir {

mkdir -p $WORKDIR/REF
mkdir $RESULT1
mkdir $RESULT2
mkdir $RESULT3
mkdir $RESULT4
mkdir $RESULT5
mkdir $RESULT6
mkdir $RESULT7
mkdir $RESULT8
}
setup_dir

#STEP_0 DOWNLOAD
wget -nc -P $WORKDIR/REF https://github.com/pr2database/pr2database/releases/download/v5.1.1/pr2_version_5.1.1_SSU_UTAX.fasta.gz

#CONFIG
$HOME/apps/conda/bin/conda init


#WORKFLOW ========================================================
#STEP 0: cleaning ARRAY
#CAPTUS CLEAN
#ALTERNATIVE: FASTP CLEAN
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

function unoise {

source $HOME/apps/conda/bin/activate metagenomics

minsize=2

JOB=$(sbatch --array=1-53,55-67 \
              --cpus-per-task=2 \
              --mem=4G \
              --job-name="unnoise" \
              $SCRIPTS/meta_functions.sh unoise $METADATA $RESULT3 $RESULT4 $minsize \
              | awk '{print $4}')
echo "Job $JOB submitted"

#these two takes more mem
JOB=$(sbatch --array=54 \
              --cpus-per-task=2 \
              --mem=10G \
              --job-name="unnoise" \
              $SCRIPTS/meta_functions.sh unoise $METADATA $RESULT3 $RESULT4 $minsize \
              | awk '{print $4}')
echo "Job $JOB submitted"

}

function chimera {

source $HOME/apps/conda/bin/activate metagenomics

minsize=2

JOB=$(sbatch --array=1-67 \
              --cpus-per-task=2 \
              --mem=1G \
              --job-name="chimera" \
              $SCRIPTS/meta_functions.sh chimera $METADATA $RESULT4 $RESULT5 $minsize \
              | awk '{print $4}')
echo "Job $JOB submitted"

}

function compile_index {  #ONLY ON COMPILED

source $HOME/apps/conda/bin/activate metagenomics


#68 is all the samples compiled
JOB=$(sbatch  --cpus-per-task=8 \
              --mem=2G \
              --job-name="compile_index" \
              $SCRIPTS/meta_functions.sh compile_index $METADATA $RESULT5 $RESULT6 \
              | awk '{print $4}')
echo "Job $JOB submitted"

}

function taxonomy { #ONLY ON COMPILED

source $HOME/apps/conda/bin/activate metagenomics

minsize=2

#68 is all the samples compiled
JOB=$(sbatch  --cpus-per-task=4 \
              --mem=2G \
              --job-name="taxonomy" \
              $SCRIPTS/meta_functions.sh taxonomy $METADATA $RESULT6/asv_sequences.fa $PR2 \
              | awk '{print $4}')
echo "Job $JOB submitted"

}

function asv_counts {

source $HOME/apps/conda/bin/activate metagenomics

minsize=2

#68 is all the samples compiled
JOB=$(sbatch --array=2-53,55,57-67 \
              --cpus-per-task=2 \
              --mem=4G \
              --job-name="asv_counts" \
              $SCRIPTS/meta_functions.sh asv_counts $METADATA $RESULT5 $RESULT6 $minsize \
              | awk '{print $4}')
JOB=$(sbatch --array=54,56 \
              --cpus-per-task=2 \
              --mem=6G \
              --job-name="asv_counts" \
              $SCRIPTS/meta_functions.sh asv_counts $METADATA $RESULT5 $RESULT6 $minsize \
              | awk '{print $4}')
echo "Job $JOB submitted"

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
#simple_QC #clean reads and trim lowQ regions with slide windows, disable adapter trimmer!!! otherwise overrepresented sequences (our targets) will be removed# omit this step
#v_reorient #> 99.9% sequences oriented. remove the wrongly oriented sequences
#cat /mnt/shared/scratch/zchen/Lichen_18S_metabarcoding/02_REORIENT/*.fa > /mnt/shared/scratch/zchen/Lichen_18S_metabarcoding/02_REORIENT/compiled_18S.fa #also put compiled as the last entry on sample csv: 68
#derep
#chimera ##68 is the compiled reads from all sample. by running #68 we also get a nonchim file for all the ASV. use that one as 'all'
#================
#compile_index # unoise with the collective set, then unoise for individual
#unoise  #use minsize of 2 to 8 to denoise
#taxonomy
#================
#asv_counts
#compile_table #grep -wFf to further filter down samples
}

main