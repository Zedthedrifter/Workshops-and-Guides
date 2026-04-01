#!/bin/bash

##########################################################################################################

#FUNCTION FASTP CLEAN

function simple_QC {

#SKIP TRIM ADAPTOR
#annotion
INDIR=$1
OUTDIR=$2
SAMPLE=$3

seqkit seq -m 100 \
           --min-qual 20 \
           --remove-gaps \
           ${INDIR}/${SAMPLE}_18S.fastq.gz > ${OUTDIR}/${SAMPLE}_18S.fq

gzip ${OUTDIR}/${SAMPLE}_18S.fq

}

#REORIENT

function v_reorient {

INDIR=$1
OUTDIR=$2
REF=$3
SAMPLE=$4

vsearch --orient ${INDIR}/${SAMPLE}_18S.fastq.gz \
        --db $REF \
        --fastaout ${OUTDIR}/${SAMPLE}_18S.fa
}

# ==================== FUNCTION: DEREPLICATION ====================
# Dereplicates oriented reads, adding size annotations.
function derep {
    INDIR=$1
    OUTDIR=$2
    SAMPLE=$3

    echo "Dereplicating $SAMPLE"
    vsearch --derep_fulllength $INDIR/${SAMPLE}_18S.fa \
            --sizeout --relabel uniq \
            --output $OUTDIR/${SAMPLE}_derep.fa
}

# ==================== FUNCTION: UNOISE3 DENOISING ====================
# Runs cluster_unoise with given minsize to infer ASVs.
function unoise {
    INDIR=$1
    OUTDIR=$2
    MINSIZE=$3
    SAMPLE=$4

    echo "Denoising $SAMPLE with minsize $MINSIZE"
    vsearch --cluster_unoise $INDIR/${SAMPLE}_derep.fa \
            --id 1 \
            --sizein --sizeout \
            --minsize $MINSIZE \
            --consout $OUTDIR/${SAMPLE}_unoise.min${MINSIZE}.fa
}

# ==================== FUNCTION: CHIMERA REMOVAL ====================
# Uses uchime3_denovo to remove chimeric ASVs.
function chimera {
    INDIR=$1
    OUTDIR=$2
    MINSIZE=$3
    SAMPLE=$4

    echo "Removing chimeras for $SAMPLE"
    vsearch --uchime3_denovo $INDIR/${SAMPLE}*.fa \
            --sizein --sizeout \
            --nonchimeras $OUTDIR/${SAMPLE}.fa
}

# ==================== FUNCTION: BUILD ASV TABLE ====================

function compile_index {
    INDIR=$1
    OUTDIR=$2
    
    cat ${INDIR}/*_nochim.min2.fa > $OUTDIR/tmp_asv_sequences.fa
    #remove duplicates
    echo "Dereplicating Compiled reads"
    vsearch --derep_fulllength $OUTDIR/tmp_asv_sequences.fa \
            --sizeout --relabel uniq \
            --output $OUTDIR/tmp_asv_sequences.dedup.fa
    awk '/^>/ {printf ">ASV_%d\n", ++c; next} {print}' $OUTDIR/tmp_asv_sequences.dedup.fa > $OUTDIR/asv_sequences.fa
    # Index the ASV sequences
    vsearch --usearch_global $OUTDIR/asv_sequences.fa \
            --db $OUTDIR/asv_sequences.fa \
            --self --id 1.0 --userout $OUTDIR/asv_index.txt \
            --userfields query+target

}

# Creates a count table of ASVs per sample.
function asv_counts {

    INDIR=$1
    OUTDIR=$2
    MINSIZE=$3
    SAMPLE=$4
    # Map each sample's non-chimeric reads to the ASV set
    vsearch --usearch_global $INDIR/${SAMPLE}_nochim.min${MINSIZE}.fa \
            --db $OUTDIR/asv_sequences.fa \
            --strand plus \
            --id 1 \
            --otutabout $OUTDIR/${SAMPLE}_tmp.txt
    sort -k1,1V -nr $OUTDIR/${SAMPLE}_tmp.txt > $OUTDIR/${SAMPLE}_asv_counts.txt
    rm -f $OUTDIR/${SAMPLE}_tmp.txt
}

# ==================== FUNCTION: TAXONOMY ASSIGNMENT ====================
# Uses sintax to assign taxonomy to the final ASV set.
# This function processes all samples together (not array).
function taxonomy {
    INFILE=$1
    REF=$2

    # Assign taxonomy
    vsearch --sintax $INFILE \
            --db $REF \
            --tabbedout ${INFILE/fa/tsv} \
            --sintax_cutoff 0.8 \
            --threads 8
    cat ${INFILE/fa/tsv}|sort -k1,1V -nr > ${INFILE/fa/txt}
}

function compile_table {

    METADATA=$1
    INDIR=$2
    TAXON=$3
    
    #get heading
    sed -i 's/\r$//' $TAXON #after manually made the file in excel
    echo 'INDEX' > $INDIR/tmp.heading.txt
    cat $TAXON|tail -n +2 |sort -k1,1V -nr >> $INDIR/tmp.heading.txt
    #get entry from each sample
    for SAMPLE in $(cat ${METADATA} |tail +2| cut -d ',' -f4)
      do
        echo ${SAMPLE} > $INDIR/tmp.${SAMPLE}_clm2.txt #collect sample name
        cat $INDIR/${SAMPLE}_asv_counts.txt|sort -k1,1V -nr |head -n -1| cut -f 2  >> $INDIR/tmp.${SAMPLE}_clm2.txt
      done
    paste -d ',' $TAXON $INDIR/tmp.*_clm2.txt > $INDIR/asv_counts_compiled.csv
    
    rm $INDIR/tmp.heading.txt $INDIR/tmp.*_clm2.txt -f
    
    #extract just the ASVs with chlorophyta
    echo 'INDEX' > Chlorophyta_ASV.txt
    grep 'Chlorophyta' $TAXON |cut -f 1|sort -k1,1V -nr >> Chlorophyta_ASV.txt
    echo "extracting $(wc -l Chlorophyta_ASV.txt) chlorophyta containing ASVs" 
    grep -wFf Chlorophyta_ASV.txt $INDIR/asv_counts_compiled.csv > $INDIR/asv_counts.chloro.csv
    
    
}

#ALTERNATIVE UNOISE

function com_UNOISE {

    INDIR=$1
    OUTDIR=$2
    MINSIZE=$3
    
    cat ${INDIR}/*.fa > $OUTDIR/tmp_asv_sequences.fa
    #remove duplicates
    echo "Dereplicating Compiled reads"
    vsearch --derep_fulllength $OUTDIR/tmp_asv_sequences.fa \
            --sizeout --relabel uniq \
            --output $OUTDIR/tmp_asv_sequences.dedup.fa
    echo "done derep"
    #UNOISE
    vsearch --cluster_unoise $OUTDIR/tmp_asv_sequences.dedup.fa \
            --id 1 \
            --sizein --sizeout \
            --minsize $MINSIZE \
            --consout $OUTDIR/tmp_asv_sequences.unoise.fa
    awk '/^>/ {printf ">ASV_%d\n", ++c; next} {print}' $OUTDIR/tmp_asv_sequences.unoise.fa > $OUTDIR/asv_sequences.fa
    echo "done unoise"
    # Index the ASV sequences
    vsearch --usearch_global $OUTDIR/asv_sequences.fa \
            --db $OUTDIR/asv_sequences.fa \
            --self --id 1.0 --userout $OUTDIR/asv_index.txt \
            --userfields query+target
    echo "done index"
    grep '>' $OUTDIR/asv_sequences.fa|cut -d '>' -f 2 > $OUTDIR/asv_keep.txt
    
    rm $OUTDIR/tmp_asv_sequences*fa
    

}


#================================================================================
#USER INPUTS
MY_FUNCTION=$1 #YOU don't need this line. just a reminder of what the first variable means
METADATA=$2
#GET SAMPLE
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${METADATA} | cut -d ',' -f4)
echo $SAMPLE
# Check command-line argument and call the function
case "$1" in
    simple_QC)
        simple_QC $3 $4 $SAMPLE # 
        ;;
    v_reorient)
        v_reorient $3 $4 $5 $SAMPLE # 
        ;;
    derep )
        derep  $3 $4 $SAMPLE # 
        ;;
    unoise )
        unoise  $3 $4 $5 $SAMPLE # 
        ;;
    chimera )
        chimera  $3 $4 $5 $SAMPLE # 
        ;;
    compile_index )
        compile_index  $3 $4 # 
        ;;
    asv_counts )
        asv_counts  $3 $4 $5 $SAMPLE # 
        ;;
    compile_table )
        compile_table $2 $3 $4 # 
        ;;
    taxonomy )
        taxonomy  $3 $4 # 
        ;;
    com_UNOISE)
        com_UNOISE $2 $3 $4 # 
        ;;
    *)
        echo "Unknown function: $1: Usage: $0 {captus_clean|captus_assembly|fastp_qc} [args]"
        exit 1
        ;;
esac