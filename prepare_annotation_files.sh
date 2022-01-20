#!/bin/bash

## this script use the reference genome to generate all the annotation files used for raw data processing
## reference genome (fasta file) and the genome annotation (gtf file) should be pre-downloaded.
## (optional) filter the gtf file to protein coding genes 

#########
# usage #
#########
USAGE="$0 <fasta-input-file> <dtf-input-file> \
	<dict-outputfile> <gene_info-outputfile> \
	<exon-output-file> <intron-output-file>
	<ambiguous-bed-output-file> <exon-unambigous-bed-output-file> \
	<intron-unambigous-bed-output-file>"

# Check for proper usage.
[ $# -eq 5 ] || { printf '%s\n' "usage: $USAGE" >&2; exit 1; } 

#################
# parameters ##
################
FASTA=$1
GTF=$2

DICT=$3
GENE_INFO=$4 ##gene_info

EXONS=$5
INTRONS=$6
AMBIGUOUS=$7
EXONS_UNAMBIGUOUS=$8
INTRONS_UNAMBIGUOUS=$9

THREADS=8
MEMORY=10G

############################
## tools & commonds
###########################
DROPSEQ_TOOLS_DIR="Drop-seq_tools-2.0.0"
BEDTOOLS="bedtools"
BEDTOOLS_INTERSECT="$BEDTOOLS intersect -sorted -s"
SORT="sort --parallel=$THREADS -S$MEMORY"
SORT_BED="$SORT -k1,1 -k2,3n"
PIGZ="pigz -p $THREADS"
COMPRESS="$PIGZ --best"
BEDTOOLS_SUBTRACT="$BEDTOOLS subtract -sorted -s"

############################
## generate DropseqTools annotation
###########################
prefix=`basename "$GTF" .sorted.gtf`
echo $prefix

### generate refFlat file for Dropseq tools from reference genome
java -jar $DROPSEQ_TOOLS_DIR//3rdParty/picard/picard.jar CreateSequenceDictionary \
	REFERENCE=$FASTA\
	OUTPUT=$DICT

### ConvertToRefFlat
$DROPSEQ_TOOLS_DIR/ConvertToRefFlat 
	ANNOTATIONS_FILE=$GTF \
	SEQUENCE_DICTIONARY=$DICT \
	OUTPUT="${prefix}.refFlat"

#######################################
### retrive gene names from annotation
######################################
## get all the genes in the filtered gtf file
cat $GTF \
	| awk 'BEGIN{FS="\t"}{split($9,a,";"); if($3~"gene") print a[1]"\t"a[3]"\t"$1"\t"$4"\t"$5"\t"a[2]"\t"$7"\t"$5-$4}' \
	| sed 's/gene_id "//' \
	| sed 's/gene_id "//' \
        | sed 's/gene_type "//' \
	| sed 's/gene_name "//' \
	| sed 's/"//g' \
	| sed "1i\Geneid\tGeneSymbol\tChromosome\tStart\tEnd\tClass\tStrand\tLength"  \
	| cut -f2,3,4,5 | sed 's/^ *//g' > $GENE_INFO


###################################################
## Convert gtf to unique exon and intron bed files
##################################################
## exons
cat $GTF \
    | awk 'OFS="\t" { if ($3=="exon") print $1, $4-1, $5, substr($20, 2, length($20) - 3), ".", $7}' \
    | sort --parallel=8 -S10G -k1,1 -k2,2n -k3,3n \
    | awk -F\\t -vOFS=\\t '$4="."' \
    | uniq \
    | pigz --best \
    > $EXONS
## introns
cat $GTF \
	| awk 'OFS="\t" { if ($3=="exon") print $1, $4-1, $5, substr($20, 2, length($20) - 3), ".", $7}' \
	| sort --parallel=8 -S10G -k1,1 -k6,6 -k4,4 -k2,2n -k3,3n \
	| awk -F\\t -vOFS=\\t '$4!=prev{cnt=0}cnt{print $1, end, $2, prev ";" cnt, ".", $6}{++cnt;prev=$4;end=$3}' \
	| sort --parallel=8 -S10G -k1,1 -k2,2n -k3,3n \
	| awk -F\\t -vOFS=\\t '$4="."' \
	| uniq \
	| grep -v ^chrM \
	| pigz --best \
	> $INTRONS

#########################################################
## clean bed files by annotating ambiguous regoins
########################################################

####################
# ambigous regions #
####################

# Identify ambiguous regions as overlaps between introns and exons.
$BEDTOOLS_INTERSECT \
  -a "$INTRONS" -b "$EXONS" \
| $SORT_BED \
| $COMPRESS \
> $AMBIGUOUS

###########################
# unambigous exon regions #
###########################

# Identify unambiguous exon regions as difference between exons and
# ambigous regions.
$BEDTOOLS_SUBTRACT \
  -a "$EXONS" -b "$AMBIGUOUS" \
| $SORT_BED \
| $COMPRESS \
> $EXONS_UNAMBIGUOUS

#############################
# unambigous intron regions #
#############################

# Identify unambiguous intron regions as difference between introns and
# ambigous regions.
$BEDTOOLS_SUBTRACT \
  -a "$INTRONS" -b "$AMBIGUOUS" \
| $SORT_BED \
| $COMPRESS \
> $INTRONS_UNAMBIGUOUS
