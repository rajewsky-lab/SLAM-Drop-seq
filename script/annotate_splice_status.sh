#!/bin/bash

# BAM to database conversion script
#
# Copyright (C) 2019  Marcel Schilling
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


#######################
# general information #
#######################

# file:        add_splice_status_tag.sh
# created:     2019-07-10
# last update: 2019-07-11
# last_update: 2019-08-02
# author:      Marcel Schilling <marcel.schilling@mdc-berlin.de>
# purpose:     add splice status tag to BAM file based on (sorted,
#              non-overlapping) exon/intron/ambiguous regions


######################################
# change log (reverse chronological) #
######################################
# 2019-08-02: Haiyue:split intron_exon reads to intron_exon and exon_ambiguous
#             if unspliced reads that > 50% overlap with intron, we call it intron_exon, otherwise exon-ambiguous.
#             The same for inron_exon_exon reads.
# 2019-07-11: added sorting of input BED files to match BAM file order
#             added exon-intergenic category
# 2019-07-10: initial version


#########
# usage #
#########

# Define usage.
USAGE="$0 <bam-file> <exon-unambiguous-bed> <intron-unambiguous-bed> \
<ambiguous-bed>"

# Check for proper usage.
[ $# -eq 4 ] || { printf '%s\n' "usage: $USAGE" >&2; exit 1; } 


##############
# parameters #
##############

ALIGNMENTS=$1

EXONS_UNSORTED=$2

INTRONS_UNSORTED=$3

AMBIGUOUS_UNSORTED=$4

TAG=SS:Z:

THREADS=8

TMP_DIR=$(mktemp \
            --directory \
            .tmp_$(basename --suffix='.sh' "$0").XXXXX)

GENOME="$TMP_DIR/genome.tsv"

EXONS="$TMP_DIR/exons.bed.gz"

INTRONS="$TMP_DIR/introns.bed.gz"

AMBIGUOUS="$TMP_DIR/ambiguous.bed.gz"

REMAINING="$TMP_DIR/remaining.bam"

SPLICED="$TMP_DIR/spliced.bam"

UNSPLICED="$TMP_DIR/unspliced.bam"


############
# commands #
############

SAMTOOLS="samtools"

SAMTOOLS_VIEW="$SAMTOOLS view -@ $THREADS"

BAM2HEADER="$SAMTOOLS_VIEW -H"

GAWK="awk"

BEDTOOLS="bedtools"

BEDTOOLS_SORT="$BEDTOOLS sort"

BAM2SAM="$SAMTOOLS_VIEW -h"

SAM2BAM="$SAMTOOLS_VIEW -b -"

BEDTOOLS_INTERSECT="$BEDTOOLS intersect -sorted -s -split"

BAM2SAM_NOHEADER="$SAMTOOLS_VIEW"

RM="rm"

REMOVE="$RM --recursive --"

SAMTOOLS_SORT="$SAMTOOLS sort -@ $THREADS"

SORT_BAM="$SAMTOOLS_SORT -"


#################
# bash settings #
#################

# Exit with error as soon as any command exists with error.
set -e

# Treat error in any part of pipeline as pipeline error.
set -o pipefail


###############
# BED sorting #
###############

# Generate genome file from SAM header.
$BAM2HEADER "$ALIGNMENTS" \
| $GAWK \
  '
    # Read/write TSV.
    BEGIN{
      FS="\t"
      OFS=FS
    }
   
    # Print the chromosome names and lengths.
    $1 == "@SQ" {
      print substr($2, 4), substr($3, 4)
    }
  ' \
> "$GENOME"

# Sort introns BED to match BAM order.
$BEDTOOLS_SORT \
  -g "$GENOME" \
  -i "$INTRONS_UNSORTED" \
> "$INTRONS"

# Sort exons BED to match BAM order.
$BEDTOOLS_SORT \
  -g "$GENOME" \
  -i "$EXONS_UNSORTED" \
> "$EXONS"

# Sort umbiguous regions BED to match BAM order.
$BEDTOOLS_SORT \
  -g "$GENOME" \
  -i "$AMBIGUOUS_UNSORTED" \
> "$AMBIGUOUS"


##############################
# unmapped / intergenic read #
##############################

{
  # Tag & print unmapped reads and keep reads assigned to genes for later.
  $BAM2SAM "$ALIGNMENTS" \
  | $GAWK \
      --assign=tag="$TAG" \
      '
        # Read/write TSV.
        BEGIN{
          FS = "\t"
          OFS = FS
  
        # Use STDERR to separate some reads for later.
          remaining = "/dev/stderr"
        }
  
        # Print header and copy for later.
        !header_done{
          if(/^@/){
              print
              print > remaining
              next
          }
          header_done = 1
        }
  
        # Keep reads assigned to genes for later.
        /	GN:Z:/{
          print > remaining
          next
        }
  
        # Tag & print unmapped & intergenic reads.
        {
          $(++NF) = tag (and($2, 4) ? "unmapped" : "intergenic")
        }
        1
    ' \
  2> >($SAM2BAM \
       > $REMAINING)
  
  
###################
# exon-only reads #
###################
  
  # Tag & print exon-only reads.
  $BEDTOOLS_INTERSECT \
    -g "$GENOME" \
    -wa -u \
    -f 1 \
    -a "$REMAINING" -b "$EXONS" \
  | $BAM2SAM_NOHEADER \
  | $GAWK \
      --assign=tag="$TAG" \
      '
        # Read/write TSV.
        BEGIN{
          FS = "\t"
          OFS = FS
        }
  
        # Tag & print exon-only reads.
        {
          $(++NF) = tag "exon_only"
        }
        1
    '
  
  
#####################
# intron-only reads #
#####################
  
  # Tag & print intron-only reads.
  $BEDTOOLS_INTERSECT \
    -g "$GENOME" \
    -wa -u \
    -f 1 \
    -a "$REMAINING" -b "$INTRONS" \
  | $BAM2SAM_NOHEADER \
  | $GAWK \
      --assign=tag="$TAG" \
      '
        # Read/write TSV.
        BEGIN{
          FS = "\t"
          OFS = FS
        }
  
        # Tag & print intron-only reads.
        {
          $(++NF) = tag "intron_only"
        }
        1
    '
  
  
########################
# ambiguous-only reads #
########################
  
  # Tag & print ambiguous-only reads.
  $BEDTOOLS_INTERSECT \
    -g "$GENOME" \
    -wa -u \
    -f 1 \
    -a "$REMAINING" -b "$AMBIGUOUS" \
  | $BAM2SAM_NOHEADER \
  | $GAWK \
      --assign=tag="$TAG" \
      '
        # Read/write TSV.
        BEGIN{
          FS = "\t"
          OFS = FS
        }
  
        # Tag & print ambiguous-only reads.
        {
          $(++NF) = tag "ambiguous"
        }
        1
    '
  
  
###################
# unspliced reads #
###################
  
  # Separate unspliced reads from spliced reads.
  $BEDTOOLS_INTERSECT \
    -g "$GENOME" \
    -wa -v \
    -f 1 \
    -a "$REMAINING" -b "$EXONS" "$INTRONS" "$AMBIGUOUS" \
  | $BAM2SAM \
  | $GAWK \
      '
        # Read/write TSV.
        BEGIN{
          FS = "\t"
          OFS = FS
  
        # Use STDERR to separate spliced reads.
          spliced = "/dev/stderr"
        }
  
        # Print header and copy for spliced reads.
        !header_done{
          if(/^@/){
              print
              print > spliced
              next
          }
          header_done = 1
        }
  
        # Separate spliced reads.
        $6 ~ /N/{
          print > spliced
          next
        }
  
        # Separate unspliced reads
        1
    ' \
  2> >($SAM2BAM \
       > $SPLICED) \
  | $SAM2BAM \
  > $UNSPLICED
  
  # Remove unseparated remaining reads file.
  $REMOVE "$REMAINING"
  
  
#############################
# ambiguous unspliced reads #
#############################
  
  # Tag & print ambiguous unspliced reads.
  $BEDTOOLS_INTERSECT \
    -g "$GENOME" \
    -wa -u \
    -a "$UNSPLICED" -b "$AMBIGUOUS" \
  | $BAM2SAM_NOHEADER \
  | $GAWK \
      --assign=tag="$TAG" \
      '
        # Read/write TSV.
        BEGIN{
          FS = "\t"
          OFS = FS
        }
  
        # Tag & print ambiguous reads.
        {
          $(++NF) = tag "ambiguous"
        }
        1
    '

  # Keep non-ambigous unspliced reads for later.
  $BEDTOOLS_INTERSECT \
    -g "$GENOME" \
    -wa -v \
    -a "$UNSPLICED" -b "$AMBIGUOUS" \
  > $REMAINING
  
  # Remove unseparated unspliced reads file.
  $REMOVE "$UNSPLICED"
  
  
#####################
# intron-exon reads #
#####################
  
  # Tag & print intron-exon reads.
  $BEDTOOLS_INTERSECT \
    -g "$GENOME" \
    -wa -u \
    -f 0.075 \
    -a "$REMAINING" -b "$INTRONS" \
  | $BAM2SAM_NOHEADER \
  | $GAWK \
      --assign=tag="$TAG" \
      '
        # Read/write TSV.
        BEGIN{
          FS = "\t"
          OFS = FS
        }
  
        # Tag & print intron-exon reads.
        {
          $(++NF) = tag "intron_exon"
        }
        1
    '

#########################
# exon-intergenic reads #
#########################
  
  # Tag & print exon-intergenic reads.
  $BEDTOOLS_INTERSECT \
    -g "$GENOME" \
    -wa -v \
    -a "$REMAINING" -b "$INTRONS" \
  | $BAM2SAM_NOHEADER \
  | $GAWK \
      --assign=tag="$TAG" \
      '
        # Read/write TSV.
        BEGIN{
          FS = "\t"
          OFS = FS
        }
  
        # Tag & print exon-intergenic reads.
        {
          $(++NF) = tag "exon_intergenic"
        }
        1
    '
  
  # Remove unseparated unspliced reads file.
  $REMOVE "$REMAINING"
  
  
###################
# exon-exon reads #
###################
  
  # Tag & print exon-exon reads.
  $BEDTOOLS_INTERSECT \
    -g "$GENOME" \
    -wa -v \
    -a "$SPLICED" -b "$INTRONS" \
  | $BAM2SAM_NOHEADER \
  | $GAWK \
      --assign=tag="$TAG" \
      '
        # Read/write TSV.
        BEGIN{
          FS = "\t"
          OFS = FS
        }
  
        # Tag & print exon-exon reads.
        {
          $(++NF) = tag "exon_exon"
        }
        1
    '
  
  # Keep non-exon-exon reads for later.
  $BEDTOOLS_INTERSECT \
    -g "$GENOME" \
    -wa -u \
    -a "$SPLICED" -b "$INTRONS" \
  > $REMAINING
  
  # Remove unseparated spliced reads file.
  $REMOVE "$SPLICED"
  
  
##########################
# intron-exon-exon reads #
##########################
  
  # Tag & print intron-exon-exon reads.
  $BEDTOOLS_INTERSECT \
    -g "$GENOME" \
    -wa -u \
    -f 0.075 \
    -a "$REMAINING" -b "$INTRONS" \
  | $BAM2SAM_NOHEADER \
  | $GAWK \
      --assign=tag="$TAG" \
      '
        # Read/write TSV.
        BEGIN{
          FS = "\t"
          OFS = FS
        }
  
        # Tag & print intron-exon-exon reads.
        {
          $(++NF) = tag "intron_exon_exon"
        }
        1
    '
  
###########################
# ambiguous spliced reads #
###########################
  
  # Tag & print ambiguous spliced reads.
  $BEDTOOLS_INTERSECT \
    -g "$GENOME" \
    -wa -v \
    -a "$REMAINING" -b "$EXONS" \
  | $BAM2SAM_NOHEADER \
  | $GAWK \
      --assign=tag="$TAG" \
      '
        # Read/write TSV.
        BEGIN{
          FS = "\t"
          OFS = FS
        }
  
        # Tag & print ambiguous spliced reads.
        {
          $(++NF) = tag "ambiguous"
        }
        1
    '

  # Remove remaining temporary data.
  $REMOVE "$TMP_DIR"

} \
| $SAM2BAM \
| $SORT_BAM
