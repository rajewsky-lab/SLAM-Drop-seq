#!/data/murphy/home/haliu/anaconda3/bin/python

import os
import sys
import fileinput
import pysam

## input
inbam = sys.argv[1]
inCellList = sys.argv[2]
outbam = sys.argv[3]


bamfile =pysam.AlignmentFile(inbam, 'rb')
mod_bamfile = pysam.AlignmentFile(outbam, mode='wb',template=bamfile)

# Best cell list
cell_list = list() 
with open(inCellList) as f: 
    for line in f: 
        cell_list. append(line.strip())
#print(cell_list)
print("The number of cells selected is: ", len(cell_list))

## get XF, gf, gn, gs tags
for read in bamfile.fetch(until_eof=True):
    cell = read.get_tag("XC")
    if cell in cell_list:
        ## check the flag if a read is unmapped or not
        if read.is_unmapped:
            ## set an TS tag to unmmaped
            type = "unmapped"
        else:
            ## check the mapq score
            if read.mapping_quality != 255:
                type = "multimapped"
            else:
                ## for uniqly mapping reads, first check the strand that reads mapped to
                if read.is_reverse:
                    strand = '-'
                else:
                    strand = '+'

                # Since we noticed that we could get different function assignment(XF) for the same mapping,
                # we igore the XF tag and assign these reads mapped to more than one features as ambiguous
                if read.has_tag('gf') == True:
                    #XF = read.get_tag('XF')
                    gf_list = read.get_tag('gf').split(',')
                    gn_list = read.get_tag('gn').split(',')
                    gs_list = read.get_tag('gs').split(',')

                    # check if the read mapped to some genes on the same strand
                    if strand in gs_list:
                        try:
                            ## match strand with gs tag
                            gene_index = [i for i,val in enumerate(gs_list) if val==strand]
                            gs_list = [gs_list[i] for i in gene_index]
                            gn_list = [gn_list[i] for i in gene_index]
                            gf_list = [gf_list[i] for i in gene_index]

                            # get the length of unique gene names
                            gene_number = len(set(gn_list))

                            if gene_number == 1:
                                type = "assigned"
                                gene = gn_list[0]
                                ## Only if the read is assigned, create Gene tag: GN
                                read.set_tag('GN',gene,'Z')
                            else:
                                type = "ambiguous"

                        except ValueError:
                            continue
                    else:
                        type = "unassigned"

                else:
                    type = "unassigned"

        ## Create assignment tag: TS
        read.set_tag('TS', type, 'Z')
        mod_bamfile.write(read)

## save bamfile    
bamfile.close()
mod_bamfile.close()

## index out bamfile
pysam.samtools.index(outbam)




