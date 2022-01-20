#!/data/murphy/home/haliu/anaconda3/bin/python

# This is a modified script to tag all the possible conversion and their positions on the reads
# v0.003
import os
import pysam
import numpy as np
import pandas as pd
import sys
inbam = sys.argv[1]
outbam = sys.argv[2]

bamfile =pysam.AlignmentFile(inbam, 'rb')
mod_bamfile = pysam.AlignmentFile(outbam, mode='wb',template=bamfile)
qual = 20

def createTag(d):
    return ''.join([''.join(key) + str(d[key]) + ';' for key in d.keys()])[:-1]

for read in bamfile.fetch(until_eof=True):
    if read.has_tag('MD') == True:
        specific_conversions = {}
        total_content = {'a' : 0, 'c' : 0, 'g' : 0, 't' : 0}
        specific_conversions[('c', 'A')] = 0
        specific_conversions[('g', 'A')] = 0
        specific_conversions[('t', 'A')] = 0
        specific_conversions[('a', 'C')] = 0
        specific_conversions[('g', 'C')] = 0
        specific_conversions[('t', 'C')] = 0
        specific_conversions[('a', 'G')] = 0
        specific_conversions[('c', 'G')] = 0
        specific_conversions[('t', 'G')] = 0
        specific_conversions[('a', 'T')] = 0
        specific_conversions[('c', 'T')] = 0
        specific_conversions[('g', 'T')] = 0
        specific_conversions[('a', 'N')] = 0
        specific_conversions[('c', 'N')] = 0
        specific_conversions[('g', 'N')] = 0
        specific_conversions[('t', 'N')] = 0
        ## add deletion events
        specific_conversions[('A', 'N')] = 0
        specific_conversions[('C', 'N')] = 0
        specific_conversions[('G', 'N')] = 0
        specific_conversions[('T', 'N')] = 0

        accept_ref_names= ['chrX','chrY']
        for i in range(1,23):
            accept_ref_names.append('chr'+str(i))

        tC_loc = []
        aG_loc = []

        tC_read_loc = []
        aG_read_loc = []

        try:
            refseq = read.get_reference_sequence().lower()
        except (UnicodeDecodeError):
            refseq=''

        for base in total_content.keys():
            total_content[base] += refseq.count(base)
        #print(total_content)

        for pair in read.get_aligned_pairs(with_seq=True):
            try:
                ## deletions
                if pair[0] is None and pair[1] is not None and pair[2] is not None:
                    specific_conversions[(pair[2],'N')] += 1
                ## not indels
                if pair[0] is not None and pair[1] is not None and pair[2] is not None:
                    if str(pair[2]).islower() and not read.query_qualities[pair[0]] < qual:
                        specific_conversions[(pair[2],read.seq[pair[0]])] += 1
                        if (pair[2],read.seq[pair[0]]) == ('t', 'C'):
                            tC_loc.append(pair[1])
                            tC_read_loc.append(pair[0])
                        if (pair[2],read.seq[pair[0]]) == ('a', 'G'):
                            aG_loc.append(pair[1])
                            aG_read_loc.append(pair[0])
            except (UnicodeDecodeError, ValueError, KeyError):
                continue
        SC_tag = createTag(specific_conversions)
        TC_tag = createTag(total_content)
        #print(tC_loc)
        if len(tC_loc) == 0:
            tC_loc.append(0)
        if len(aG_loc) == 0:
            aG_loc.append(0)
        if len(tC_read_loc) == 0:
            tC_read_loc.append(0)
        if len(aG_read_loc) == 0:
            aG_read_loc.append(0)

        #tags = SC_tag, TC_tag, tC_loc, aG_loc
        #print(tags)
        read.set_tag('SC',SC_tag,'Z')
        read.set_tag('TC',TC_tag,'Z')
        read.set_tag('TL',tC_loc)#,'B')
        read.set_tag('AL',aG_loc)#,'B')
        #tag: CL: TC position on reads
        #tag: GL: AG position on reads
        read.set_tag('tC', tC_read_loc)
        read.set_tag('aG', aG_read_loc)

        mod_bamfile.write(read)
    else: 
        mod_bamfile.write(read)
bamfile.close()
mod_bamfile.close()
print('Wrote tags to {}'.format(outbam))
