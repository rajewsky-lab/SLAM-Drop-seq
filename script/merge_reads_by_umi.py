#!/data/murphy/home/haliu/anaconda3/bin/python
# coding: utf-8
import pysam
import os
import sys
import pandas as pd
import numpy as np
import string
import multiprocessing as mp
import time
import re

def merge_reads_with_coverage(inbam, gene_input, chrom_input, gene_start_input, gene_end_input):
    bamfile = pysam.AlignmentFile(inbam, 'rb')
    ## initiate some dict/list to store values
    ## keys are cb_umi_gene
    d_ref_pos = {}    ## overlapping reference positon
    d_ref_seq = {}    ## overlapping reference sequence
    d_tc_pos = {}     ## tC conversion position on the reference
    d_ag_pos = {}     ## aG conversion position on the reference
    d_gene = {}       ## annotated genes
    d_ss = {}         ## annotated splice status

    d_read_count = {}   ## reads number for each UMI
    d_strand = {}       ## overlapping strand

    ## record the unique cb_umi_gene/transcript number
    cb_umi_gene_list = []  ##list
    umi_count = 0
    read_count = 0

    chrom = chrom_input
    gene_start = gene_start_input
    gene_end = gene_end_input
    
    for read in bamfile.fetch(chrom, gene_start, gene_end):
        if read.has_tag("GN"):
            gene = read.get_tag("GN")
            if gene == gene_input:
                read_count += 1
                ##############################
                ## CB_UMI_gene identity
                ##############################
                cb_umi_gene = read.get_tag('XC') + '_' + read.get_tag('XM') + "_" + read.get_tag('GN')

                #############################
                ## TC and AG locs && strand
                #############################
                ## take AG on '-' strand and 'TC' on'+' strand
                if read.is_reverse:
                    tc_locs = read.get_tag('AL').tolist()
                    ag_locs = read.get_tag('TL').tolist()
                    strand = '-'
                else:
                    tc_locs = read.get_tag('TL').tolist()
                    ag_locs = read.get_tag('AL').tolist()
                    strand = '+'

                ##################################
                ## Splice status of read
                ##################################
                ss = read.get_tag('SS')

                ##################################
                ## GENE ASSIGNED
                ##################################
                gene = read.get_tag('GN')


                ## if the cb_umi is seen for the first time, create a new key.
                if cb_umi_gene not in cb_umi_gene_list:
                    umi_count += 1
                    cb_umi_gene_list.append(cb_umi_gene)

                #############################################################
                ## reads with reference sequence and corresponding positions
                #############################################################
                    pos = []
                    ref = ""
                    ## get the read information: postions and sequence on the reference
                    for pair in read.get_aligned_pairs(with_seq=True):
                        #print(pair)
                        try:
                            if pair[2] is not None:
                                pos.append(pair[1])
                                ref += str(pair[2]).upper()
                            else:
                                continue

                        except (ValueError,KeyError):
                            print("Error")
                    #print(pos)
                    d_ref_pos[cb_umi_gene] = pos
                    d_ref_seq[cb_umi_gene] = ref
                    d_gene[cb_umi_gene] = [gene]
                    d_ss[cb_umi_gene] = [ss]
                    d_read_count[cb_umi_gene] = [1]
                    d_strand[cb_umi_gene] = [strand]

                    if tc_locs[0] != 0:
                        d_tc_pos[cb_umi_gene] = tc_locs
                    else:
                        d_tc_pos[cb_umi_gene] = []

                    if ag_locs[0] != 0:
                        d_ag_pos[cb_umi_gene] = ag_locs
                    else:
                        d_ag_pos[cb_umi_gene] = []

                ## if the cb_umi is areadty seen in the keys, add the value to the key.
                else:
                    pos = []
                    ref = ""
                    for pair in read.get_aligned_pairs(with_seq=True):

                        try:
                            if pair[2] is not None:
                                ## record overlapping positions
                                pos.append(pair[1])
                                if pair[1] not in d_ref_pos[cb_umi_gene]:
                                    ref += str(pair[2]).upper()  ## only record ref sequence from non-overlapping postion
                                else:
                                    continue
                            else:
                                continue
                        except (ValueError,KeyError):
                            print("Error")

                    d_ref_pos[cb_umi_gene].extend(pos)
                    d_ref_seq[cb_umi_gene] += ref
                    d_read_count[cb_umi_gene].extend([1])
                    d_strand[cb_umi_gene].extend([strand])

                    if gene not in d_gene[cb_umi_gene]:
                        d_gene[cb_umi_gene].extend([gene])

                    if ss not in d_ss[cb_umi_gene]:
                        d_ss[cb_umi_gene].extend([ss])

                    if tc_locs[0] != 0:
                        d_tc_pos[cb_umi_gene].extend(tc_locs)

                    if ag_locs[0] != 0:
                        d_ag_pos[cb_umi_gene].extend(ag_locs)

                    if read_count%1e6==0:
                        print("Processing read", read_count, "from UMI", umi_count)


    ####    ############################################################
    ## count #Ts and #As on the reference for each cb_umi_gene
    ################################################################
    d_t = {}
    d_a = {}
    for key, item in d_ref_seq.items():
        if d_strand[key][0] == '+':
            d_t[key] = item.count('T')
            d_a[key] = item.count('A')
        else:
            d_t[key] = item.count('A')
            d_a[key] = item.count('T')

    ###################################################
    ## count tC number and aG per cb_umi_gene
    ###################################################
    d_tc = {key:len(set(tc)) for key,tc in d_tc_pos.items()}
    d_ag = {key:len(set(ag)) for key,ag in d_ag_pos.items()}


    #########################################################################
    ## count tC & aG coverage on each converted position for each cb_umi_gene
    #########################################################################
    d_tc_cover = {}
    d_ag_cover = {}
    for key, item in d_tc_pos.items():
        #print(key, item)
        ttt=(pd.Index(item).value_counts())
        d_tc_cover[key] = ttt.to_dict()

    for key, item in d_ag_pos.items():
        ttt=(pd.Index(item).value_counts())
        d_ag_cover[key] = ttt.to_dict()

    ###########################################################################
    ## count reads coverage of each tC converted position for each cb_umi_gene
    ###########################################################################
    d_tc_ref_cover = {}
    for cb_umi_gene, element in d_tc_cover.items():
        d_tc_ref_cover[cb_umi_gene] = {}
        if len(element) != 0:
            for pos, count in element.items():
                d_tc_ref_cover[cb_umi_gene][pos] = d_ref_pos[cb_umi_gene].count(pos)
        else:
            continue

    ### tC coverage to dataframe using pandas series
    ind = []
    lst = []
    for key, item in d_tc_cover.items():
        for pos, value in item.items():
            ind.append(key + "_" + str(pos))
            lst.append(value)
    # create Pandas Series with define indexes
    df_tc = pd.Series(lst, index = ind)
    df_tc = pd.DataFrame(df_tc).reset_index()
    df_tc.columns = ['id', 'tc_cover']
    df_tc.head()

    ### ref coverage to dataframe using pandas series
    ind = []
    lst = []
    for key, item in d_tc_ref_cover.items():
        for pos, value in item.items():
            #print(key + "_" + str(pos), value)
            ind.append(key + "_" + str(pos))
            lst.append(value)
    # create Pandas Series with define indexes
    df_tc_ref = pd.Series(lst, index = ind)
    df_tc_ref = pd.DataFrame(df_tc_ref).reset_index()
    df_tc_ref.columns = ['id', 'ref_cover']
    df_tc_ref.head()

    ## combine tc coverage and ref converage
    data_tc_cover=pd.concat([df_tc, df_tc_ref["ref_cover"]], axis=1, ignore_index=False)

    #data_tc_cover.to_csv(tc_outcsv, index=False)
    #print(data_tc_cover)

    ###############################################################################
    ## count reads coverage of each aG converted position for each cb_umi_gene
    ###############################################################################
    d_ag_ref_cover = {}
    for cb_umi_gene, element in d_ag_cover.items():
        d_ag_ref_cover[cb_umi_gene] = {}
        if len(element) != 0:
            for pos, count in element.items():
                d_ag_ref_cover[cb_umi_gene][pos] = d_ref_pos[cb_umi_gene].count(pos)
        else:
            continue



    ### aG coverage to dataframe using pandas series
    ind = []
    lst = []
    for key, item in d_ag_cover.items():
        for pos, value in item.items():
            ind.append(key + "_" + str(pos))
            lst.append(value)
    # create Pandas Series with define indexes
    df_ag = pd.Series(lst, index = ind)
    df_ag = pd.DataFrame(df_ag).reset_index()
    df_ag.columns = ['id', 'ag_cover']
    df_ag.head()


    ### ref coverage to dataframe using pandas series
    ind = []
    lst = []
    for key, item in d_ag_ref_cover.items():
        for pos, value in item.items():
            ind.append(key + "_" + str(pos))
            lst.append(value)
    # create Pandas Series with define indexes
    df_ag_ref = pd.Series(lst, index = ind)
    df_ag_ref = pd.DataFrame(df_ag_ref).reset_index()
    df_ag_ref.columns = ['id', 'ref_cover']
    df_ag_ref.head()

    ## combine ag coverage and ref converage
    data_ag_cover=pd.concat([df_ag, df_ag_ref["ref_cover"]], axis=1, ignore_index=False)
    #data_ag_cover.to_csv(ag_outcsv, index=False)

    ###############################################
    ## Merging SS tag and output table
    ###############################################
    ## create a empty datafram
    Outputdf = pd.DataFrame(columns=['cb_umi_gene', 'reads_count', 'length', 'ss', 't', 'tC', 'a', 'aG'])

    print("Deside the splice status of each cb_umi for gene: %s" % gene_input)

    ## decide the ss for each cb_umi
    index = 0
    for key, ss in d_ss.items():
        index += 1
        if len(ss) == 1:
            if ss == 'ambiguous':
                new_ss = 'ambiguous'
            elif ss == 'intron_only' or ss == 'intron_exon' or ss == 'intron_exon_exon':
                new_ss = 'unspliced'
            else:
                new_ss = 'spliced'

        else:
            if 'intron_exon' in ss or 'intron_only' in ss or 'intron_exon_exon' in ss:
                new_ss = 'unspliced'
            else:
                new_ss = 'spliced'

        Outputdf.loc[index] = [key, len(d_read_count[key]), len(d_ref_seq[key]), str(new_ss), d_t[key], d_tc[key], d_a[key], d_ag[key]]

    Outputdf.reset_index(drop=True)
    #Outputdf.to_csv(merged_outcsv, index=False)
    return(data_tc_cover, data_ag_cover, Outputdf)


if __name__ == "__main__":

    inbam = sys.argv[1]
    inGeneList = sys.argv[2]
    cores_to_use = sys.argv[3]
    merged_table = sys.argv[4]
    tc_table = sys.argv[5]
    ag_table = sys.argv[6]

    # Gene list
    #gene_list = list()
    #with open(inGeneList) as f:
    #    for line in f:
    #        gene_list.append(line.strip())

    #bamfile = pysam.AlignmentFile(inbam, "rb")
    filename = inbam.split(".")[0]
    print(filename)
    
    gene_list = list()
    chrom_list = list()
    gene_start_list = list()
    gene_end_list = list()
    with open(inGeneList) as f:
        for line in f:
            gene_list.append(re.split(r'\t+', line)[0])
            chrom_list.append(re.split(r'\t+', line)[1])
            gene_start_list.append(re.split(r'\t+', line)[2])
            gene_end_list.append(re.split(r'\t+', line)[3].strip())
    ## create files for each sample
    #if not os.path.isfile(os.path.join(path, filename + "_tc.csv")):
    #    tc_table = os.path.join(path, filename + "_tc.csv")
    #if not os.path.isfile(os.path.join(path, filename + "_ag.csv")):
    #    ag_table = os.path.join(path, filename + "_ag.csv")
    #if not os.path.isfile(os.path.join(path, filename + "_merged.csv")):
    #    merged_table = os.path.join(path, filename + "_merged.csv")
        
    ## return tables from merge
    #results = [pool.apply(howmany_within_range, args=(row, 4, 8)) for row in data]
    #results = [pool.apply(merge_reads_with_coverage, args=(inbam, gene)) for gene in gene_list]
    #print(results.head())
    #sub_tc_table, sub_ag_table, sub_merged_table = merge_reads_with_coverage(inbam, gene)
    
    ##########################
    print("multiprocessing...")
    t1 = time.time()
    #cores = int(mp.cpu_count())
    pool = mp.Pool(processes=int(cores_to_use))

    tasks = []
    for i in range(1, len(gene_list)):
        tasks.append((inbam, gene_list[i], chrom_list[i], int(float(gene_start_list[i])), int(float(gene_end_list[i]))))
    #result_tc, result_ag, result_merged = pool.starmap(merge_reads_with_coverage, tasks)

    results = pool.starmap(merge_reads_with_coverage, tasks)
    pool.close()
    
    #print(type(results))
    
    result_tc = []
    result_ag = []
    result_merged = []
    for i in range(1, len(gene_list)):
        print(gene_list[i])
        result_tc.append(results[i-1][0])
        result_ag.append(results[i-1][1])
        result_merged.append(results[i-1][2])
        print("\n")

    ## create files for each sample
    #if not os.path.isfile(os.path.join(path, filename + "_tc.csv")):
    #    tc_table = os.path.join(path, filename + "_tc.csv")
    #if not os.path.isfile(os.path.join(path, filename + "_ag.csv")):
    #    ag_table = os.path.join(path, filename + "_ag.csv")
    #if not os.path.isfile(os.path.join(path, filename + "_merged.csv")):
    #    merged_table = os.path.join(path, filename + "_merged.csv")

    pd.concat(result_tc).to_csv(tc_table, index=False)
    pd.concat(result_ag).to_csv(ag_table, index=False)
    pd.concat(result_merged).to_csv(merged_table, index=False)

    t2 = time.time()
    print("The total running time is:", (t2 - t1).__str__() + "\n")

