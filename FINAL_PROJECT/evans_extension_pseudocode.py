#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 16:17:16 2018

@author: nathaniel j evans
"""

# Read in RNA-seq data
RNA_seq = open_and_parse('./data/rna_seq.fastq') # array of 100 bp nt seqa 

RNA_seq_subset = RNA_seq.select_subset(percentage=0.01, seed=0) # select 0.01 % of the subsets randomly  

#load subset of genome 
chr1_subset = open_and_parse('./data/chr1', chunk=100000) # select a chunk of the genome to align to

# write to fasta files 
for rna in RNA_seq_subset: 
    with open ('rna%s.fasta' %rna, 'w') as f: 
        f.write(header)
        f.write(chr1_subset) 
        f.write(rna) 

for path in os.listdir('*.fasta'): 
    mypackage.runClustalW(path)
 
variation = []
for path in os.listdir('*.aln'): 
    alignment = open_and_read(path)
    alignment = check_global_optimal_alignment() # there will be lots of data that doesn't align at all, these must be discarded. 
    variation.append(alignment.get_edit_distance() ) 
    
    
T = np.mean(variation) + 4*np.std(variation)

print('Suggested T based on given RNA-seq data: %d' %T)

