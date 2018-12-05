#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 16:17:16 2018

@author: nathaniel evans
@class: BMI650 
@HW: FINAL PROJECT
"""

# Read in RNA-seq data
RNA_seq = open_and_parse('./data/rna_seq.fastq') # array of 100 bp nt seqa 

RNA_seq_subset = RNA_seq.select_subset(percentage=0.01, seed=0) # select 0.01 % of the subsets randomly  

# load genome into memory
chr1_subset = open_and_parse('./data/chr1')

# write to fasta files 
for rna in RNA_seq_subset: # prepare data for clustalW alignment 
    with open ('rna%s.fasta' %rna, 'w') as f: 
        f.write(header)
        f.write(chr1_subset) 
        f.write(rna) 

for path in os.listdir('*.fasta'): # run an optimal global pairwise alignment between B6 genome and PWK RNA-seq read  
    mypackage.runClustalW(path)
 
variation = []
for path in os.listdir('*.aln'): # calculate the variation in our alignments 
    alignment = open_and_read(path) 
    variation.append(alignment.get_edit_distance() ) 
    
T = np.mean(variation) + 4*np.std(variation) # Set threshold ...or calculate the T from sensitivity calculations as done in my main algorithm. 

print('Suggested T based on given RNA-seq data: %d' %T)

