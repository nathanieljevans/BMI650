# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 09:10:45 2018

@author: nathaniel evans 
@class: ALGORITHMS 
@HW: 1 
"""

'''
3. Create the pseudocode (70 pts) that does the following:

Reads in a FASTA file of nucleotide (DNA) sequences 
Carry out the translation to a protein sequence in all possible frames 
Finds  putative cleavage sites which are “double basic" (e.g., KK, KR, RK, RR) in the protein sequence 

Outputs the resulting protein sequences and the possible locations in a FASTA file that carries the description field from the nucleotide entry. The locations of the cleavage sites should be described in the description field 
'''
# open file (fasta)
# read file and split into dictionary where key (gene) -> list(values)
# for each nucleotide sequence: 

    # def generator function that takes gene seq and returns the next codon in the list
    # init 3 generator functions starting at each reading frame 
    # create empty dictionary of lists
    # while the generators aren't empty, loop through and convert each codon into its amino acid - How should we deal with stop codons? 
    
    #for each reading frame list of amino acids, loop through and compare previous amino acid with current aa searching for matches within {KK, KR, RK, RR}, 
    # when a double basic is found, record it in a list, including index info (from beginning of aa seq) and reading frame. 
    # open new fasta file for output
    # convert aa lists to strings
    # print aa strings, with labeled reading frames to fasta file 
    # print double basic clevage site locations to description field

# FASTA format goal (something like) : 
''' 
> gene_pa1_rdfm01 (name) | description (from ns) | double basic cleavage sites | amino acid string
> ... 

''' 




'''
4. Create a python script that then implements each of the above tasks (point breakdown below -Note: you can still pass this assignment if your programming is rough given point breakdown -emphasis is  on ordered, concrete tasks !!):

Reads in a FASTA file of nucleotide (DNA) sequences (5 pts)
Carry out the translation to a protein sequence in all possible frames (10 pts) 
Finds  putative cleavage sites which are “double basic" (e.g., KK, KR, RK, RR) in the protein sequence (5 pts)

Outputs the resulting protein sequences and the possible locations in a FASTA file that carries the description field from the nucleotide entry. The locations of the cleavage sites should be described in the description field (5 pts) 

Please test your code on the supplied FASTA file “pa1.fasta” and supply the out put file along with your python code. Note: it should run on any FASTA file :) 
'''

# relative 
FASTA_PATH = "pa1.fasta" 

class gene: 
    
    def __init__(self, name, seq, desc):
        self.geneName = name
        self.seq = seq
        self.desc = desc 
    
    def translate(): 
        cod_gen = __codon_gen(self.seq)
        self.aa_rd0 = [] 
        self.aa_rd1 = [] 
        self.aa_rd2 = [] 
        
        for codons in gen: 
            aa_rd0.append(__codon2aa(codons[0]))
            aa_rd1.append(__codon2aa(codons[1]))
            aa_rd2.append(__codon2aa(codons[2]))
            
    def __find_double_basics(aas):
        clev_loc = []
        double_basic = {"KK", "KR", "RK", "RR"}
        last = aas[0]
        for i,aa in enumerate(aas[1:]): 
            twofer = last + aa
            last = aa 
            if twofer in double_basic: 
                clev_loc.append(i)
                
        return clev_loc 
    
    def get_cleavage_locations(): 
        self.clv_loca_rd0 = __find_double_basics(self.aa_rd0)
        self.clv_loca_rd1 = __find_double_basics(self.aa_rd1)
        self.clv_loca_rd2 = __find_double_basics(self.aa_rd2)
        
    def get_fasta_string(): 
        faa_st0 = '>' + self.geneName + '|reading frame 0|'
                        + self.desc + '|' 
                        + ''.join(self.aa_rd0) + '|' 
                        + ''.join(self.clev_loca_rd0)
                        
            
        
            

    def __codon2aa(codon): 
        
        # taken from https://stackoverflow.com/questions/19521905/translation-dna-to-protein
        codontable = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',}
        
        return codontable[codon]
        
            

    def __codon_gen(seq): 
        i = 0
        while (i < len(seq)-2)
            rd0 = seq[i:3]
            rd1 = seq[i+1:4]
            rd2 = seq[i+2:5]
            yield (rd0, rd1, rd2)
   
if __name__ == __main__: 
    
    print('finished')

