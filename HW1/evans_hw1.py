# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 09:10:45 2018

@author: nathaniel evans 
@class: BMI650 ALGORITHMS 
@HW: 1 
"""

'''
PSEUDOCODE 
'''

# open file (fasta) or pipe in fasta file from bash command line 
# read file and save into individual gene info to data structure
# for each nucleotide sequence: 

    # transcribe the DNA seq: flip the kbs of seq 
    # create 3 reading frames with initial base at, i0, i0+1, i0+2, such that i0 = kb[0]
    # for each reading frame
        # increment through seq, selecting appropriate codons for each reading frame 
            # at each codon, convert to amino acid and store in appropriate data structure maintaining aa order 
    
    # for each sequence of amino acids (aa)
        # loop through aa sequence and compare previous aa with current aa searching for matches within {KK, KR, RK, RR}, 
        # when a double basic is found, record it, including index info (from beginning of aa seq) 
    # open new fasta file for output
    # convert aa lists to appropriate text format
    # print aa strings, with labeled reading frames to fasta file 
    # print double basic clevage site locations to description field

# FASTA format goal (something like) : 
''' 
> gene_pa1 (name) | reading frame 0X | description (from ns) | double basic cleavage sites | amino acid string rd 0 
> ... 

''' 

'''
PYTHON CODE 

To run this file, pipe a fasta file into this script from the command line using the format (bash): 
        $cat your-file.fasta | python evans_hw1.py 
    
    output will be saved in the same directory as this script, named: "evans_output.fasta" 
    
    Additional notes: 
        Docstring does include doctesting which can be run by specifying additional argument from command line (bash): 
            $python evans_hw.py --doctest
'''

import sys

class gene: 
    '''
   This class represents a gene sequence and has methods to emulate concepts of trascription, translation and observation of double basic cleavage location. 
   
    >>> gene("my_gene_name", "ATCGCTGCATGCATCGA", "my description | here | and here|").seq
    'ATCGCTGCATGCATCGA'
    >>> g = gene("my_gene_name", "ATCGCTGCATGCATCGA", "my description | here | and here|")
    >>> g.transcribe()
    >>> g.mRNA
    'UAGCGACGUACGUAGCU'
    >>> g = gene("my_gene_name", "ATCGCTGCATGCATCGA", "my description | here | and here|")
    >>> g.transcribe()
    >>> g.translate()
    >>> [g.aa_rd0, g.aa_rd1, g.aa_rd2]
    [['_', 'R', 'R', 'T', '_'], ['S', 'D', 'V', 'R', 'S'], ['A', 'T', 'Y', 'V', 'A']]

    Note:
        During translation, stop and start codons are not taken into account, rather, the translated protein is contiguous from beginning of mRNA to end. 

    Attributes:
        geneName (str): name of gene
        desc (str): description line of fasta file, exempt of gene name
        seq (str): nucleotide sequence represented by string of A,T,C,G

    '''
    
    def __init__(self, name, seq, desc):
        '''
            init method of gene class
        Args:
            name (str) : seq of name
            seq (str) : sequence of kbs
            desc (str) : gene description 
        Returns:
           None
        '''
        self.geneName = name
        self.seq = seq
        self.desc = desc 
    
    def transcribe(self): 
        '''
            Creates new attribute (mRNA) with the complementary bases of nucleotide seq with T -> U
        Args:
            None
        Returns:
           None
        '''
        self.mRNA = ''
        rna_comp = {'A':'U', 'T':'A', 'C':'G', 'G':'C'}
        for kb in self.seq: 
            self.mRNA += rna_comp[kb]
    
    def translate(self): 
        '''
            creates 3 new attributes representing aa sequence coded for by mRNA at each possible reading frame
        Args:
            None
        Returns:
           None
        '''
        cod_gen = self.__codon_gen(self.mRNA)
        self.aa_rd0 = [] 
        self.aa_rd1 = [] 
        self.aa_rd2 = [] 
        
        for codons in cod_gen: 
            if (len(codons[0]) == 3): 
                self.aa_rd0.append(self.__codon2aa(codons[0]))
            if (len(codons[1]) == 3): 
                self.aa_rd1.append(self.__codon2aa(codons[1]))
            if (len(codons[2]) == 3): 
                self.aa_rd2.append(self.__codon2aa(codons[2]))
            
    def __find_double_basics(self,aas):
        '''
            finds the location of specific sub strings within the aa sequence representing double base cleavage points and creates a list of loci. 
        Args:
            aas (str) 
        Returns:
           clev_loc (list) : list representing all double base cleavage indices 
        '''
        clev_loc = []
        double_basic = {"KK", "KR", "RK", "RR"}
        last = aas[0]
        for i,aa in enumerate(aas[1:]): 
            twofer = last + aa
            last = aa 
            if twofer in double_basic: 
                clev_loc.append(i)
                
        return clev_loc 
    
    def find_cleavage_locations(self): 
        '''
            sets new attributes representing each reading frames double base cleavage locations by index from start of amino acid sequence 
        Args:
            None
        Returns:
           None
        '''
        self.clv_loca_rd0 = self.__find_double_basics(self.aa_rd0)
        self.clv_loca_rd1 = self.__find_double_basics(self.aa_rd1)
        self.clv_loca_rd2 = self.__find_double_basics(self.aa_rd2)
        
    def get_fasta_string(self): 
        '''
            generates an appropriate fasta formated text representing one gene
        Args:
            None
        Returns:
           (str) fasta gene 
        '''              
        
        faa_st0 = '>' + self.geneName + '|reading frame 0|' + self.desc + '|' + ','.join(str(x) for x in self.clv_loca_rd0) + '\n' + ''.join(str(x)+'\n' if (i%80==0 and i>0) else x for (i,x) in enumerate(self.aa_rd0))
        
        faa_st1 = '>' + self.geneName + '|reading frame 1|' + self.desc + '|' + ','.join(str(x) for x in self.clv_loca_rd1) + '\n' + ''.join(str(x)+'\n' if (i%80==0 and i>0) else x for (i,x) in enumerate(self.aa_rd1))
                        
        faa_st2 = '>' + self.geneName + '|reading frame 2|' + self.desc + '|' + ','.join(str(x) for x in self.clv_loca_rd2) + '\n' + ''.join(str(x)+'\n' if (i%80==0 and i>0) else x for (i,x) in enumerate(self.aa_rd2))
                        
        return faa_st0 + '\n' + faa_st1 + '\n' + faa_st2 + '\n'
                        
    def __codon2aa(self,codon): 
        '''
            converts given codon (mRNA) to amino acid
        Args:
            codon (str) : three characters representing a codon
         Returns:
           (str) : single character representing amino acid
        '''
        assert len(codon) == 3, 'improper codon length'
        # taken from https://stackoverflow.com/questions/19521905/translation-dna-to-protein
        # but converted T -> U for mRNA
        codontable = {
        'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M',
        'ACA':'U', 'ACC':'T', 'ACG':'T', 'ACU':'T',
        'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',
        'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
        'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
        'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
        'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
        'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
        'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
        'UAC':'Y', 'UAU':'Y', 'UAA':'_', 'UAG':'_',
        'UGC':'C', 'UGU':'C', 'UGA':'_', 'UGG':'W',}
        
        return codontable[codon]
        
    def __codon_gen(self,seq):
        '''
            This is a generator function that returns the next codon in a given mRNA seq at each possible reading frame
        Args:
            seq (str) : mRNA sequence
        Returns:
           (tuple<str>) : 3-len tuple where each value is a three character string representing the next codon
        '''
        i = 0
        while (i < len(seq)-2): 
            rd0 = seq[i:i+3]
            rd1 = seq[i+1:i+4]
            rd2 = seq[i+2:i+5]
            i += 3
            yield (rd0, rd1, rd2)
   
# this handles file input, highlevel processing order and file output
if (__name__ == '__main__'): 
    # doctesting 
    if (len(sys.argv) > 1): 
        if (sys.argv[1] == "--doctest"): 
            import doctest 
            doctest.testmod()
            print('doctest complete')
        else: 
            print("did you mean to use --doctest ?")
    
    # running on a .fasta file with kb sequences 
    else:
        gene_store = [] 
        fasta_data = sys.stdin.read() 
        
        genes = fasta_data.split('>')
        failed_gene_parse = 0 
        
        for _gene in genes: 
            if (len(_gene) > 0) : 
                
                try: 
                    pcs = _gene.splitlines()
                    seq = ''.join(pcs[1:]) #this is super messy -ick 
                    pcs2 = pcs[0].split('|')
                    geneName = pcs2[0]
                    desc = '|'.join(pcs2[1:])
                    assert len(seq) > 0 , 'nucleotide sequence must be nonzero string'
                    gene_store.append(gene(name=geneName, seq=seq, desc=desc))
                except: 
                    failed_gene_parse += 1
             
        with open('evans_output.fasta', 'w') as f: 
            for i,g in enumerate(gene_store): 
                print('gene # ', i)
                print('transcribing...')     
                g.transcribe()    
                print('translating...')
                g.translate()
                print('analyzing for clevage locations...')
                g.find_cleavage_locations()
                print('writing output to fasta file...')
                f.write(g.get_fasta_string()) 
    
        print('Gene analysis finished. Number of genes that failed to parse: ', failed_gene_parse)
        print('Output file written to local directory as \"evans_output.fasta\"')

