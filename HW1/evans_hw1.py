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
# read file and save into individual gene info in to data structure
# for each nucleotide sequence: 

    # transcribe the DNA seq: flip the kbs of seq 
    # 
    # create 3 reading frames (i0, i0+1, i0+2) such that i0 = kb[0]
    # increment through seq, selecting appropriate codons for each reading frame 
        # at each codon, convert to amino acid and store in appropriate data structure maintaining aa order 
    
    # for each sequence of amino acids
        # loop through and compare previous aa with current aa searching for matches within {KK, KR, RK, RR}, 
        # when a double basic is found, record it, including index info (from beginning of aa seq) 
    # open new fasta file for output
    # convert aa lists to appropriate text format
    # print aa strings, with labeled reading frames to fasta file 
    # print double basic clevage site locations to description field

# FASTA format goal (something like) : 
''' 
> gene_pa1_rdfm01 (name) | description (from ns) | double basic cleavage sites | amino acid string rd 0 | "rd 1 | "rd 2
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
    
    def transcribe(self): 
    
        self.mRNA = ''
        rna_comp = {'A':'U', 'T':'A', 'C':'G', 'G':'C'}
        for kb in self.seq: 
            self.mRNA += rna_comp[kb]
    
    def translate(self): 
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
        self.clv_loca_rd0 = self.__find_double_basics(self.aa_rd0)
        self.clv_loca_rd1 = self.__find_double_basics(self.aa_rd1)
        self.clv_loca_rd2 = self.__find_double_basics(self.aa_rd2)
        
    def get_fasta_string(self): 
        #print(self.clv_loca_rd1)               
        
        faa_st0 = '>' + self.geneName + '|reading frame 0|' + self.desc + '|' + ','.join(str(x) for x in self.clv_loca_rd0) + '\n' + ''.join(str(x)+'\n' if (i%80==0 and i>0) else x for (i,x) in enumerate(self.aa_rd0))
        
        faa_st1 = '>' + self.geneName + '|reading frame 1|' + self.desc + '|' + ','.join(str(x) for x in self.clv_loca_rd1) + '\n' + ''.join(str(x)+'\n' if (i%80==0 and i>0) else x for (i,x) in enumerate(self.aa_rd1))
                        
        faa_st2 = '>' + self.geneName + '|reading frame 2|' + self.desc + '|' + ','.join(str(x) for x in self.clv_loca_rd2) + '\n' + ''.join(str(x)+'\n' if (i%80==0 and i>0) else x for (i,x) in enumerate(self.aa_rd2))
                        
        return faa_st0 + '\n' + faa_st1 + '\n' + faa_st2 + '\n'
                        
    def __codon2aa(self,codon): 
        
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
        i = 0
        while (i < len(seq)-2): 
            rd0 = seq[i:i+3]
            rd1 = seq[i+1:i+4]
            rd2 = seq[i+2:i+5]
            i += 3
            #print( (i, len(seq)) )
            #print( (rd0, rd1, rd2) )
            yield (rd0, rd1, rd2)
   
if (__name__ == '__main__'): 
    
    gene_store = [] 
    with open(FASTA_PATH, 'r') as f: 
        fasta_data = f.read() 
        genes = fasta_data.split('>')
        failed_gene_parse = 0 
        
        for _gene in genes: 
            if (len(_gene) > 0) : 
                
                try: 
                    pcs = _gene.splitlines()
                    #print(pcs)
                    seq = ''.join(pcs[1:]) #this is super messy -ick 
                    pcs2 = pcs[0].split('|')
                    geneName = pcs2[0]
                    desc = '|'.join(pcs2[1:])
                    #print(desc)
                    #print((geneName, desc, seq))
                    assert len(seq) > 0 , 'nucleotide sequence must be nonzero string'
                    gene_store.append(gene(name=geneName, seq=seq, desc=desc))
                except: 
                    failed_gene_parse += 1
         
    with open('evans_output.fasta', 'w') as f: 
        for i,g in enumerate(gene_store): 
            if ( i == 61): 
                print(g.seq)
            # print("analyzing gene %d of %d", % (i, len(gene_store))) # not sure why this wont work
            print('gene # ', i)
            print('transcribing...')     
            g.transcribe()    
            print('translating...')
            g.translate()
            print('analyzing for clevage locations...')
            g.find_cleavage_locations()
            print('writing output to fasta file...')
            f.write(g.get_fasta_string()) 


    
    print('finished. Number of genes that failed to parse: ', failed_gene_parse)

