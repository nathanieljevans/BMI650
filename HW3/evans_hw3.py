# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 19:58:59 2018

@author: nathaniel evans 
@class: BMI650 ALGS
@HW: 3 

Instructions: 
For this homework, please consider the attached alignment file (hw3.txt) which contains 3 sequences. 

If a column has a gap, for this assignment ignore the column.

Provide well commented code (python or otherwise) that can take this file and calculate two scores:

entropy  (50 pts)
sum of pairs  (50 pts)
 
Please utilize a BLOSUM62 substitution matrix for (ii)

"""

import re
import numpy as np


class substitution_matrix: 
    
    
    def __init__(self, file): 
        with open(file) as f: 
            
            temp_table = []
            header = True
            for line in f.readlines(): 
                if (line.strip()[0]=='#' or len(line.strip())==0):
                    pass
                else: 
                    if(header): 
                        self.header = re.sub('[ *\n]','', line) # remove spaces
                        header = False 
                    else: 
                        temp_table.append( [int(scr) for scr in re.split(' *', line[1:].strip()) ])
            
            self.table = np.array(temp_table)
            print('header:', self.header)
            print('header length:', len(self.header))
            
            
    def __get_aa_index__(self, aa): 
        return self.header.find(aa)
    
    def get_score(self, aa1, aa2): 
        return self.table[self.__get_aa_index__(aa1), self.__get_aa_index__(aa2)]
    

class sequence_alignment: 
    
    def __init__(self, sub_mat, seqs): 
        
        assert max([len(x) for x in seqs]) == min([len(x) for x in seqs]), "sequence lengths must be the same length to key starting index, add indels to offset lengths to proper alignment" 
              
        self.score_matrix = sub_mat
        self.seqs = seqs
        self.len = len(seqs[0])
        self.num_seqs = len(seqs)
        
        self.__generate_entropy_matrix__()
        
        with open('./data/entropy_matrix.txt', 'w') as f: 
            np.set_printoptions(precision=3,threshold=np.nan)
            f.write(np.array2string(self.__ent_mat__, max_line_width=np.inf))
        
        self.__score_alignment__()
        
        
        
        
    def __generate_entropy_matrix__(self): 
        ''' numpy matrix representing entropy matrix such that row correlates to header amino acid order 
        i  0   1   ...   n   
        A .2   0
        R .4   1
        N .2   0
        D .2   0
        .
        :
        *  0   0
        
        TODO: I think this is printed cols as rows, arbitrary, just note. Also, HEADER IS WRONG NUM OF CHARACTERS, len = 23, should be 22 
        '''
        self.__ent_mat__ = np.zeros([self.len, len(self.score_matrix.header)])
        
        for i in range(self.len): 
            for j in range(self.num_seqs): 
                aa = self.seqs[j][i]
                self.__ent_mat__[i][self.score_matrix.__get_aa_index__(aa)] += 1/self.num_seqs
                
        print(self.num_seqs)
                
                
        
    def __score_alignment__(self): 
        pass
    
    
                     

if __name__ == '__main__' : 
    
    smat = substitution_matrix("./data/BLOSUM62.txt") 
    #print( smat.get_score('S','A') )
    
    seqs = []
    with open('./data/hw3.txt') as f: 
        for line in f.readlines(): 
            if (line.strip()):
                seqs.append(line.strip())
                print(len(line.strip()))
                print(line.strip())

        
    seq_align = sequence_alignment(smat, seqs) 
                

        
        
        
        
        
        
        
        