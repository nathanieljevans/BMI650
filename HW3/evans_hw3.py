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

if columns have a gap, then ignore it 



#  Matrix made by matblas from blosum62.iij
#  * column uses minimum score
#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
#  Entropy =   0.6979, Expected =  -0.5209
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 



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
                        self.header = re.sub('[ \n]*','', line) # remove spaces
                        header = False 
                    else: 
                        temp_table.append( [int(scr) for scr in re.split(' *', line[1:].strip()) ])
                        
            self.table = np.array(temp_table)
            
            with open('./data/substitution_matrix.txt', 'w') as f: 
                np.set_printoptions(precision=3,threshold=np.nan)
                f.write(np.array2string(self.table, max_line_width=np.inf))
            
            
    def __get_aa_index__(self, aa): 
        return self.header.find(aa)
        
    
    def get_score(self, aa1, aa2, score_indel=True): 
        #print( (aa1, aa2) )
        if ( not score_indel and (aa1=='*' or aa2=='*')) :
            return 0
            
        else:
            return self.table[self.__get_aa_index__(aa1), self.__get_aa_index__(aa2)]

        
    
    

class sequence_alignment: 
    
    def __init__(self, sub_mat, seqs): 
        
        assert max([len(x) for x in seqs]) == min([len(x) for x in seqs]), "sequence lengths must be the same length to key starting index, add indels to offset lengths to proper alignment" 
              
        self.score_matrix = sub_mat
        self.seqs = np.array( [[aa for aa in seq] for seq in seqs] )
        self.len = len(seqs[0])
        self.num_seqs = len(seqs)
        
        self.__generate_entropy_matrix__()
        
#        with open('./data/entropy_matrix.txt', 'w') as f: 
#            np.set_printoptions(precision=3,threshold=np.nan)
#            f.write(np.array2string(self.__ent_mat__, max_line_width=np.inf))
#        
#        with open('./data/entropy_matrix2.txt', 'w') as f: 
#            np.set_printoptions(precision=3,threshold=np.nan)
#            f.write(np.array2string(self.__ent_mat__, max_line_width=np.inf))
        
    def score(self, method='entropy', include_gaps=True): 
        assert method == 'entropy' or method == 'sum_of_pairs', "unrecognized scoring method" 
        
        if (method=='entropy') : 
            self.__score_entropy_matrix__(include_indels=include_gaps)
            return self.entropy_sum 
        
        else: # sum of pairs 
            return self.__score_sum_of_pair__(score_indel=include_gaps)
        
    def __pair_sum__(self, seqs_col, score_indel): 
        score = 0
        if len(seqs_col) > 1 :
            aa1 = seqs_col[0]
            for i,aa2 in enumerate(seqs_col[1:]): 
                score += self.score_matrix.get_score(aa1,aa2,score_indel) + self.__pair_sum__(seqs_col[i+1:], score_indel)
            return score
        else: 
            return 0
                
    
    # this also includes cols with indels if there are non indel pairs to compare 
    def __score_sum_of_pair__(self, score_indel): 
        
        score = 0
        for i in range( self.len ): 
            col = self.seqs[:,i]
            score += self.__pair_sum__(col, score_indel)
            
        return score
    
    
                # make sure to talk about including the indel with one gap in scoring (even with include_indel == false). Entropy created by indel is NOT included, but entropy created from the two other mismatches aas ARE included. 
    def __score_entropy_matrix__(self, include_indels=True): 
        
        ent_vec = np.vectorize(self.__entropy__)
        
        if(not include_indels): 
            # remove indel col (last )
            self.__ent_mat__ = np.delete(self.__ent_mat__, -1, 1)
            # normalize each row 
            #x = np.array([LA.norm(v,ord=1) for v in X]
            norm = np.transpose( np.tile((np.linalg.norm(self.__ent_mat__, axis=1, ord=1)) , (23,1) ) )
            
            self.__ent_mat__ = self.__ent_mat__ / norm # normalize rows of the indel-less entropy matrix
        
        self.entropy_sum = np.sum( ent_vec(self.__ent_mat__) ) 

        
    def __entropy__(self, pi): 
        if (pi > 0 and pi < 1):
            return -pi*np.log2(pi)
        else:
            return 0 
        
        
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
        
        '''
        self.__ent_mat__ = np.zeros([self.len, len(self.score_matrix.header)])
        
        for i in range(self.len): 
            for j in range(self.num_seqs): 
                aa = self.seqs[j,i]
                self.__ent_mat__[i, self.score_matrix.__get_aa_index__(aa)] += 1/self.num_seqs
                
                     

if __name__ == '__main__' : 
    
    smat = substitution_matrix("./data/BLOSUM62.txt") 
    #print( smat.get_score('S','A') )
    
    seqs = []
    with open('./data/hw3.txt') as f: 
        for line in f.readlines(): 
            if (line.strip()):
                test = re.sub('-', '*', line.strip())
                seqs.append(re.sub('-', '*', line.strip()))
                

        
    seq_align = sequence_alignment(smat, seqs) 
                
    es_g = seq_align.score(method='entropy', include_gaps=True)
    print('entropy score with gaps included:', es_g)
    
    es_ng = seq_align.score(method='entropy', include_gaps=False) # this actually alters the object, future calls will not be able to include gaps without reinitializing the object. 
    print('entropy score, gaps not included:',es_ng)
    
    sop = seq_align.score(method='sum_of_pairs', include_gaps=True) 
    print('sum of pairs gaps included:',sop )
    
    sop_ng = seq_align.score(method='sum_of_pairs', include_gaps=False) 
    print('sum of pairs gaps not included:',sop_ng )
        
        
        
        
        
        
        