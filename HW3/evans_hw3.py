# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 19:58:59 2018

@author: nathaniel evans 
@class: BMI650 ALGS
@HW: 3 

The following program can be run for the purpose of scoring any number and length(*see time/mem complexity discussion) of amino acid sequence alignments. The sequence to be scored must be in the directory, relative to this script, ./data/hw3.txt and contain only characters: [ARNDCQEGHILKMFPSTWYVBZX-]. Each sequence should have the same length. 

Scoring can be calculated in two ways, and with the option of including or excluding gap scores, which produces four total output scores. 

The entropy score is produced by using shannon's entropy equation (-sigma: pi*log2(pi) ) for each column and summing over the entire sequence length. This can be done including or excluding indel values BUT columns are not dropped when indels are present, rather, they are excluded in the entopy equation but a normalized column is still used. 
    e.g. (excluding indels) given the column: [K - T], the entropy score would be calculated as: 
            2*(0.5 * log2(0.5))
            This provides a more accurate entropy score by including all information where sequences overlap. 
            
Likewise, sequence score can also be calculated using the sum of pairs method, and this requires an input substitution matrix. This file should be located in ./data/BLOSOM62.txt. For syntax on valid input files, please see example at: https://www.ncbi.nlm.nih.gov/Class/BLAST/BLOSUM62.txt. Naming conventions should be adapted if using non-default substitution matrix. This method can also be calculated including or excluding indels, however, like before this method does not drop the column but instead treats a indel-aminoacid pair as zero scoring, thereby maximizing the amount of information retained in the calculated score. 

These scores are written to console and an output file located: ./data/evans_output.txt

Additionaly, the entropy matrix is printed to file, located in ./data/entropy_matrix.txt which can be a useful tool to visulize the sequence alignment as a frequency distribution. Columns represent amino acids (in the order: [ARNDCQEGHILKMFPSTWYVBZX-]) and rows represent sequence alignment index. As an example, the first row, first column represents the proportion of A's that occur in the first position of the sequence alignment. Please note, if you write the entropy matrix to file AFTER the entropy calculation has been preformed, the indel column may be dropped, but the rows will still be normalized to sum to 1. 

For the end results of this hw assignment, please consider the values 
    entropy score: 66.26
    sum of pairs: 1003 
Please note, this value will be slightly different from many of my classmates because I chose to include the relevant entropy and sum of pair values when only one sequence has an indel. In this case, entropy is calculated as if it was a two pair sequence and only non-indel pairs are summed. For more details on specifics, please inquire with author: evansna@ohsu.edu. 
    
## Time and Memory Complexity discussion ##     
This method of scoring does store sequences as frequency matrices and therefore the memory complexity increases as O(n), if n = seq_length and therefore, memory requirements are not neglible for large sequences. 

This method does have some advantages in time complexity however, since matrix operations are often faster. Time complexity should follow O(n*m^2) or O(n*2^m)-need to think about it a bit more-, if m = # seqs aligned 


"""
import re
import numpy as np

class substitution_matrix: 
    '''
    description
    
    Args:
        input
    Returns:
        output
    '''
    
    def __init__(self, file): 
        '''
        description
        
        Args:
            input
        Returns:
            output
        '''
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
        '''
        description
        
        Args:
            input
        Returns:
            output
        '''
        return self.header.find(aa)
        
    
    def get_score(self, aa1, aa2, score_indel=True): 
        '''
        description
        
        Args:
            input
        Returns:
            output
        '''
        #print( (aa1, aa2) )
        if ( not score_indel and (aa1=='*' or aa2=='*')) :
            return 0
            
        else:
            return self.table[self.__get_aa_index__(aa1), self.__get_aa_index__(aa2)]

        
    
    

class sequence_alignment: 
    '''
    description
    
    Args:
        input
    Returns:
        output
    '''
    
    def __init__(self, sub_mat, seqs): 
        '''
        description
        
        Args:
            input
        Returns:
            output
        '''
        
        assert max([len(x) for x in seqs]) == min([len(x) for x in seqs]), "sequence lengths must be the same length to key starting index, add indels to offset lengths to proper alignment" 
              
        self.score_matrix = sub_mat
        self.seqs = np.array( [[aa for aa in seq] for seq in seqs] )
        self.len = len(seqs[0])
        self.num_seqs = len(seqs)
        
        self.__generate_entropy_matrix__()
    
    def write_ent_matrix_to_file(self, name='entropy_matrix.txt'): 
        '''
        description
        
        Args:
            input
        Returns:
            output
        '''
        with open('./data/' + name, 'w') as f: 
            np.set_printoptions(precision=3,threshold=np.nan)
            f.write(np.array2string(self.__ent_mat__, max_line_width=np.inf))
        
    def score(self, method='entropy', include_gaps=True): 
        '''
        description
        
        Args:
            input
        Returns:
            output
        '''
        assert method == 'entropy' or method == 'sum_of_pairs', "unrecognized scoring method" 
        
        if (method=='entropy') : 
            self.__score_entropy_matrix__(include_indels=include_gaps)
            return self.entropy_sum 
        
        else: # sum of pairs 
            return self.__score_sum_of_pair__(score_indel=include_gaps)
        
    def __pair_sum__(self, seqs_col, score_indel): 
        '''
        description
        
        Args:
            input
        Returns:
            output
        '''
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
        '''
        description
        
        Args:
            input
        Returns:
            output
        '''
        
        score = 0
        for i in range( self.len ): 
            col = self.seqs[:,i]
            score += self.__pair_sum__(col, score_indel)
            
        return score
    
    
                # make sure to talk about including the indel with one gap in scoring (even with include_indel == false). Entropy created by indel is NOT included, but entropy created from the two other mismatches aas ARE included. 
    def __score_entropy_matrix__(self, include_indels=True): 
        '''
        description
        
        Args:
            input
        Returns:
            output
        '''
        
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
        '''
        description
        
        Args:
            input
        Returns:
            output
        '''
        if (pi > 0 and pi < 1):
            return -pi*np.log2(pi)
        else:
            return 0 
        
        
    def __generate_entropy_matrix__(self):
        '''
        numpy matrix representing entropy matrix such that row correlates to header amino acid order 
        i  0   1   ...   n   
        A .2   0
        R .4   1
        N .2   0
        D .2   0
        .
        :
        *  0   0
    
        Args:
            input
        Returns:
            output
        '''

        self.__ent_mat__ = np.zeros([self.len, len(self.score_matrix.header)])
        
        for i in range(self.len): 
            for j in range(self.num_seqs): 
                aa = self.seqs[j,i]
                self.__ent_mat__[i, self.score_matrix.__get_aa_index__(aa)] += 1/self.num_seqs
                
                     
'''
    description
    
    Args:
        input
    Returns:
        output
    '''
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
    seq_align.write_ent_matrix_to_file()
                
    es_g = seq_align.score(method='entropy', include_gaps=True)
    print('entropy score with gaps included:', es_g)
    
    es_ng = seq_align.score(method='entropy', include_gaps=False) # this actually alters the object, future calls will not be able to include gaps without reinitializing the object. 
    print('entropy score, gaps not included:',es_ng)
    
    sop = seq_align.score(method='sum_of_pairs', include_gaps=True) 
    print('sum of pairs gaps included:',sop )
    
    sop_ng = seq_align.score(method='sum_of_pairs', include_gaps=False) 
    print('sum of pairs gaps not included:',sop_ng )
        
    with open('./data/evans_output.txt','w') as f: 
        f.write('entropy score with gaps included:' + str(es_g))
        f.write('entropy score, gaps not included:' + str(es_ng))
        f.write('sum of pairs gaps included:' + str(sop)) 
        f.write('sum of pairs gaps not included:' + str(sop_ng ))
        
        
        
         