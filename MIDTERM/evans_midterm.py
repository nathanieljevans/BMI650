# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 17:08:23 2018

@author: natha


"pseduocode" for BMI650 midterm 

"""
import numpy as np 
import re
import networkx as nx 
from matplotlib import pyplot as plt

S1 = 'GWWPDT' 
S2 = 'WRRKHY'

    
 #### PROBLEM 1 ####   
#----------------------------------------------------------

def generate_BW_transform(seq): 
    '''
    Slightly modified from HW4
    This function generates a Burrows Wheeler transform data structure, which represents the last column of a BW matrix for a given sequence. 
    
    input
        seq <string> sequence to be transformed 
        
    output 
        BW <list> last column of the BW matrix. 
    '''
    
    mixed = [seq[i:] + seq[0:i] for i in range(len(seq))]
    mixed.sort()
    BW = np.array([[l for l in s] for s in mixed]) 
    print(BW)
    return BW[:,0], BW[:,-1]

#use generator to efficiently rebuild, keeps function state between calls
def  build_row(row_index, FC, LC): 
    nxt = FC[row_index] 
    yield nxt

    b_rank = ''.join(FC)[0:row_index].count(nxt) 
    while (nxt != '$'): 
        nxt, b_rank = find_next(nxt, b_rank, FC, LC) 
        yield nxt 
    
    while (True):
        yield '$'
		
def find_next(cur, b_rank, FC, LC): 
    '''
    Slightly modified from HW4
    finds the next character in the given suffix (eg row), as defined by cur and b_rank. 
    
    inputs
        cur <char> single char representing the last char in suffix 
        b_rank <int> the occurence of the cur char in the FC array 
        FC <list> list of chars representing first column of burrows-wheeler matrix 
        LC <list> list of chars representing last column of burrows-wheeler matrix 
        
    ouputs 
        nxt <char> next char in given suffix
        new_b_rank <int> the occurence count of nxt in the FC 
    '''
    
    nxt_row = [m.start() for m in re.finditer(cur,''.join(LC))][b_rank] # find b_rank-ith occurence of cur in the LC 
    nxt = FC[nxt_row]
    new_b_rank = ''.join(FC)[0:nxt_row].count(nxt) 
    
    return nxt, new_b_rank
	
def check_contains_unique_string(rows): 
    not_unique = set()
    for i, s1 in enumerate(rows[:-1]): 
        is_unique = True

        for j, s2 in enumerate(rows[i+1:]): 
            if ( s1 == s2 or s1 in not_unique or ('$' in s1)): 
                is_unique = False
                not_unique.add(s1)
                break
            
        if (is_unique): 
            return True
    return False 
				

def search_for_shortest_unique_substring_length(S) : 
    FC, LC = generate_BW_transform(S) 
    row_generators = [build_row(i, FC, LC) for i in range(len(S))] # generators for each row
    rows = [next(gen) for gen in row_generators] # first letter of each row stored in list, index -> row index
    
    for i in range(len(S)): 
        if (check_contains_unique_string(rows)) : 
            return i + 1
        
        for j in range(len(rows)): 
            rows[j] += next(row_generators[j])
    
    return (len(S)) # no shortest substring found 
#-------------------------------------------------------
    #### PROBLEM 6 #### 
#-------------------------------------------------------------------------------
    
def is_circular_rotation(s1,s2): 
    for i in range(len(s1)): 
        if (s1[i:] + s1[:i] == s2): 
            return True 
        
    return False 

#-------------------------------------------------------------------------------
#Runs in O(n) where n is the length of S1/S2 
#-------------------------------------------------------------------------------
   #### PROBLEM 7 ####
'''
def rec_search(node, tree): 
	count = 0 
	If (node.has_no_children()): 
		return 1 
	
	For child in node.get_children(): 
		count += rec_search (child, tree) 
	
	return count 
-------------------------------------------------------------------------------
num_distinct_substrings = rec_search(root, tree) 
'''	

#-------------------------------------------------------------------------------



#-------------------------------------------------------

if __name__ == '__main__' : 

    shortest_unique_str_len = search_for_shortest_unique_substring_length('ATGATGATG$')
    print('answer:', shortest_unique_str_len)   
    
    s1 = 'TTGATC' 
    s2 = 'ATCTTG'
    
    print('circ rot?', is_circular_rotation(s1,s2))
    
    G = nx.Graph() 
    G.add_node(0)
    G.add_node(1) 
    G.add_edge(0,1, label='NA')
    G.add_node(2)
    G.add_edge(1,2, label='$')
    G.add_node(3)
    G.add_edge(1,3, label='$')
    #nx.draw(G, with_labels=True, edge)

    
    nx.draw_networkx(G)
# Draw the node labels
    nx.draw_networkx_labels(G)
# Draw the edges
    nx.draw_networkx_edges(G)
# Draw the edge labels
    nx.draw_networkx_edge_labels(G, edge_labels=nx.get_edge_attributes(G, 'label'))
    # Remove the axis
    plt.axis('off')
    # Show the plot
    plt.show()
        
#---------------------------------------------------------
        
        
        


