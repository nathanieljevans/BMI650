# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 17:08:23 2018

@author: natha


"pseduocode" for BMI650 midterm 

"""
import numpy as np 
import re
import networkx as nx 
from networkx.drawing.nx_agraph import write_dot, graphviz_layout

#----------------------------------------------------------   
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
            if (s1 == s2 or s1 in not_unique or ('$' in s1)): 
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
    #### PROBLEM 7 ####
    
def build_suffix_trie(S): 
    S += '$' 
    G = nx.DiGraph() 
    G.add_node(0, label='root')
    
    suffix = '' 
    node_id = 1
    labels = {0:'root'}
    for i in range(1, len(S)+1): 
        cur_node = 0
        suffix = S[(len(S) - i):]
        print(suffix)
        
        for c2 in suffix: 
            nxt_edge = has_edge_char(G, cur_node, c2)
            if (nxt_edge): 
                cur_node = nxt_edge
                
            else: 
                G.add_node(node_id, label = c2)
                G.add_edge(cur_node, node_id)
                labels[node_id] = c2 
                cur_node = node_id
                node_id += 1 
                
    return G, labels

# TODO -------------------------------------------
def convert_trie_to_tree(G): 
    
    G2 = nx.DiGraph() 
    id_ = object()
    id_.val  = 0
    cur_node = 0
    
    # traverse to end, fall back to branch pt and make node with all
        
    # There is a clever way to rebuild this recursively but it's eluding me right now
def rebuild(G, node, _id, last_node, s): 
    ney = [x for x in G.neighbors(node)]
    if (len(ney) == 0): 
        return G.nodes[node]['label']
    elif (len(ney) == 1): 
        return G.nodes[node]['label'] + rec_rebuild(G, ney[0], _id)
    else: 
    
        for node in ney: 
            _id.val+=1
            s = rec_rebuild(G, node, _id)
# ------unfinished ---------------------------------            
    

def has_edge_char(G, node, c): 
    for node in G.neighbors(node): 
        if (G.nodes[node]['label'] == c) : 
            return node 
            
    return None 

def rec_search(node, tree): 
    if (node != 0): # don't include root
        count = len([x for x in tree.nodes[node]['label'] if x != '$']) # don't include terminal characters in count, either in it's own node, or with others 
    else: 
        count = 0
	
    for child in G.neighbors(node): 
        count += rec_search (child, tree) 
	
    return count




#-------------------------------------------------------

if __name__ == '__main__' : 

    shortest_unique_str_len = search_for_shortest_unique_substring_length('ATGATGATG$')
    print('answer:', shortest_unique_str_len)   
    
    s1 = 'TTGATC' 
    s2 = 'ATCTTG'
    
    print('circ rot?', is_circular_rotation(s1,s2))
    
    G, labels = build_suffix_trie('AATATT')
    
    num_distinct_substrings = rec_search(0, G) 
    print('num distinct substrings in suffix:', num_distinct_substrings)

    #pos= graphviz_layout(G, prog='dot')
    nx.draw(G, labels = labels, with_label=True)
'''

'''     
#---------------------------------------------------------
        
        
        


