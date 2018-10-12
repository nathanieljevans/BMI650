# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 19:29:30 2018

@author: natha
"""
import numpy as np

#seq 2 across top
SEQ2 = "GWWPDT"
SEQ1 = "WRRKHY" 
SHAPE = (len(SEQ1), len(SEQ2))

HP = {'G','W','T','K','H','Y'}

def make_mat(shape): 
    # even indices represent edges (starting at 0), odd represent nodes 
    # dim has to be double seq lens then
    
    return np.zeros( (shape[0]*2 + 1, shape[1]*2 + 1) )



'''i 0   1    2   3 
 0   *0* |-2|*-2*|-2|...
 1   -2  |-5| -2 |-2|... 
 2   *-2*|-2|*-4*|-2|... 
 j
 
 nodes are when i & j are even 
 edges are all other times 
     indel when i=odd, j=even
     indel when i=even, j=odd 
     diag when i=odd, j=odd 
 
 
'''
def get_score(i, j, seq1, seq2): 
    
    try: 
        s1 = seq1[int((i-1)/2) ]
        s2 = seq2[int((j-1)/2) ]
    except: 
        return -999
        raise
    
    if (s1 == s2): 
        print([i,j, s1,s2])
        return 5 # match
    elif(s1 in HP and s2 in HP): 
        return 1 # both HP
    elif(s1 in HP or s2 in HP): 
        return -5 # only 1 in HP 
    else: 
        return 0 # neither are HP

# cartesian distance of a tuple - for sorting
def dist(tup): 
    return ( (tup[0]**2 + tup[1]**2)**(0.5) ) 

if __name__ == "__main__" :
    print(SHAPE)
    board = make_mat(SHAPE)
    
    # populate indels match/mismatch - diag transitions
    # -2 gap penalty 
    for i in range(0, board.shape[0], 1): 
        for j in range(0, board.shape[1], 1): 
            if (i==0 and j==0):
                board[i,j]=0
            # indel position 
            
            elif ( (i%2==1 and j%2==0) or (i%2==0 and j%2==1) ):
                board[i][j] = -2
            
            elif (i%2==1 and j%2==1): 
                board[i][j] = get_score(i, j, SEQ1, SEQ2) # diags

    #print(np.array2string(board, max_line_width=np.inf))             
    #print()
        # create node indices, then sort by distance from 0 to create proper ordering 
    node_order = []
    for i in range(0,board.shape[0],2): 
        for j in range(0, board.shape[1], 2): 
            node_order.append( (i,j) )
    
    # sort by distance from origin
    node_ordered = sorted(node_order, key=dist)
    # remove 0,0 
    node_ordered.remove((0,0))
    
    # populate nodes in proper order 
    for node in node_ordered: 
        i,j = node
        max_ls = []
        if (j>1): 
            max_ls.append( board[i][j-2] + board[i][j-1] )
        if (i>1): 
            max_ls.append( board[i-2][j] + board[i-1][j] ) 
        if (j>1 and i>1): 
            max_ls.append( board[i-2][j-2] + board[i-1][j-1] )
        board[i][j] = np.max(max_ls) 
     
    print(np.array2string(board, max_line_width=np.inf))
    
    # populate match/mismatch - diag transitions
    
    
    
    
    
    