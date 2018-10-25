# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 20:43:09 2018

@author: nathaniel evans 
@class: BMI650
@HW: 4
"""

import numpy as np
import re

DNA = 'ACTGCTCGGCT$'


def generate_BW_transform(seq): 
    
    mixed = [seq[i:] + seq[0:i] for i in range(len(seq))]
    mixed.sort()
    BW = np.array([[l for l in s] for s in mixed]) 
    print(BW)
    return BW[:,-1]
        
def generate_suffix_array(seq):
    
    mixed = [seq[i:] + seq[0:i] for i in range(len(seq))]
    mixed.sort() 
    
    SFX = []
    for s in mixed: 
        i = s.index('$')
        offset = len(seq) - i
        suffix = s[0:i]
        gap = '' 
        if offset < 10 : 
            gap = '0'
        SFX.append( (offset, suffix) ) 
        print(gap + str(offset) + ' || ' + suffix)
        
    return SFX

def search(BW, SFX, match): 
    LC = BW 
    FC = np.sort(BW)
    print(FC)
    print(match)
    
    b_rank = 0
    matched_i = 1
    for row, start in enumerate(FC): 
        if (start == match[0]): 
            nxt = start
            print('row', row)
            while (True): 
                nxt, b_rank = find_next(nxt,b_rank, FC, LC)   
                print('compare:', (nxt, match[matched_i]))
                
                if (matched_i > len(match)): 
                    print('longer than match length')
                    break
                
                if nxt == match[matched_i]: 
                    if matched_i is len(match)-1: 
                        print('A FULL MATCH, row:', row)
                        break
                    else: 
                        print('its a match! passing!')
                        pass # keep searching 
                    matched_i += 1 
                else: 
                    print('letter does not match, not end')
                    break # not the string!  
        
            
            
            
    
    print(FC)

def find_next(cur, b_rank, FC, LC): 
    
    #cur = FC[row]
    #b_rank = ''.join(FC)[0:row].count(cur) # could pass b_rank around instead of row to make it computationally cheaper, but intuitvely it makes more sense to me... for conceptual value 
    nxt_row = [m.start() for m in re.finditer(cur,''.join(LC))][b_rank] # find b_rank-ith occurence of cur <-------------------- this may not be working quite right, or might need to add +1 ... def not working right. 
    
    print('next row',nxt_row)
    nxt = FC[nxt_row]
    
    new_b_rank = ''.join(FC)[0:nxt_row].count(nxt) 
    print('nxt, new_b_rank', (nxt, new_b_rank) )
    
    return nxt, new_b_rank


if __name__ == '__main__' : 
    
    BW = generate_BW_transform(DNA)
    
    SFX = generate_suffix_array(DNA)
    
    print(BW)
    
    search(BW, SFX, 'GCT')
    