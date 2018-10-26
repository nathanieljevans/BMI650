# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 20:43:09 2018

@author: nathaniel evans 
@class: BMI650
@HW: 4

This script is used to generate a BW transform and suffix array for a given string as well as show a rudimentary method of subsequence search of a BW transform. It can be run by navigating to appropriate directory and running: 
    $python evans_hw4.py 
    
    To adjust input sequence and search string, modify the global constants at the top of this script.

After producing this, there are a variety of optimizations and tradeoffs that could justify modification. As discussed in class, memory vs computational tradeoffs play a significant component. 

This search method will return a given string independant of whether it is a suffix, prefix or midstring. The produced offset represents the index (1-indexed) where the matched substring occurs. 

Future enhancements: 
    Generate FC/LC b-rank array to avoid counting occurences at each next character step 
    Use suffix array index array to get offsets 
    Hybridization approach to partially precomputed arrays above 
    Cleaner code approach

"""
# input sequence and substring search values 
DNA = 'ACTGCTCGGCT$'
SEARCH = 'GCT'


import numpy as np
import re


def generate_BW_transform(seq): 
    '''
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
    return BW[:,-1]
        
def generate_suffix_array(seq):
    ''' 
    This function generates a suffix array of offsets for a given sequence
    
    input 
        seq <string> sequence to generate suffix array from
    
    output
        SFX <list> suffix array offsets 
    '''
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
    ''' 
    This function searches a BW transform data structure for a given string, if one or more is found, the offset is printed to terminal.
    
    input 
        BW <list> list of characters representing the last column of a Burrows Wheeler transform. 
        SFX <list> list of ints representing suffix offsets, indexing matches BW. Unused in this case, but should it be, it provides a computationally more efficient method of finding the offset (when the matched string isn't a true suffix eg mid-string or prefix)
        match<string> the input string to search for 
    output 
        None 
    '''
    LC = BW 
    FC = np.sort(BW)
    
    b_rank = 0
    matched_i = 1
    for row, start in enumerate(FC): 
        if (start == match[0]): 
            nxt = start
            b_rank2=b_rank # b_rank keeps track of starting char b_rank, while b_rank2 is used in the find next trace back
            matched_i = 1 #index of char we're trying to match next 
            
            while (True): # messy but functional
                nxt, b_rank2 = find_next(nxt,b_rank2, FC, LC)   
                
                if (matched_i > len(match)): 
                    break # longer than match length, stop searching
                
                if nxt == match[matched_i]: 
                    if matched_i is len(match)-1: 
                        i_plus = get_index_of_terminal_char(nxt, b_rank2, FC,LC)
                        offset = len(FC) - (len(match) + i_plus)
                        print('match at offset:', offset)
                        break
                    else: 
                        pass # keep searching 
                    matched_i += 1 
                else: 
                    break # not the string!  
        
            b_rank += 1
            
    
def get_index_of_terminal_char(cur, b_rank, FC, LC): 
    '''
    This is used to get the index of the terminal character in a suffix, which allows us to find the offset of a matched string without using the suffix array. Computationally more expensive, but memory wise cheaper. 
    
    inputs
        cur <char> single char representing the last char in suffix 
        b_rank <int> the occurence of the cur char in the FC array 
        FC <list> list of chars representing first column of burrows-wheeler matrix 
        LC <list> list of chars representing last column of burrows-wheeler matrix 
    
    outputs
        i<int> the number of characters after cur to the terminal character. 
    '''
    cur, b_rank = find_next(cur, b_rank, FC, LC) # if next char is terminal, return 0
    i = 0

    while cur != '$' : 
        cur, b_rank = find_next(cur, b_rank, FC, LC)
        i+=1 
    
    return i

def find_next(cur, b_rank, FC, LC): 
    '''
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


if __name__ == '__main__' : 
    
    print('Sequence to generate BW transform and Suffix array from:', DNA)
    print('Subsequence to search for:', SEARCH)
    print('\n') 
    
    BW = generate_BW_transform(DNA) # BW is the LC of the BW matrix 
    print('\n')
    SFX = generate_suffix_array(DNA) 
    print('\n')
    search(BW, SFX, SEARCH) # didn't use the SFX array for offset, but could/should have. Would have been cleaner. 
    