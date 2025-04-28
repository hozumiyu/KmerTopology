# -*- coding: utf-8 -*-
"""
Created on Sun Apr 27 15:16:19 2025

@author: yutah
"""

from itertools import product
import numpy as np

def get_kmers_como(kmers_size):
    #obtain all the k-mers combination
    #return the list of kmers
    #ex: if kmers_size = 2, kmers_list = [AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT]
    temp_list  = list(product('ACGT', repeat= kmers_size))
    kmers_list = []
    for t in temp_list:
        kmers_list.append(''.join(t))
    kmers_list.sort()
    return kmers_list

def find_kmers_position(sequence, kmers_size):
    '''
        Extract all the positions of the k-mers

    Parameters
    ----------
    sequence : str
        Nucleotide sequence. Please convert all uracil (U) to T, and make the sequence a capital letter
        for now, invalid nucleotides are ignored.

    Returns
    -------
    None.

    '''
    kmers_list = get_kmers_como(kmers_size)
    length_seq = len(sequence) - kmers_size + 1
    
    
    kmers_pos = { k: [] for k in kmers_list  }
    for idx in range(length_seq):
        s = sequence[idx:idx+kmers_size]
        if s in kmers_list:
            kmers_pos[s].append(idx+1)
        else:
            print(s, 'not valid dna')  #<- comment this out if you dont want this message. Won't affect your result.
        
    for kmers in kmers_list:
        kmers_pos[kmers] = np.array(kmers_pos[kmers]).astype(float)
    return kmers_list, kmers_pos