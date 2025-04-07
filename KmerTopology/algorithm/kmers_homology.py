# -*- coding: utf-8 -*-
"""
Created on Sun Feb 11 15:03:42 2024

@author: yutah
"""

import ripser
import gudhi as gd
import numpy as np
from gudhi.representations.vector_methods import BettiCurve
from itertools import product
from scipy import sparse

'''
This code is used to generate persistent homology feature
You will need ripser and Gudhi to compute this
'''



class KmersHomology():
    def __init__(self, kmers_size = 2, thresh  = 0):
        '''
        Parameters
        ----------
        kmers_size : int, optional
            The k-mers size. You will have a topotal of 4^k different k-mer of size k. The default is 2.
        max_dimension : TYPE, optional
            DESCRIPTION. The default is 0.
        thresh : int, optional
            Maximum Filtration radius. The default is -1. If you sit it to -1, It will go to 4^(k-1)*50 steps.
            If you have this set to -1, maximum filtration is going to be the length of you DNA sequence (not recommended)

        Returns
        -------
        None.

        '''
        self.kmers_size = kmers_size
        self.kmers_list = self.get_kmers_como()
        if thresh == 0:
            self.thresh = 4**(kmers_size-1)*50
        elif thresh > 0:
            self.thresh = thresh
        elif thresh == -1:
            self.thresh = -1
            
            
    def get_kmers_como(self):
        #obtain all the k-mers combination
        #return the list of kmers
        #ex: if kmers_size = 2, kmers_list = [AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT]
        
        temp_list  = list(product('ACGT', repeat= self.kmers_size))
        kmers_list = []
        for t in temp_list:
            kmers_list.append(''.join(t))
        kmers_list.sort()
        return kmers_list

    def find_kmers_position(self, sequence):
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
        length_seq = len(sequence) - self.kmers_size + 1
        
        
        self.kmers_pos = { k: [] for k in self.kmers_list  }
        for idx in range(length_seq):
            s = sequence[idx:idx+self.kmers_size]
            if s in self.kmers_list:
                self.kmers_pos[s].append(idx+1)
            else:
                print(s, 'not valid dna')  #<- comment this out if you dont want this message. Won't affect your result.
            
        for kmers in self.kmers_list:
            self.kmers_pos[kmers] = np.array(self.kmers_pos[kmers]).astype(float)
        return 
    
    
    
 
    def compute_kmers_persistent_diagram(self, kmers = 'AT'):
        '''
            Compute the persistnet diagram for a given kmer
            
        Parameters
        ----------
        kmers : str
            The k-mer you want. This needs to match the kmers_size you initialize with

        Returns
        -------
        pd_kmers : np.array
            persistent diagram. the matrix is M by 2, where M is the number of kmer.
            pd_kmers[:, 0] is the birth time, pd_kmers[:,1] is the death time. 
            Note that pd_kmers[:, 0] are all 0 because we are only doing 0-th order.
        '''
            
        
        pos = self.kmers_pos[kmers][self.kmers_pos[kmers] > 0]  #only keep the nonzero entry
        
        #dis = sparse.coo_array((pos.shape[0], pos.shape[0]))
        row = []; col = []; val = []
        if pos.shape[0] > 1:
            if self.thresh == -1:
                thresh = np.max(pos[1:] - pos[:-1]) + 1
            else:
                thresh = self.thresh
            
            for current_idx in range(pos.shape[0]):
                current_position = pos[current_idx].copy()  #get the current position
                
                current_distance = np.abs(current_position - pos)#distance btween current position and all the position
                index_threshold = np.where(current_distance < thresh)[0]   #extract all the index satisfying the threshold criteria
                
                for index in index_threshold:
                    if index != current_idx:
                        row.append(current_idx)
                        col.append(index)
                        val.append(current_distance[index])
            row = np.array(row); col = np.array(col); val = np.array(val)
            dis = sparse.coo_array( (val, (row, col)), shape = (pos.shape[0], pos.shape[0]))
            pd_kmers = ripser.ripser(dis, thresh = thresh, distance_matrix = True)['dgms'][0]
        elif pos.shape[0] == 1:
            pd_kmers = np.array([[0, np.inf]])
        elif pos.shape[0] == 0:
            pd_kmers = np.array([[0, 0]])
        return pd_kmers
    
    

    
    def compute_kmers_betti(self, pd_kmers):
        '''
        

        Parameters
        ----------
        pd_kmers : np.array
            Persistent diagram

        Returns
        -------
        betti_curve : 1d vector of np.array
            Betti curve.

        '''
        if pd_kmers.ndim == 1:
            if pd_kmers[1] == 0:
                betti_curve = np.array([0,0])
            else:
                betti_curve = np.array([1,1])
        else:
            sorted_kmers = pd_kmers[:, 1].copy()
            sorted_kmers[sorted_kmers == np.inf] = 0
            max_filtration = np.max(sorted_kmers)
            bc = BettiCurve(predefined_grid = np.arange(max_filtration+2))
            betti_curve = bc.fit_transform([pd_kmers]).reshape(-1)
        return betti_curve

