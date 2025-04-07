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





class KmersCircularLaplacian():
    def __init__(self, kmers_size = 2, steps = np.array([0,1,2])):
        #kmers_size: for 2mers, enter 2
        # max_dimension: the dimension of the simplicial complex
        # thresh: maximum edge length, default value is 4**kmersize * 10
        self.kmers_size = kmers_size
        self.kmers_list = self.get_kmers_como()
        self.steps = steps
        
    def get_kmers_como(self):
        #obtain all the kmers
        temp_list  = list(product('ACGT', repeat= self.kmers_size))
        kmers_list = []
        for t in temp_list:
            kmers_list.append(''.join(t))
        kmers_list.sort()
        return kmers_list

    def find_kmers_position(self, sequence):
        #extract the position of the nucleotide, starting at 1
        #set the position to 0 if the nucleotide is not 'A', for example
        #return the positon of each nucleotide, and the masking vector as a dictionary
        length_seq = len(sequence) - self.kmers_size + 1
        self.length_seq = length_seq
        
        self.kmers_pos = { k: [] for k in self.kmers_list  }
        for idx in range(length_seq):
            s = sequence[idx:idx+self.kmers_size]
            if s in self.kmers_list:
                self.kmers_pos[s].append(idx+1)
            else:
                print(s, 'not valid dna')
            
        for kmers in self.kmers_list:
            self.kmers_pos[kmers] = np.array(self.kmers_pos[kmers]).astype(float)
        return 
    
    
    
 
    def compute_kmers_Distance(self, kmers = 'AT'):
        '''
        Compute the distance between the kmers

        Parameters
        ----------
        kmers : str, optional
            Kmers type. Make sure this matches the kmers_size you initialized with. The default is 'AT'.

        Returns
        -------
        D : np.array
            distance matrix for kmers

        '''
        pos = self.kmers_pos[kmers][self.kmers_pos[kmers] > 0]  #only keep the nonzero entry
        #define the threshold
        
        D = np.zeros([pos.shape[0], pos.shape[0]])
        for idx in range(pos.shape[0]):
            for idx2 in range(idx+1, pos.shape[0]):
                p1 = pos[idx2] - pos[idx]
                p2 =  self.length_seq - pos[idx2] + pos[idx]
                d = np.min( [ p1,p2] )  # check which direction yields the smallest distance
                D[idx,idx2]= d
                D[idx2, idx] = d
        return D
    
    
    
    
    def computeLaplacianEig(self, D):
        '''
            Compute all the eigenvalues associated with the filtration (step)

        Parameters
        ----------
        D : np.array
            distance matrix for kmers

        Returns
        -------
        eigenval_list : list
            list all the eigenvalues. The first few should be 0, and the number corresponds to the Betti 0, or the number of connected components
            

        '''
        eigenval_list = []
        for idx, s in enumerate(self.steps):
            if s == 0:
                eig = np.zeros(D.shape[0])
                eigenval_list.append(eig)
            else:
                A = D.copy()
                A[A>s] = 0
                A[A!=0] = 1
                L = np.diag(np.sum(A,axis = 0)) - A
                eig = np.linalg.eigvalsh(L)
                eigenval_list.append(eig)
        return eigenval_list
        
    

    
    def compute_kmers_betti_min(self, eigenval_list):
        '''
        Compute the betti number and the smallest nonzero eigenvalue

        Parameters
        ----------
        eigenval_list : list
            list all the eigenvalues. The first few should be 0, and the number corresponds to the Betti 0, or the number of connected components

        Returns
        -------
        Betti curve: length is determined by the filtration size (steps)
        Minimal_eigenvalue_curve: length is determined by the filtration size (steps)
            Note that if all your points are isolated, then the value will be set to 0 by default (radius = 0 will always yield 0)

        '''
        betti = []
        eig_min = []
        for eig in eigenval_list:
            numzero = np.where(eig <1e-6)[0].shape[0]
            betti.append(numzero)
            if numzero == len(eig):
                eig_min.append(0)
            else:
                eig_min.append(eig[numzero])
        return np.array(betti), np.array(eig_min)



'''
#test_dna = 'CGGATAACGTCCAGCAGTCAGTGATCGCATATCTTGAC'
np.random.seed(1)
length_seq = 1000
sequence = np.random.choice(['A', 'C', 'G', 'T'], length_seq)
sequence = ''.join(sequence)

step = np.arange(50)
kmers_size = 1
kmers = 'A'

myKmersLaplacian = KmersLaplacian(kmers_size = kmers_size, steps= step)
myKmersLaplacian.find_kmers_position(sequence)
kmers_list = myKmersLaplacian.kmers_list


import time
start = time.time()
D = myKmersLaplacian.compute_kmers_Distance(kmers)
eigenval_list = myKmersLaplacian.computeLaplacianEig(D)
betti_curve_full, eig_min = myKmersLaplacian.compute_kmers_betti_min(eigenval_list)
end = time.time()
print('Full ripser:', end - start)

import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(step, betti_curve_full)
ax2 = ax.twinx()
ax2.plot(step, eig_min)

'''