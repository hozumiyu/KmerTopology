# -*- coding: utf-8 -*-
"""
Created on Sun Apr 27 15:03:58 2025

@author: yutah
"""

import numpy as np
from scipy import sparse
from gudhi.representations.vector_methods import BettiCurve
import ripser

def compute_filtration_topology(positions, step_size, max_step):
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
    positions = positions[positions> 0]  #only keep the nonzero entry
    D = np.zeros([positions.shape[0], positions.shape[0]])
    for idx in range(positions.shape[0]):
        for idx2 in range(idx+1, positions.shape[0]):
            D[idx, idx2] = positions[idx2] - positions[idx]
            D[idx2, idx] = positions[idx2] - positions[idx]
    betti = []
    eig_min = []
    filtration = np.linspace(0, step_size*max_step, max_step)
    for idx, s in enumerate(filtration):
        if s == 0:
            eig = np.zeros(D.shape[0])
        else:
            A = D.copy()
            A[A>s] = 0
            A[A!=0] = 1
            L = np.diag(np.sum(A,axis = 0)) - A
            eig = np.linalg.eigvalsh(L)
        numzero = np.where(eig <1e-6)[0].shape[0]
        betti.append(numzero)
        if numzero == len(eig):
            eig_min.append(0)
        else:
            eig_min.append(eig[numzero])
    return np.array(betti), np.array(eig_min)



def compute_kmers_persistent_diagram(positions, step_size, max_step):
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
        
    
    positions = positions[positions> 0]  #only keep the nonzero entry
    
    #dis = sparse.coo_array((pos.shape[0], pos.shape[0]))
    row = []; col = []; val = []
    if positions.shape[0] > 1:
        for i in range(positions.shape[0]): #compute the neighboring positions' distance
            if i ==0:
                row.append(i); col.append(i+1); val.append(positions[i+1] - positions[i])
            elif i == positions.shape[0] - 1:
                row.append(i); col.append(i-1); val.append(positions[i] - positions[i-1])
            else:
                row.append(i); col.append(i-1); val.append(positions[i] - positions[i-1])
                row.append(i); col.append(i+1); val.append(positions[i+1] - positions[i])
        row = np.array(row); col = np.array(col); val = np.array(val)
        dis = sparse.coo_array( (val, (row, col)), shape = (positions.shape[0], positions.shape[0]))
        pd_kmers = ripser.ripser(dis, thresh = step_size*max_step, distance_matrix = True)['dgms'][0]
    elif positions.shape[0] == 1:
        pd_kmers = np.array([[0, np.inf]])
    elif positions.shape[0] == 0:
        pd_kmers = np.array([[0, 0]])
    return pd_kmers



    
def compute_kmers_betti(pd_kmers, step_size, max_step):
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
            betti_curve_kmers = np.array([0,0])
        else:
            betti_curve_kmers = np.array([1,1])
    else:
        sorted_kmers = pd_kmers[:, 1].copy()
        sorted_kmers[sorted_kmers == np.inf] = 0
        filtration = np.linspace(0, step_size*max_step, max_step)
        
        bc = BettiCurve(predefined_grid = filtration)
        betti_curve_kmers = bc.fit_transform([pd_kmers]).reshape(-1)
    return betti_curve_kmers