# -*- coding: utf-8 -*-
"""
Created on Sun Apr 27 16:04:37 2025

@author: yutah
"""

from sklearn.metrics import pairwise_distances
import numpy as np

def single_scale_distance(X):
    D = pairwise_distances(X)
    D = D / (np.mean(D) + 1e-6)
    return D

def multiscale_distance(X, kmers_size):
    max_step = int(X.shape[1] / (4**kmers_size))
    D = np.zeros([X.shape[0], X.shape[0]])
    for i in range(max_step):
        index = np.array( [i + max_step*k for k in range(int(4**(kmers_size)))])
        D_temp = pairwise_distances(X[:, index])
        D = D  + D_temp / np.mean(D_temp + 1e-6) 
    return D


def weighted_distance(D_list, weight):
    D = np.zeros_like(D_list[0])
    for i, DD in enumerate(D_list):
        D = D + weight[i] * DD
    return D