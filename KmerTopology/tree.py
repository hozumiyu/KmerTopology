# -*- coding: utf-8 -*-
"""
Created on Sun Apr 27 16:16:51 2025

@author: yutah
"""
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor, DistanceMatrix
import matplotlib.pyplot as plt

def UPGMA(D, labels, outfile):
    DD = []
    for idx in range(D.shape[0]):
        DD.append(list(D[idx, :idx+1]))
        
    
    
    
    #$color_dict = { unique_group[i]: color_list[i] for i in range(len(unique_group))}
    
    distMatrix = DistanceMatrix(names = labels, matrix = DD)
    constructor = DistanceTreeConstructor()
    UPGMATree = constructor.upgma(distMatrix)
    
    
    
        
    for idx, clade in enumerate(UPGMATree.find_clades()):
        if clade.name not in labels:
            clade.name = ''
    fig = plt.figure()   
    axes = fig.add_subplot(111)
    Phylo.draw(UPGMATree, axes = axes, show_confidence=True, do_show=False)
    axes.set_xticks([])
    axes.set_yticks([])
    axes.set_xlabel('')
    axes.set_ylabel('')
    
    Phylo.write(UPGMATree, outfile,"newick")
    return fig, axes
        
        
        