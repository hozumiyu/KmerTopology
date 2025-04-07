# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 00:29:20 2025

@author: yutah
"""
from algorithm.kmers_homology import KmersHomology
import os
import numpy as np

'''
This file is used to save the Betti features, persistent diagrams, and loading the persistent diagram
Because the Persistent diagrams all start at 0, we do only save the death time to save space.
'''

def write_Betti(outfile, betti_dict, kmers_list):
    '''
    Save the betti numbers
    The betti_dict is a dictionary that contains all the betti numbers for each k-mer.
    Each line of the file is the betti number for the k-mers (sorted alphabetically)

    Parameters
    ----------
    outfile : str
        outfile name
    betti_dict : dict
        Dictionary containing the betti numbers.
    kmers_list : list
        list of the kmers. For 2 mers, [AA, AC, AG, AT ...].

    Returns
    -------
    None.

    '''
    file = open(outfile, 'w')
    for kmers in kmers_list:
        betti = betti_dict[kmers]
        s = ''
        for b0 in betti:
            s = s + str(b0) + ','
        s = s[:-1] + '\n'
        file.writelines(s)
    file.close()
    return

def obtain_betti_features(infile, step_size = 2, max_step = 50):
    '''
    For the Betti file that we generate from above, extract the Betti numbers we want.
    step_size and max_step should be tuned for your problem. Too small step size may not be meaningful. 
    For example, if you have 5mers, you have 1024 different 5 mers. The average occurace would only be 1 in 1024, so step size of 1 would not make sense.
    If you pick too large of a max_step, your tail end will all be 1.

    Parameters
    ----------
    infile : str
        the betti feature file we generate
    step_size : int, optional
        The step size for the filtration (default is 4^(k-1)). The default is 2.
    max_step : int, optional
        Number of step size. The default is 50.

    Returns
    -------
    betti : TYPE
        DESCRIPTION.

    '''
    file = open(infile)
    lines = file.readlines()
    file.close()
    betti = []
    for line in lines:
        line = line.strip()
        line = line.split(',')
        line = np.array(line).astype(float)
        len_line = len(line)
        feature_temp = np.ones(max_step)
        for idx in range(max_step):
            if idx*step_size < len_line:
                feature_temp[idx] = line[idx*step_size]
            else:
                feature_temp[idx] = line[-1]
        betti.append(feature_temp)
    betti = np.concatenate(betti)
    return betti

def write_pd0(outfile, accession, pd0_dict, kmers_list):
    '''
    Save the persistent diagram file
    We anly save the death time. Also, there is alot of repeated values. 
    Therefore, the format is the following
    >accession|pd0|kmers
    100:10
    20:5
    1:100
    This means 10 values have death time of 10, 5 has death time of 20, etc.
    If the kmers does not exist, it is stored as 0:1 (This means you should decrease your k because the choice of k is too big for the length of your sequence)
    
    Parameters
    ----------
    outfile : str
        outfile name.
    accession : str
        name of the sequence. By default, it will take on the name of the fasta file. >accession. This is to be consistent with some of the other software that deal with fasta file.
    pd0_dict : dict
        persistent diagram.
    kmers_list : list
        list of all the k-mers.

    Returns
    -------
    None.

    '''
    file = open(outfile, 'w')
    for kmers in kmers_list:
        file.writelines('>' + accession + '|pd0|' + kmers + '\n')
        pd0 = pd0_dict[kmers]
        if pd0.ndim == 1:
            if pd0[1] == np.inf:
                file.writelines('inf:1\n')
            else:
                file.writelines('0:1\n')
        else:
            unique_val = np.unique(pd0[:, 1])
            unique_val.sort()
            for unique in unique_val:
                count = np.where(pd0[:, 1] == unique)[0].shape[0]
                if unique != np.inf:
                    file.writelines('%d:%d\n'%(unique, count))
                else:
                    file.writelines('inf:%d\n'%count)
        
    file.close()
    return

def load_pd0(infile):
    '''
    Load the persistent diagram file

    Parameters
    ----------
    infile : str
        file name

    Returns
    -------
    pd0_dict : dict
        persistent diagram disctionary.

    '''
    file = open(infile)
    lines = file.readlines()
    file.close()
    pd0_dict = {}
    for idx, line in enumerate(lines):
        line=line.strip()
        if line[0] == '>':
            kmers = line.split('|')[-1]
            pd_kmers_death = []
        else:
            death, count = line.split(':')
            if death != 'inf':
                pd_kmers_death.append(np.repeat(float(death), count))
            else:
                pd_kmers_death.append(np.repeat(np.inf, count))
            
        if idx+1 == len(lines) or lines[idx+1][0] == '>':
            pd_kmers_death = np.concatenate(pd_kmers_death)
            pd_kmers = np.zeros([len(pd_kmers_death), 2])
            pd_kmers[:, 1] = pd_kmers_death
            pd0_dict[kmers] = pd_kmers
    return pd0_dict


def writeSpectrum(spectrum, kmers_list, outfile):
    file = open(outfile, 'w')
    for kmers in kmers_list:
        file.writelines('>%s'%kmers + '\n')
        eigen_list = spectrum[kmers]
        for eigen in eigen_list:
            s = ''
            for e in eigen:
                s = s + str(e) + ','
            s = s[:-1]
            file.writelines(s + '\n')
    file.close()
    return

def loadSpectrum(infile):
    file = open(infile)
    lines = file.readlines()
    file.close()
    numlines = len(lines)
    spectrum = {}
    start = True
    for idx, line in enumerate(lines):
        if line[0] == '>' and start == True:
            start = False
            kmers = line[1:].strip()
            eigen_val = []
        
        elif start == False:
            line = line.strip()
            line = line.split(',')
            line = np.array(line).astype(float)
            eigen_val.append(line)
            if idx+1 == numlines:
                eigen_val = np.array(eigen_val)
                spectrum[kmers] = eigen_val
            else:
                if lines[idx+1][0] == '>':
                    eigen_val = np.array(eigen_val)
                    spectrum[kmers] = eigen_val
                    start = True
    return spectrum