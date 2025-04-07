# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 00:29:20 2025

@author: yutah
"""
import os
import numpy as np

'''
This file is used to save the Betti features, persistent diagrams, and loading the persistent diagram
Because the Persistent diagrams all start at 0, we do only save the death time to save space.
'''



def writeSpectrum(spectrum, kmers_list, outfile):
    '''
    Save the specta data of the kmers
    Ex:
        >AA
        0,0,1,2,3
    This would mean 2 connected components (Betti = 2), and the minimum eigenvalue is 1
    Each line corresponds to the filtration value you specified.
    Parameters
    ----------
    spectrum : dict
        spectra info
    kmers_list : list
        list of kmers
    outfile : str
        output file name

    Returns
    -------
    None.

    '''
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
    '''
    Load the spectrum file

    Parameters
    ----------
    infile : str
        input file name

    Returns
    -------
    spectrum : dict
        Dictionary of the spectrums. Keys are the kmers

    '''
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