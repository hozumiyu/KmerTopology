# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 00:32:16 2025

@author: yutah
"""

import argparse, os, time
from Bio import SeqIO
from algorithm.kmersCircularLaplacian import KmersCircularLaplacian
from algorithm.auxilary import writeSpectrum, loadSpectrum
import numpy as np

def makeFolder(outpath):
    try:
        os.makedirs(outpath)
    except:
        return
    return

parser = argparse.ArgumentParser(prog='Monomer_homology',
                    description='Construct betti curve for the homology')
parser.add_argument('--infile', default='OQ852526.1_S.fasta',type=str)
parser.add_argument('--name', default=None ,type=str)
parser.add_argument('--kmers_size', default=2,type=int)
parser.add_argument('--step_size', default=None,type=int)
parser.add_argument('--max_step', default=None,type=int)
args = parser.parse_args()


all_sequences = SeqIO.parse(args.infile,  'fasta')
for fasta in all_sequences:
    accession, sequence = fasta.id, str(fasta.seq)
    
if args.name == None:
    args.name = accession

if args.step_size == None:
    args.step_size = int(4**(args.kmers_size - 1))
if args.max_step == None:
    args.max_step = 50
steps = [args.step_size * i for i in range(args.max_step)]


outpath = args.name + '_circular_Laplacian/' ; makeFolder(outpath)
spectrum_file = outpath + args.name + '_%d_circular_spectrum.txt'%args.kmers_size
betti0_file = outpath + args.name + '_%d_circular_betti0_feature.npy'%args.kmers_size
eigmin_file = outpath + args.name + '_%d_circular_eigmin_feature.npy'%args.kmers_size

start = time.time()

myKmerLap = KmersCircularLaplacian(kmers_size = args.kmers_size, steps = steps)
myKmerLap.find_kmers_position(sequence)
kmers_list = myKmerLap.kmers_list
kmers_list.sort()

if not os.path.exists(spectrum_file): 
    print(spectrum_file, 'does not exist. Compute homology features')
    
    print('Beginning homology computation for ', accession)
    spectrum = {}
    eigen_minimum = {}
    betti_dict = {}
    for kmers in kmers_list:
        print(accession, kmers)
        D = myKmerLap.compute_kmers_Distance(kmers)
        eigenval_list = myKmerLap.computeLaplacianEig(D)
        b, e = myKmerLap.compute_kmers_betti_min(eigenval_list)
        spectrum[kmers] = eigenval_list
        betti_dict[kmers] = b
        eigen_minimum[kmers] = e
    writeSpectrum(spectrum, kmers_list, spectrum_file)
        
    
else:
    spectrum = loadSpectrum(spectrum_file)
    betti_dict = {}
    eigen_minimum = {}
    for kmers in kmers_list:
        b, e = myKmerLap.compute_kmers_betti_min(spectrum[kmers])
        betti_dict[kmers] = b
        eigen_minimum[kmers] = e

betti_feature = [betti_dict[kmers] for kmers in kmers_list ]; betti_feature = np.concatenate(betti_feature)
eigen_feature = [eigen_minimum[kmers] for kmers in kmers_list ]; eigen_feature = np.concatenate(eigen_feature)

np.save(betti0_file, betti_feature)
np.save(eigmin_file, eigen_feature)