# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 00:32:16 2025

@author: yutah
"""

import argparse, os, time
from Bio import SeqIO
from algorithm.kmers_homology import KmersHomology
from algorithm.auxilary import obtain_betti_features, write_Betti, load_pd0, write_pd0
import numpy as np

def makeFolder(outpath):
    try:
        os.makedirs(outpath)
    except:
        return
    return

parser = argparse.ArgumentParser(prog='Monomer_homology',
                    description='Construct betti curve for the homology')
parser.add_argument('--infile', default='test.fasta',type=str)
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
    
outpath = args.name + '_homology/' ; makeFolder(outpath)
pd0_file = outpath + args.name + '_%d_pd0.txt'%args.kmers_size
betti0_file = outpath + args.name + '_%d_betti0.txt'%args.kmers_size
feature_file = outpath + args.name + '_%d_betti0_feature.npy'%args.kmers_size


start = time.time()

myKmersHomology = KmersHomology(kmers_size = args.kmers_size, thresh = -1)
myKmersHomology.find_kmers_position(sequence)
kmers_list = myKmersHomology.kmers_list
kmers_list.sort()

if not os.path.exists(pd0_file): 
    print(pd0_file, 'does not exist. Compute homology features')
    
    print('Beginning homology computation for ', accession)

    betti_dict = {}
    pd0_dict = {}
    for kmers in kmers_list:
        print(accession, kmers)
        pd_kmers = myKmersHomology.compute_kmers_persistent_diagram(kmers)
        pd0_dict[kmers] = pd_kmers
    
        betti_dict[kmers] = myKmersHomology.compute_kmers_betti(pd_kmers)
    write_pd0(pd0_file, accession, pd0_dict, kmers_list)
    write_Betti(betti0_file, betti_dict, kmers_list)
    
else:
    pd0_dict = load_pd0(pd0_file)
    betti_dict = {}
    for kmers in kmers_list:
        betti_dict[kmers] = myKmersHomology.compute_kmers_betti(pd0_dict[kmers])
    write_Betti(betti0_file, betti_dict, kmers_list)
print('Finished computing all homology for ', accession)
outfile_betti = '%s_%dmers_betti0.txt'%(accession, args.kmers_size)

end = time.time()
print('Time: ', end - start)


if args.step_size == None:
    args.step_size = int(4**(args.kmers_size - 1))
if args.max_step == None:
    args.max_step = 50
betti = obtain_betti_features(betti0_file, step_size = args.step_size, max_step = args.max_step )

np.save(feature_file, betti)