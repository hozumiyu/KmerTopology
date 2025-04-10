# -*- coding: utf-8 -*-
"""
Created on Tue Apr  8 15:44:16 2025

@author: yutah
"""

import argparse
from algorithm.auxilary import makeFolder, read_fasta

parser = argparse.ArgumentParser(prog='K-mer Topology',
                    description='Construct the Topological curves used for alignment-free sequence comparation')
parser.add_argument('--infile', default='CY117027.1.fasta',type=str, 
                    help='outfile will be stored in the folder test/ if not specified')
parser.add_argument('--outpath', default=None, type=str,
                    help='Where to store the outfile. Default will just be the infile name')

parser.add_argument('--mode', choices = ['homology', 'topology'], default='homology')
parser.add_argument('--kmers_size', default=2, type=int,
                    help='Which k to consider')
parser.add_argument('--step_size', default=None,type=int,
                    help='filtration size. Default is 4^{k-1}')
parser.add_argument('--max_step', default=50, type=int,
                    help='Maximum number of steps. Default is 50 steps')
args = parser.parse_args()


#create the outpath
if args.outpath == None:
    outpath = args.infile.replace(".fasta", "") ; makeFolder(outpath)
else:
    outpath = args.outpath; makeFolder(outpath)
    
    
sequence_id, sequence = read_fasta(args.infile)