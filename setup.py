# -*- coding: utf-8 -*-
"""
Created on Sun Apr 27 16:11:55 2025

@author: yutah
"""
from setuptools import setup, find_packages
setup(
    name='KmerTopology',
    version='0.0',
    packages=find_packages(),
    install_requires=[
  	"numpy", 
 	"gudhi", 
  	"Cython",
  	"Ripser", 
  	"matplotlib",
  	"scikit-learn", 
  	"scipy",
        "biopython",
	"pandas"
    ],
    author='Yuta Hozumi',
    author_email='yhozumi3@gatech.edu',
    description='Kmer Topology extracts topological information about the distribution of kmers in a nucleotide sequence',
)
