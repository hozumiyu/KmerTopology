# K-mer Topopology

## Overview
K-mer topology is an alignment-free sequence comparison method that leverage tools from topological data analysis and spectral graph theory to uncover the shape of the genome. It is a k-mer based method, and utilizes the position of the k-mer to generate topological features. Currently, only DNA sequence is supported, but we will work on protein sequence in the future.

## Dependency
K-mer topology is a python-based program, and requires Gudhi and Ripser. To run our tutorial and code, you will also need biopython (loading the fasta seq), numpy, panda, matplotlib, seaborn (and dependencies on these packages)

You can create a virutal environment and install the packages using the following:
```
conda create -n KmerTopology
pip install -r requirement.txt
```

The inidivudal packages that we utilized is listed below:
```
numpy=1.26.4
pandas=2.2.2
scikit-learn=1.5.1
scipy=1.14.0
gudhi=3.10.1
ripser=0.6.10
biopython=1.8.4
matplotlib=3.9.1
seaborn=0.13.2
```

## Tutorial and Usage
We have a detailed tutorial in ``tutorial.jpy``. Please refer to the example to use the code.

In the current version, the code can only take 1 input sequence. If you have a fasta file with more than 1 sequence, you will need to split it into individual file, or modify the code.



## Respository File description
