# K-mer Topopology

## Overview
K-mer topology is an alignment-free sequence comparison method that leverage tools from topological data analysis and spectral graph theory to uncover the shape of the genome. It is a k-mer based method, and utilizes the position of the k-mer to generate topological features. Currently, only DNA sequence is supported, but we will work on protein sequence in the future.

## Installation
K-mer topology is a python-based program, and requires Gudhi and Ripser. 

You can create a virutal environment and install the packages using the following:
```
conda create -n KmerTopology pip
pip install .
```




## Basic Usage
* If you only want the Betti Numbers, use KmerHomology
```
from KmerTopology.kmer_topology import KmerHomology
betti_features = KmerHomology(sequence, kmers_size, step_size, max_step)
```

* If you only want the Betti Numbers and the connectivity, using KmerTopology
```
from KmerTopology.kmer_topology import KmerTopology
betti_features, eigmin_features = KmerTopology(sequence, kmers_size, step_size, max_step)
```

* If you consolidate your features, you can compute the distance between sequences given a kmers size:
```
from KmerTopology.distance import single_scale_distance
D = single_scale_distance(X)
```

* If you consolidate your features, you can compute the multiscale distance distance between sequences given a kmers size:
```
from KmerTopology.distance import multiscale_distance
D = multiscale_distance(X, kmers_size)
```

* If you have a series of distance matrices, you can compute the weighted distance by:
```
from KmerTopology.distance import weighted_distance
D = weighted_distance(D_list, weight)
```


## Citation
Please cite us with

Hozumi, Yuta, and Guo-Wei Wei. "Revealing the Shape of Genome Space via K-mer Topology." arXiv preprint arXiv:2412.20202 (2024).

```
@article{hozumi2024revealing,
  title={Revealing the Shape of Genome Space via K-mer Topology},
  author={Hozumi, Yuta and Wei, Guo-Wei},
  journal={arXiv preprint arXiv:2412.20202},
  year={2024}
}
```



## Required Packages
The inidivudal packages that we utilized is listed below:
```
numpy
scikit-learn
scipy
gudhi
cypthon
ripser
biopython
matplotlib
```

## Respository File description
* KmerTopology/ - All the code
* example/ - contains a tutorial
* data/ - all the data used in the papar.