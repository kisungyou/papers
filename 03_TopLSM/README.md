# TopLSM

### Motivation
In network analysis, a class of [**latent space models**](https://sites.stat.washington.edu/raftery/Research/latent.html) is often used to learn latent representation of nodes given an instance of of a graph/network. This repo aims at demonstrating how population-level analysis given multiple latent space embeddings can benefit from tools of [topological data analysis (TDA)](https://en.wikipedia.org/wiki/Topological_data_analysis).




There are two files `example_ER.R` and `example_SBM.R` for replicating 
simulated examples. Both scripts contain installation of dependencies, 
path setting, data generation, inference, and visualization. The 'src.R' 
file contains auxiliary functions that are called from two example scripts. 
We strongly recommend to run these scripts in Rstudio.
