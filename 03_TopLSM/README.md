# Comparing multiple latent space embeddings using topological analysis

### 1. Motivation
In network analysis, a class of [**latent space models**](https://sites.stat.washington.edu/raftery/Research/latent.html) is often used to learn latent representation of nodes given an instance of of a graph/network. This repo aims at demonstrating how population-level analysis given multiple latent space embeddings can benefit from tools of [topological data analysis (TDA)](https://en.wikipedia.org/wiki/Topological_data_analysis).

### 2. Download

In order to download this folder only, please use the following commands.
```
git clone --depth 1 --filter=blob:none --no-checkout https://github.com/kisungyou/papers
cd papers/
git checkout master -- 03_TopLSM
```

### 3. Workflow

There are two main scripts - `example_ER.R` and `example_SBM.R` - that replicate 
examples from the manuscript. Both deliver end-to-end analysis, including

- generating random graphs by [Erdős–Rényi model][1] and [Stochastic block model][2],
- fitting the standard LSM model to the generated graphs, 
- learning topological representations from estimated embeddings, and 
- performing hypothesis testing and clustering.

### 4. Notes

* In order to control simulation settings such as numbers of networks, link 
probabilities, and bounds for network sizes, change the values in `parameters` section.
* The `src.R` file contains auxiliary functions that are called from two example scripts.
* We strongly recommend to run these scripts in **Rstudio** for path setting.


[2]: https://en.wikipedia.org/wiki/Stochastic_block_model
[1]: https://en.wikipedia.org/wiki/Erd%C5%91s%E2%80%93R%C3%A9nyi_model
