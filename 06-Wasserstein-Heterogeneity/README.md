# README

This repository contains `R` notebooks and compiled `md` documents for the paper

> ["On Heterogeneity in Wasserstein Space"](https://arxiv.org/abs/2603.14815)

written by Kisung You (CUNY Baruch College).

### Composition

This sub-repository includes two demonstration documents:

- [`demo-simulation.md`](demo-simulation.md), and
- [`demo-MNIST.md`](demo-MNIST.md).

Each file is compiled from a coresponding Quarto (`.qmd`) document within the `code` folder. Together, 
they illustrate the main components of the proposed framework:

- estimation of Wasserstein heterogeneity via pairwise distances,
- large-sample inference using empirical eccentricities,
- robustness through transform choice, and
- interpretability via identification of influential observations.

The simulation example provides a controlled setting with known structure, while the MNIST example demonstrates how the method operates on real image data.

### Notes

- If you want to reproduce the results from the source `.qmd` files, you will need to render the documents within the `code` folder.
- The MNIST example assumes precomputed pairwise distance matrices (`pdmat`) for each digit class, which are also included.
- Only within-class heterogeneity is illustrated in the notebook. Additional comparisons and extended analyses are omitted for brevity.
- The focus of these demos is clarity and reproducibility rather than computational optimization.

### Download the repo

```
git clone --depth 1 --filter=blob:none --no-checkout https://github.com/kisungyou/papers
cd papers/
git checkout master -- 06-Wasserstein-Heterogeneity
```