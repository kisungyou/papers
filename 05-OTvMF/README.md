# README

This repo contains `R` notebooks and compiled `html` documents for the following paper:

> ["Learning over von Mises--Fisher distributions via a Wasserstein-like-geometry"](
https://doi.org/10.48550/arXiv.2504.14164)

written jointly by Kisung You (CUNY Baruch College), Dennis Shung (Mayo Clinic), and Mauro Giuffrè (Yale University, University of Trieste). 


### Composition

There are 3 html files 

- `demo-interpolation.html`, 
- `sim-4-types.html`, and 
- `sim-reduce.html`

in the top level, which are compiled quarto documents in the style of typical notebooks. Each html file 
was generated from running `qmd` file under the corresponding folder of the same name. We note that 
real data examples are excluded for brevity and polyglot workflow that can't be easily handled 
in a single quarto file without making the code cryptic. Still, these three notebooks use 
all major functionalities listed in the table below.

### Execution

If you want to run `qmd` files, you will need to have `R` installed along with the **[maotai](https://cran.r-project.org/web/packages/maotai/index.html)** package with version `0.3.0` or higher available on CRAN. The package contains following functions 
that implement core routines proposed in the paper:

| Function        | Description |
| :- | :---        |
| `WLpdist`       | Pairwise Wasserstein-like Distance between two vMF distributions |
| `WLbarycenter`  | Barycenter of vMF Distributions Under a Wasserstein-Like Geometry |
| `WLmedian`      | Geometric Median of vMF Distributions Under a Wasserstein-Like Geometry |
| `movMF_convert` | Convert `movMF` object |
| `movMF_info`    | Extract meaningful information from the vMF mixture model |
| `movMF_reduce_greedy` | vMF mixture model reduction - Greedy method | 
| `movMF_reduce_partitional` | vMF mixture model reduction - partitional method | 

### Download the repo

If you want to download this repo only, please use the following commands in your CLI.
```
git clone --depth 1 --filter=blob:none --no-checkout https://github.com/kisungyou/papers
cd papers/
git checkout master -- 05-OTvMF
```



