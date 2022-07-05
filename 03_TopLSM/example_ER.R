#' Simulated Example with Erdos-Renyi (ER) Models
#'
#' This script contains an end-to-end example from generating random graphs 
#' according to what is described in the paper from model fitting through 
#' inference. To run this script without hassle, it is strongly recommended 
#' to run this script on RStudio for setting up the path. 
#' 
#' We consider two classes of networks from ER models. In order to control 
#' simulation settings such as numbers of networks, link probabilities, and 
#' bounds for network sizes, change the values in 'parameters' section.

packs = c("knitr","igraph","network","TDAkit","T4cluster","latentnet",
          "rstudioapi","ggplot2","ggforce")
rpack = packs[!(packs %in% installed.packages()[,"Package"])]
if (length(rpack) > 0){
  install.packages(rpack)
}

library(knitr)
library(igraph)
library(network)
library(TDAkit)
library(T4cluster)
library(latentnet)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("src.R")


# parameters --------------------------------------------------------------
ngraph = 50    # number of generated graphs per class
prob_1 = 0.25  # link probability for class 1
prob_2 = 0.05  # link probability for class 2
n_min  = 80    # lower bound for network size
n_max  = 120   # upper bound for network size


# step 1. data generation -------------------------------------------------
# generate graphs from Erdos-Renyi models
list_all = list()
for (i in 1:ngraph){
  nsize1 = round(stats::runif(1, min=n_min, max=n_max))
  nsize2 = round(stats::runif(1, min=n_min, max=n_max))

  graph_er     = as.matrix(igraph::as_adjacency_matrix(igraph::sample_gnp(nsize1, prob_1)))
  graph_bi     = as.matrix(igraph::as_adjacency_matrix(igraph::sample_gnp(nsize2, prob_2)))
  
  list_all[[i]]        = network::as.network(graph_er, directed=FALSE)
  list_all[[i+ngraph]] = network::as.network(graph_bi, directed=FALSE)
}

# visualize two graphs from each class
x11()
par(mfrow=c(1,2))
plot(list_all[[1]],        main="graph from class 1")
plot(list_all[[2*ngraph]], main="graph from class 2")

# step 2. latent space embedding ------------------------------------------
n_all = length(list_all)
fit_latents = vector("list", n_all)
fit_effect  = rep(0, n_all)
fit_labels  = as.factor(rep(c("class1","class2"), each=ngraph))

for (i in 1:n_all){
  graph_i = list_all[[i]]
  fitts_i = latentnet::ergmm(graph_i~euclidean(d=2), tofit="mle")
  
  fit_effect[i]    = fitts_i$mle$beta
  fit_latents[[i]] = fitts_i$mle$Z
  print(paste0("step 2 : iteration ",i,"/",n_all," complete..."))
}


# step 3. compute persistent homology -------------------------------------
# compute persistence diagrams using Rips filtration
list_diagrams = list()
for (i in 1:length(list_all)){
  list_diagrams[[i]] = TDAkit::diagRips(fit_latents[[i]], maxdim=1)
}

# convert into persistence landscapes of orders 0 and 1
list_land0 = list()
list_land1 = list()
for (i in 1:length(list_diagrams)){
  PDi = list_diagrams[[i]]
  list_land0[[i]] = TDAkit::diag2landscape(PDi, dimension=0)
  list_land1[[i]] = TDAkit::diag2landscape(PDi, dimension=1)
  print(paste0("step 3 : iteration ",i,"/",n_all," complete..."))
}

# pairwise distance for future use
land0dist  = TDAkit::fsdist(list_land0)
land1dist  = TDAkit::fsdist(list_land1)

# step 4. network features ------------------------------------------------
network_class1 = list()
network_class2 = list()
for (i in 1:ngraph){
  A = list_all[[i]]
  B = list_all[[i+ngraph]]
  
  network_class1[[i]] = igraph::graph_from_adjacency_matrix(network::as.matrix.network.adjacency(A), mode="undirected")
  network_class2[[i]] = igraph::graph_from_adjacency_matrix(network::as.matrix.network.adjacency(B), mode="undirected")
}
features_class1 = cbind(fit_effect[1:ngraph], network_features(network_class1)) 
features_class2 = cbind(fit_effect[(ngraph+1):(2*ngraph)], network_features(network_class2)) 
features_print  = rbind(colMeans(features_class1), colMeans(features_class2))

# print on the console
colnames(features_print) = c("intercept","density","AD","ASD","diameter","CC")
rownames(features_print) = c("class 1", "class 2")
print(knitr::kable(features_print, caption = "Network summary statistics"))


# step 5. hypothesis testing ----------------------------------------------
# run the tests
table_SHT = array(0,c(2,2))
table_SHT[1,1] = TDAkit::fseqdist(list_land0, fit_labels, method = "original")$p.value
table_SHT[1,2] = TDAkit::fseqdist(list_land0, fit_labels, method = "disco")$p.value
table_SHT[2,1] = TDAkit::fseqdist(list_land1, fit_labels, method = "original")$p.value
table_SHT[2,2] = TDAkit::fseqdist(list_land1, fit_labels, method = "disco")$p.value

# print on the console
colnames(table_SHT) = c("k-sample","DISCO")
rownames(table_SHT) = c("order 0", "order 1")
print(knitr::kable(table_SHT, caption="Empirical p-values"))


# step 6. cluster analysis ------------------------------------------------
# run algorithms with k=2 (the number of clusters)
dim0_sc = TDAkit::fssc05Z(land0dist, k=2)
dim0_km = TDAkit::fskmedoids(land0dist, k=2)
dim0_kg = TDAkit::fskgroups(land0dist, k=2)
dim1_sc = TDAkit::fssc05Z(land1dist, k=2)
dim1_km = TDAkit::fskmedoids(land1dist, k=2)
dim1_kg = TDAkit::fskgroups(land1dist, k=2)

# compute clustering accuracy with Rand index
rand_idx = array(0,c(2,3))
rand_idx[1,1] = T4cluster::compare.rand(dim0_sc, fit_labels)
rand_idx[1,2] = T4cluster::compare.rand(dim0_km, fit_labels)
rand_idx[1,3] = T4cluster::compare.rand(dim0_kg, fit_labels)
rand_idx[2,1] = T4cluster::compare.rand(dim1_sc, fit_labels)
rand_idx[2,2] = T4cluster::compare.rand(dim1_km, fit_labels)
rand_idx[2,3] = T4cluster::compare.rand(dim1_kg, fit_labels)

# print on the console
colnames(rand_idx) = c("spectral clustering", "k-medoids","k-groups")
rownames(rand_idx) = c("order 0", "order 1")
print(knitr::kable(rand_idx, caption="Clustering accuracy for k=2"))

# step 7. metric MDS ------------------------------------------------------
embed_dim0 = TDAkit::fsmds(land0dist, method="metric")
embed_dim1 = TDAkit::fsmds(land1dist, method="metric")

# visualize
x11()
par(mfrow=c(1,2), pty="s")
plot(embed_dim0, xlab="MDS1", ylab="MDS2", col=fit_labels, 
     pch=19, main = "Distribution of PLs (order 0)")
plot(embed_dim1, xlab="MDS1", ylab="MDS2", col=fit_labels, 
     pch=19, main = "Distribution of PLs (order 1)")
