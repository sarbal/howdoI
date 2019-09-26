---
title: 'How-to: run EGAD' 
---

EGAD is an R package for the functional analysis of gene networks. It contains a series of tools to calculate functional properties in networks based on the guilt-by-association principle. 

The functions implemented here can be applied to gene networks constructed from a range of data types (e.g., protein-protein interactions, expression, etc) across a subset of species with available functional annotations (e.g., human, mouse, zebrafish, worm, fly and yeast).

# Installation
You will need the latest version of [R](https://cran.r-project.org/). Then install via [bioconductor](http://bioconductor.org/). 
```{}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("EGAD")

library(EGAD)
```
# Set up your annotations and network  
## Protein-protein interaction network 
```{}
genelist <- make_genelist(biogrid)
gene_network <- make_gene_network(biogrid,genelist)
```

## Aggregate co-expression network 
```{}
netfile="blood.rerank.Rdata"
label="blood.rerank"
nettype= label
load(netfile)
gene_network = diag(length(genes.t))
bottom = row(gene_network) > col(gene_network)
colnames(gene_network) = genes.t
rownames(gene_network) = genes.t
gene_network[bottom] = temp
gene_network = gene_network + t(gene_network)
diag(gene_network) = 1
```

## GO annotations 
```{}
gogenes <- unique(GO.human[,2])
goterms <- unique(GO.human[,3])
annotations <- make_annotations(GO.human[,c(2,3)],gogenes,goterms)
```

# Running GBA 
```{}
aurocs_GO <- run_GBA(gene_network, annotations)
```

## Comparing AUROCs
```{}
aurocs_GO_nv = aurocs_GO[[1]][,1]
aurocs_GO_nd = aurocs_GO[[1]][,3]
plot_density_compare(auc_GO_nv, auc_GO_nd) 
```

