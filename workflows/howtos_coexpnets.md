---
title: 'How-to: build co-expression networks'
---
 
Gene expression profiles that co-vary suggest co-regulation, co-function or co-localization. This pipeline involves the construction of individual gene co-expression networks, the aggregation of these networks, and then an assessment via neighbor-voting. Differences in bulk and single-cell expression data mean care needs to be taken when generating these networks. We provide guidelines for both. 

# Co-expression - bulk RNA-seq
## Data 
Read in expression data. This could be counts or normalized data. A minimum of 20 samples per experiment is recommended. And genes with expression in at least 80% of these samples, with at least 1 CPM.  
```{}
counts = read.table("counts.txt")
cpm = calc_cpm(counts)
nsamp = dim(counts)[2] 
keep = rowSums(cpm >=1) > (0.8 * nsamp) 
```

##  Building coexpression networks 
```{}
net = EGAD::build_coexp_network(cpm, genes[keep])
nd = EGAD::node_degree(net)
EGAD::plot_distribution(nd, xlab="Node degrees")
```
Or 
```{}
gene.corr = cor( t(cpm), method="s")) 
bottom = row(gene.corr) > col(gene.corr) 
gene.ranked = rank(gene.corr[bottom] , na.last = "keep", ties.method = "average")   
gene.ranked = gene.ranked/max(gene.ranked) 
net = gene.corr * 0 
net[bottom] = gene.ranked
net = net + t(net) 
diag(net) = 1
```


##  Aggregating 
```{} 
agg = diag(sum(keep))
 
for (i in 1:n_experiments) {
  id = experiment_ids[i]
  f.c = which(samples[,1] == id)

  sub = exprs[keep,f.c]
  net = build_coexp_network(sub, genes[keep])
  med = median(net, na.rm=T)
  net[is.na(net)] = med
      
  # Save individual  networks 
  # save(net, file=paste("coexp",i, "Rdata", sep=".") )
  if(i==1 ) {
    agg = net
  } else {
    agg = net + agg
  }
}

agg.rank =  matrix(rank(agg, na.last = "keep", ties.method = "average"), nrow=dim(agg)[1], ncol=dim(agg)[2] )
rownames(agg.rank) = rownames(agg)
colnames(agg.rank) = rownames(agg)
agg.rank = agg.rank/max(agg.rank, na.rm=T)
#save(agg.rank, file="coexp.agg.rank.Rdata")

```


##  Assessing  
```{}
aurocs = run_GBA(agg.rank, annotations)
aurocs[[2]] = ""
EGAD::plot_density_compare(aurocs[[1]][,1], aurocs[[1]][,3])
```

 
#  Co-expression - single-cell RNA-seq
## Data 
```{}
library(Seurat)
data <- Read10X(data.dir = "10x/v3_chemistry/10k_pbmcs_healthy_donor/filtered_feature_bc_matrix/")
pbmc <- CreateSeuratObject(counts = data, project = "pbmc", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
pbmc <- subset(pbmc, subset = nFeature_RNA > 200   & percent.mt < 20)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
DimPlot(pbmc, reduction = "pca")
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunTSNE(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "tsne")
save(pbmc, file="pbmc.Rdata")
```


##  Building coexpression networks 
```{}
exprs = pbmc@assays$RNA@scale.data
net = build_coexp_network(exprs, genes )
save(net, file="coexp.scaled.Rdata")

exprs = as.matrix(pbmc@assays$RNA@data)
net = build_coexp_network(exprs, genes )
save(net, file="coexp.norm.Rdata")

exprs = as.matrix(pbmc@assays$RNA@counts)
net = build_coexp_network(exprs, genes )
save(net, file="coexp.raw.Rdata")


exprs = calc_cpm( as.matrix(pbmc@assays$RNA@counts))
net = build_coexp_network(exprs, genes )
save(net, file="coexp.cpm.Rdata")

```

##  Aggregating 
```{} 
agg = diag(sum(keep))
 
for (i in 1:n_experiments) {
  id = experiment_ids[i]
  f.c = which(samples[,1] == id)

  sub = exprs[keep,f.c]
  net = build_coexp_network(sub, genes[keep])
  med = median(net, na.rm=T)
  net[is.na(net)] = med
      
  # Save individual  networks 
  # save(net, file=paste("coexp",i, "Rdata", sep=".") )
  if(i==1 ) {
    agg = net
  } else {
    agg = net + agg
  }
}

agg.rank =  matrix(rank(agg, na.last = "keep", ties.method = "average"), nrow=dim(agg)[1], ncol=dim(agg)[2] )
rownames(agg.rank) = rownames(agg)
colnames(agg.rank) = rownames(agg)
agg.rank = agg.rank/max(agg.rank, na.rm=T)
#save(agg.rank, file="coexp.agg.rank.Rdata")

```


## Assessing  
```{}
aurocs = run_GBA(agg.rank, annotations)
aurocs[[2]] = ""
EGAD::plot_density_compare(aurocs[[1]][,1], aurocs[[1]][,3])
```

 
