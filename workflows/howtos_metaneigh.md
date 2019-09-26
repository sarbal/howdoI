---
title: 'How-to: test for cell type replicability using MetaNeighbor'
---

#  Install MetaNeighbor 
In R: 
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("MetaNeighbor")
library(MetaNeighbor) 
```

# Load data 
```
load("/data/XCGD/SCGC-GILL-JG-03_2/mono_seo.Rdata")
seo1 = seo
load("/data/XCGD/SCGC-GILL-JG-04_2/mono_seo.Rdata")
seo2 = seo
```

# Merge data sets 
```
colData1 = seo1@colData
colData2 = seo2@colData

colData = rbind(colData1, colData2)
data1 = seo1@assays[[1]]
data2 = seo2@assays[[1]]


m = match( rownames(data1), rownames(data2))
f.d1 = !is.na(m)
f.d2 = m[f.d1]

data1 = data1[f.d1,]
data2 = data2[f.d2,]

cd1 = paste0(colnames(data1), "_1" )
cd2 = paste0(colnames(data2), "_2" )


rownames(colData) = c(cd1,cd2)
colnames(colData)[5] = "cell_type"

study_id = c(rep(1, length(cd1)), rep(2, length(cd2))  )


names(study_id) = rownames(colData)
colData = cbind(colData, study_id)
sample_id = rownames(colData)
names(sample_id) = rownames(colData)
colData = cbind(colData, sample_id)

data = cbind(data1, data2)


colnames(data) = names(sample_id)
data = as.matrix(data)
seo = SummarizedExperiment(assays=list(counts=data), colData=colData)
```


# Run analysis 
```
var_genes = variableGenes(dat = seo, exp_labels = seo$study_id)
celltype_NV = MetaNeighborUS(var_genes = var_genes,
                             dat = seo,
                             study_id = seo$study_id,
                             cell_type = seo$cell_type)
                             
top_hits = topHits(cell_NV = celltype_NV,
                   dat = seo,
                   study_id = seo$study_id,
                   cell_type = seo$cell_type,
                   threshold = 0.9)
                   
cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)
gplots::heatmap.2(celltype_NV,
                  margins=c(8,8),
                  keysize=1,
                  key.xlab="AUROC",
                  key.title="NULL",
                  trace = "none",
                  density.info = "none",
                  col = cols,
                  breaks = breaks,
                  offsetRow=0.1,
                  offsetCol=0.1,
                  cexRow = 2,
                  cexCol = 2,
                  ColSideCol = col.samp[samples[,1]],
                  RowSideCol = col.samp2[samples2[,4]]
                  )

```

