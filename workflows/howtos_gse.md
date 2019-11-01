---
title: 'How-to: perform gene set enrichment analysis'
---


# Set up environment 
```
library(EGAD)

gene_set_enrichment <- function(genes, genes.labels, voc){

        genes.names = rownames(genes.labels)
        labels.names = colnames(genes.labels)
        genes.counts = rowSums(genes.labels)
        labels.counts = colSums(genes.labels)                                   # p

        m = match ( genes, genes.names )
        filt.genes  = !is.na(m)
        filt.labels = m[filt.genes]


        labels.counts.set = rep( sum(filt.genes), length(labels.counts) )       # g

        m = match (labels.names, voc[,1])
        v.f = !is.na(m)
        v.g = m[v.f]

        universe = rep ( dim(genes.labels)[1], dim(genes.labels)[2])
        if(  length(filt.labels) == 1 ) { genes.counts.set = genes.labels[filt.labels,] }
        else { genes.counts.set = colSums(genes.labels[filt.labels,]) }            

        test =  cbind( (genes.counts.set -1) , labels.counts, universe-labels.counts, labels.counts.set)
        pvals = phyper(test[,1], test[,2], test[,3], test[,4], lower.tail=F)
        sigs = pvals < ( 0.05/length(pvals) )
        pvals.adj = p.adjust( pvals, method="BH")

        results = cbind(voc[v.g,1:2], test[v.f,c(1,2)], pvals[v.f], pvals.adj[v.f], sigs[v.f])
        colnames(results) = c("term", "descrp","p", "q", "pvals", "padj", "sig" )

        return (results)

}

# See also DE page for other functions
calc_DE <- function(cpm, filt.row, filt.col, group){
  X = cpm[filt.row,filt.col]
  group = group[filt.col]
  if( sum(group==1) < 2  ) {
    m.X1 = (X[,group==1])
  } else {
    m.X1 = rowMeans(X[,group==1])
  }
  if( sum(group==2) < 2  ) {
    m.X2 = (X[,group==2])
  } else {
    m.X2 = rowMeans(X[,group==2])
  }
  m.X = rowMeans(X)
  fc = log2(m.X1/m.X2)
  X.ps_g = sapply(1:dim(X)[1], function(k) wilcox.test(X[k,group==1], X[k,group==2], alt="g")$p.val)
  X.padj_g = p.adjust(X.ps_g , method = "BH")
  
  X.ps_l = sapply(1:dim(X)[1], function(k) wilcox.test(X[k,group==1], X[k,group==2], alt="l")$p.val)
  X.padj_l = p.adjust(X.ps_l , method = "BH")
  
  de = cbind(m.X, fc, X.ps_g, X.padj_g, X.ps_l, X.padj_l, m.X1, m.X2)
  return(de)
}

```
# Get data 
For gene set enrichment, you need a list of genes of interest (gene list), and a knowledge base (annotations) you want to assess it with. The function here also uses a vocabulary (to help understand the annotations), but is not necessary. Code can be tweaked to remove it. 
## DE gene list 
One type of gene list can be from a differential expression experiment. Here is an example. 
```
load("cpm.Rdata")

degs = calc_DE(cpm, filt.row, filt.col, groups)
filt.col = groups > 0
n = sum(filt.col)
filt.row = rowSums(cpm > 0) > (0.8*n)

fcs = deg[,2]
qval = deg[,4]
filt.sig = abs(fcs) >= 2 & qval <= 0.05 

gene_list = rownames(deg)[filt.sig]
```
## Other gene lists 
But of course, other gene lists can be used. 
```
gene_list = read.table("essential_genes") 
gene_list = read.table("hi_genes")
gene_list = read.table("synap_genes")
```
# Set up annotation matrix 
## GO annotations
These files/variables are loaded from EGAD or generated in this [tutorial](/workflows/howtos_go.md).
```
gogenes <- unique(GO.human[,1])
goterms <- unique(GO.human[,3])
annotations <- make_annotations(GO.human[,c(2,3)],gogenes,goterms)
voc <- GO.voc 
```
## Annotations from other databases
You can use any other type of data, as long as it can be tranformed into a binary matrix (indicating blah) in the format of genes by function/pathway/term. 
```
```

# Run
```
enrichments = gene_set_enrichment(gene_list, annotations, voc)
```



