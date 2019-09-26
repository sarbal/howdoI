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
```
# Get data 
```
gene_list = run_de(exprs, conditions)

gogenes <- unique(GO.human[,1])
goterms <- unique(GO.human[,3])
annotations <- make_annotations(GO.human[,c(2,3)],gogenes,goterms)
```

# Run
```
enrichments = gene_set_enrichment(gene_list, annotations, voc)
```



