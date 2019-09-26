---
title: 'How-to: run a differential expression analysis'
---
# Running DESeq2 
For comprehensive tutorials, see here: 
# Running edgeR
# Running basic DE 

```
calc_DE <- function(X, f.a, filt, group){
	X = X[f.a,filt]
	group = group[filt]
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

	X.ps = sapply(1:dim(X)[1], function(k) wilcox.test(X[k,group==1], X[k,group==2])$p.val)
	X.padj = p.adjust(X.ps , method = "BH")

	de = cbind(m.X, fc, X.ps_g, X.padj_g, X.ps_l, X.padj_l, X.ps, X.padj, m.X1, m.X2)
	return(de)
}

deg = calc_DE(cpm, keep, subset, groups  )

```




## Setting up
## DE between cases and controls 
## DE between multiple conditions 
## DE with low sample numbers 
## DE with unbalanced samples 




