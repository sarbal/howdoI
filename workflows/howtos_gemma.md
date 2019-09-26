---
title: 'How-to: download data from gemma'
---

https://www.ncbi.nlm.nih.gov/pubmed/22782548
https://gemma.msl.ubc.ca/home.html
https://github.com/PavlidisLab/gemmaAPI.R

# Install the wrapper 
```
library(devtools)
devtools::install_github('PavlidisLab/gemmaAPI.R')
```

# Download data
```
data = datasetInfo('GSE19804', request='data', filter = TRUE, return = TRUE, file = NULL )

```
# Get metadata
```
meta = compileMetadata('GSE19804', outputType = "list")  
```

