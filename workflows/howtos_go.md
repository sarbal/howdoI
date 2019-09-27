---
title: 'How-to: build the gene ontology'
---

# Get ontology
```
wget http://purl.obolibrary.org/obo/go.obo
```

# Get annotation files
```
wget http://geneontology.org/gene-associations/goa_human.gaf.gz
wget http://geneontology.org/gene-associations/gene_association.sgd.gz
wget http://geneontology.org/gene-associations/gene_association.mgi.gz
```

# Get gene info
```
wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz
```
## Make gene id 2 name file, for each species
```
zcat gene_info.gz | cut -f1-3 > temp
grep '^9606'   -w temp | cut -f2-3 > 9606.ID2name
grep '^10090'  -w temp | cut -f2-3 > 10090.ID2name
grep '^559292' -w temp | cut -f2-3 > 559292.ID2name
```

# Parse ontology (do only once)
```
perl ~/bin/ontology_parser/parse_GO_ontology.pl go.obo
perl ~/bin/ontology_parser/make_desc_list.pl go.obo.rel GOdesc
```

# Filter and run on human
```
zcat goa_human.gaf.gz | cut -f3,5,7,11  | grep -v ^! > goa_human.split
perl ~/bin/ontology_parser/parse_goa_with_iia.pl 9606.ID2name goa_human.split
perl ~/bin/ontology_parser/get_goa_with_desc_ext.pl goa_human.split.parsed GOdesc > goa_human.split.parsed.desc
```
 


