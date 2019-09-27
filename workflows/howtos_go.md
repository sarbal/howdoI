---
title: 'How-to: build the gene ontology'
---

# Get ontology and annotation files 
```
wget http://purl.obolibrary.org/obo/go.obo
wget http://geneontology.org/gene-associations/goa_human.gaf.gz 
```
Other species annotations (ontology is the same) 
```
wget http://geneontology.org/gene-associations/gene_association.sgd.gz
wget http://geneontology.org/gene-associations/gene_association.mgi.gz
```

# Get gene info files to convert between IDs 
```
wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz
```
## Make gene id to name file, for each species
```
zcat gene_info.gz | cut -f1-3 > temp
grep '^9606'   -w temp | cut -f2-3 > 9606.ID2name  # Human
grep '^10090'  -w temp | cut -f2-3 > 10090.ID2name  # Mouse 
grep '^559292' -w temp | cut -f2-3 > 559292.ID2name # Yeast 
```

# Parse ontology (do only once)
```
perl /bin/parse_GO_ontology.pl go.obo
perl /bin/make_desc_list.pl go.obo.rel GOdesc
```

# Filter and build on human data 
```
zcat goa_human.gaf.gz | cut -f3,5,7,11  | grep -v ^! > goa_human.split
perl /bin/parse_goa_with_iia.pl 9606.ID2name goa_human.split
perl /bin/get_goa_with_desc_ext.pl goa_human.split.parsed GOdesc > goa_human.split.parsed.desc
```

# Parse and process in R 
```
make_network_from_data <- function(data, listA, listB){
        nr = length(listA)
        nc = length(listB)

        net = matrix(0, ncol=nc, nrow=nr)

        m = match( (data[,1]), listA  )
        p1 = !is.na(m)
        m1 = m[p1]

        m = match( (data[p1,2]), listB )
        p2 = !is.na(m)
        m2 = m[p2]

        net[cbind(m1,m2)] = 1

        rownames(net) = listA
        colnames(net) = listB
        return(net)
}

GO.human = read.table("goa_human.split.parsed.desc", sep="\t")
filt = GO.human[,4] != "IEA"
GO.human.nonIEA = make_network_from_data( GO.human[filt,2:3], unique(GO.human[,2]), unique(GO.human[,3]) )
genes = unique(GO.human[,2])
goterms = unique(GO.human[,3])
voc =  read.table("go.obo.voc", sep="\t", quote="")
save( GO.human.nonIEA, GO.human, genes, goterms,voc, make_network_from_data, file="GO.human.Rdata")
```




