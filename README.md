
<!-- README.md is generated from README.Rmd. Please edit that file -->
XGSA: a statistical method for cross-species gene set analysis
==============================================================

Introduction
------------

XGSA is an R package that facilitates cross(X)-species Gene Set Analysis as described in - XGSA: a statistical method for cross species gene set analysis, 2016, Djordjevic D, Ho JWK, Kusumi K, Bioinformatics, DOI: 10.1093/bioinformatics/btw428

XGSA accesses Ensembl through the Biomart portal to automatically deal with homology mapping between species. Because of this we use Ensembl IDs to represent our genes. Ensembl / Biomart occassionally goes down. This can result in strange errors such as "Extra content at the end of the document". You can check the status of the Ensembl web services here <https://wtsi-status.blogspot.com.au/> and here <http://www.ensembl.info/blog/category/web/web-status/> . Currently the only solution is to wait for Ensembl to come fully back online and try again.

XGSA was written by Djordje Djordjevic - <djordje.evic@gmail.com>

Installation
------------

Open an R session.

XGSA depends on 2 packages, slam and biomaRt. Constructing the Gene Ontology relies on another 4 packages AnnotationDbi, GO.db, graph and igraph. Make sure these are available to you.

``` r
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
biocLite("AnnotationDbi")
biocLite("GO.db")
biocLite("graph")


install.packages(c("slam", "igraph"))
```

Install XGSA using:

Make sure you have devtools installed, then install XGSA from github:

``` r
install.packages('devtools')
devtools::install_github('VCCRI/XGSA')
```

Then, load the package using:

``` r
library(xgsa)
```

That's it! You can run the entire following example using:

``` r
example(xgsa)
```

Example
-------

### Load experimental data

First let's load some data. Remember we need to end up with Ensembl gene IDs. In this example we will extract a list of mouse cardiac development genes from a previous published study.

``` r
# In this example we will extract a list of mouse cardiac development genes from a previous published study.
cardiac.perturbation.data <- read.table("http://cardiaccode.victorchang.edu.au/data/Cardiaccode_in_vivo_evidence_2014_06_10.txt", sep="\t", header=TRUE, quote="\"")
mouse.cardiac.genes <- unique(unlist(cardiac.perturbation.data[cardiac.perturbation.data$Species == "Mus musculus", c("Regulator","Target")]))
#
```

XGSA harnesses the Ensembl homology pipeline, and so we need to convert all of our data sets to Ensembl gene IDs We will convert the mouse gene symbols to Ensembl IDs using the XGSA helper function "get\_ENSEMBL\_symbol\_map".

``` r
# XGSA harnesses the Ensembl homology pipeline, and so we need to convert all of our data sets to Ensembl gene IDs
# We will convert the mouse gene symbols to Ensembl IDs using the XGSA helper function "get_ENSEMBL_symbol_map".
mouse.ensembl.symbol.map <- get_ENSEMBL_symbol_map(species = 'mmusculus')
#> Loading required package: biomaRt
mouse.cardiac.ensembl.symbols <- mouse.ensembl.symbol.map$ensembl_gene_id[mouse.ensembl.symbol.map$external_gene_name %in% mouse.cardiac.genes]
#
```

Now we have a list of Ensembl IDs for mouse cardiac genes, we will turn it into an XGSA data set. Note that the input data MUST be a named list, which allows for multiple gene sets in the same data set. As we don't have a defined gene universe for this data set we will use all mouse Ensembl IDs that have an external gene symbol.

``` r
# Now we have a list of Ensembl IDs for mouse cardiac genes, we will turn it into an XGSA data set.
# Note that the input data MUST be a named list, which allows for multiple gene sets in the same data set.
# As we don't have a defined gene universe for this data set we will use all mouse Ensembl IDs that have an external gene symbol.
mouse.data <- new_XGSA_dataset(species = 'mmusculus', data = list(mouseCardiacGenes = mouse.cardiac.ensembl.symbols), type = 'genesetlist', name = 'MouseCardiacGenes', universe = unique(mouse.ensembl.symbol.map$ensembl_gene_id))
#
```

### Load reference dataset (i.e. Gene Ontology)

In this example we will compare to the zebrafish Gene Ontology using "direct" evidence only - this means the annotations are NOT transferred between species. We will use another XGSA helper function "get\_GO" to retrieve the latest Gene Ontology information from Ensembl. The gene universe we will use is all of the zebrafish biological process genes that we are testing.

``` r
# In this example we will compare to the zebrafish Gene Ontology using "direct" evidence only - this means the annotations are NOT transferred between species.
# We will use another XGSA helper function to retrieve the latest Gene Ontology information from Ensembl "get_GO_list_from_ontologies_with_evidence_codes".
# The gene universe we will use is all ofthe zebrafish biological process genes that we are testing. 
zebrafish.GO <- get_GO('drerio', ontologies = "biological_process")
#> [1] "retrieved GO"
zebrafish.GO <- zebrafish.GO[lapply(zebrafish.GO, length) > 10 & lapply(zebrafish.GO, length) < 500]
zebrafish.GO.data <- new_XGSA_dataset(species = "drerio", data = zebrafish.GO, type = 'genesetlist', name = "ZebrafishGO", universe = unique(unlist(zebrafish.GO)))
#
```

### Run the XGSA test!

Now we can compare the mouse cardiac genes to the zebrafish gene ontology.

``` r
# Now we can compare the mouse cardiac genes to the zebrafish gene ontology.
mouse.cardiac.vs.zebrafish.GO.results <- run_XGSA_test(mouse.data, zebrafish.GO.data)
#
```

### Examining the results

We need to separate the pvalues and the overlapping gene IDs, because XGSA returns both. The p-values are stored in the first element of each result, and the overlapping genes are stored in the second element.

``` r
# We need to separate the pvalues and the overlapping gene IDs, because XGSA returns both.
# The p.values are stored in the first element of each result, and the overlapping genes are stored in the second element.
resulting.pvals <- lapply(mouse.cardiac.vs.zebrafish.GO.results, function(X){ X[["pvals"]] })
resulting.overlap.genes <- lapply(mouse.cardiac.vs.zebrafish.GO.results, function(X){ X[["genes"]] })
#
```

Now we perform Benjamini-Hochberg multiple hypothesis testing correction to the p-values.

``` r
# Now we perform Benjamini Hochberg multiple hypothesis testing correction to the pvalues.
adjusted.pvals <- p.adjust(unlist(resulting.pvals), method = "BH")
#
```

We need to make the names of our results interpretable for humans, so we extract the GO Term IDs We can use another XGSA helper function to find out the GO term names, and finally we get interpretable names

``` r
# We need to make the names of our results interpretable for humans, so we extract the GO Term IDs
names(adjusted.pvals) <- unlist(lapply(strsplit(names(adjusted.pvals) ,"\\."), function(X){return(X[[2]])}))
# We can use another XGSA helper function to find out the GO term names.
zebrafish.GO.names <- get_GO_names('drerio')
# And finally we get interpretable names
names(adjusted.pvals) <- zebrafish.GO.names[match( names(adjusted.pvals), zebrafish.GO.names$go_id),"name_1006"]
#
```

Now let's look at the 10 most significant GO term results

``` r
# Now let's look at the 10 most significant GO term results
significant.GO.Terms <- adjusted.pvals[which(adjusted.pvals < 0.05)]
print(head(sort(significant.GO.Terms),10))
#>                             cell fate commitment 
#>                                     8.363504e-09 
#>                             mesoderm development 
#>                                     2.305281e-07 
#>                          cell fate specification 
#>                                     2.734734e-07 
#>                 embryonic heart tube development 
#>                                     9.498289e-07 
#>                       cardiocyte differentiation 
#>                                     4.236265e-06 
#>                         diencephalon development 
#>                                     4.506317e-06 
#>                              heart morphogenesis 
#>                                     9.151835e-06 
#>                             endoderm development 
#>                                     3.238483e-05 
#>                                heart development 
#>                                     3.238483e-05 
#> regulation of endodermal cell fate specification 
#>                                     3.238483e-05
par(mar=c(10,5,4,2))
barplot(-log10(head(sort(significant.GO.Terms),10)), ylab = "- log 10 p-value", las=2)
```

![](README-unnamed-chunk-14-1.png)

Zebrafish cardiac development terms are significantly enriched in mouse cardiac development genes, and vice-versa.

Please now try your own analysis!


Trouble Shooting
------------

When Ensembl move servers (rarely) or a particular server is down (less rare) you can use the following command to change servers once the xgsa package has been loaded:

assignInNamespace(x = "hostID", ns = "xgsa", value = "asia.ensembl.org")

Or for a stable archive version:

assignInNamespace(x = "hostID", ns = "xgsa", value = "aug2020.archive.ensembl.org")
