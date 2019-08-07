# netGO
R/Shiny package for network-integrated pathway enrichment analysis

### Introduction
netGO is an R/Shiny package for network-integrated pathway enrichment analysis. It provides gene-set enrichment analysis results evaluated from netGO and simple Fisher’s exact test. It also provides visualization to let the users easily compare the results derived from two methods. The R/Shiny package is available at Github ( https://github.com/unistbig/netGO ). Currently, netGO provides network and gene-set data of four species including human, mouse, yeast and arabidopsis. Such accessory data are available at another repository (https://github.com/unistbig/netGO-Data/)

### Install and Launch
#### * Prerequisites : Please install following R packages
- devtool, Rcpp, shiny, shinyjs, DT, doParallel, foreach, parallel, htmlwidgets, googleVis, V8, shinyCyJS

```r
install.packages(c('devtools', 'Rcpp', 'shinyjs', 'DT','doParallel', 'foreach', 'parallel', 'htmlwidgets', 'googleVis', 'V8'))
library(devtools)
install_github('unistbig/shinyCyJS') # R >=3.5.2
```

#### Example Run

```r
library(devtools) # load devtools to use ‘install_github’ function
install_github('unistbig/netGO') # install netGO
library(netGO) # load netGO
DownloadExampleData() # Download and load example datasets contatining brca, genesets, PPI and genesetV
obj = netGO(genes = brca[1:20], genesets = genesets, network = network, genesetV = genesetV) # run netGO
netGOVis(obj = obj, genes = brca[1:20], genesets = genesets, R = 50, network = network) # Visualize the result

```
### Main functions

<hr>

#### 1. netGO
This function returns a data frame of gene-set p-values derived from netGO and Fisher’s exact test.
#### Input arguments

-	genes: A character vector of input genes (e.g., DE genes). <br>
-	genesets: A list of gene-sets consisting of groups of genes.<br>
-	network: A numeric matrix of network data. The network score range is [0,1].<br>
-	genesetV: A numeric matrix of pre-calculated interaction data between gene and gene-sets. The matrix dimension must be [ {# of genes} X {# of gene-sets}]. It can be built using BuildGenesetV function with network and gene-set objects as input arguments.
<br>

```r
genesetV = BuildGenesetV(network, genesets)
```
-	alpha (optional): A numeric parameter reflecting the influence of the network data (e.g., PPI network). The value is between 0~1 and default is 0.5.<br>
-	nperm: The number of permutations.<br>

**Notice** that, member of genes should be gene symbols when using the default STRING and mSigDB data. Other types of gene names are also available if the corresponding customized data (network and gene-set data) are used.
<hr>

#### 2. netGOVis (for visualization)
This function visualizes the result on the web browser (google chrome is recommended). The result graphs and table are downloadable from the web browser.<br>

#### Input arguments

-	obj: A result data frame derived from ‘netGO’ function. It consists of three columns including 1) gene-set name and p-values evaluated from 2) netGO (netGOP) and 3) Fisher’s exact test (FisherP).<br>
-	R (optional): Gene-set rank threshold, default is 50 (Top 50 gene-sets in either method will be shown).<br>
-	Q (optional): Gene-set Q-value threshold, default is 1. (Gene-sets with Q-value <= 1 will be shown)<br>
-	genes, genesets, network: same as in the ‘netGO’ function.<br>
<hr>
