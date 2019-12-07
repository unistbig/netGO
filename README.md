# netGO
<img src = 'https://user-images.githubusercontent.com/6457691/70369320-d30da980-18fa-11ea-9b9e-d3eaf7400a4e.png'></img>


### Introduction
netGO is an R/Shiny package for network-integrated pathway enrichment analysis.<br>
It also provides the conventional Fisher’s exact test.<br> 
Specifically, it provides user-interactive visualization of enrichment analysis results and related networks.<br>
Currently, netGO provides network and annotation gene-set data for 4 species<br> 
Including Human, Mouse, Yeast, and Arabidopsis.<br> 
These data are available from this repository [netGO-Data](https://github.com/unistbig/netGO-Data/)<br>

### Install and Example codes

#### Package Dependency : Please install following R packages before use netGO
- devtools, doParallel, doSNOW, DT, foreach, googleVis, htmlwidgets, parallel, shiny, shinyCyJS, shinyJS, V8

```r
install.packages(
  c('devtools','doParallel','doSNOW','DT',
  'foreach','googleVis','htmlwidgets',
  'parallel','shiny','shinyjs','V8')
)
library(devtools)
install_github('unistbig/shinyCyJS') 
```

#### Example Run

```r
library(devtools) # load devtools to use ‘install_github’ function
install_github('unistbig/netGO') # install netGO
library(netGO) # load netGO
DownloadExampleData() # Download and load example breast tumor datasets

obj = netGO(genes = brca[1:20], genesets, network, genesetV) # run netGO
# or you can load pre calculated result using this command 
# load("brcaresult.RData") 
netGOVis(obj, genes = brca[1:20], genesets, network) # visualize netGO's result
```
This example run takes around 10 min on Desktop, and 15 - 25 min on Laptop ( not including downloading data )<br>

### Main functions
<hr>

#### 1. netGO
This function returns a data frame of gene-set Q-values derived from netGO+, netGO and Fisher’s exact test.<br>

#### Input arguments

-	genes: A character vector of input genes (e.g., DE genes).<br>

-	genesets: A list of gene-sets consisting of groups of genes.<br>

-	network: A numeric matrix of network data. The network score range is [0,1].<br>

-	genesetV: A numeric matrix of pre-calculated interaction data between gene and gene-sets.<br>
The matrix dimension must be [ {# of genes} X {# of gene-sets}]. <br>
It can be built using BuildGenesetV function with network and gene-set objects as input arguments.

```r
genesetV = BuildGenesetV(network, genesets)
```

<img src = 'https://user-images.githubusercontent.com/6457691/70369353-43b4c600-18fb-11ea-8318-f37aaafb9b17.png'></img>

-	alpha (optional): A numeric parameter weights how much network score will be effected. <br>
The value is positive numeric value with > 1 and the default is 20.<br>

- beta (optional): A numeric parameter balancing the weights between the relative and absolute network scores.<br> 
The value is between 0 and 1 and the default is 0.5.<br>

-	nperm: The number of permutations.<br>

<img src = 'https://user-images.githubusercontent.com/6457691/70369534-7e6c2d80-18fe-11ea-877e-3aa6c79cd4f1.png'></img>

**Notice** member of genes should be gene symbols when using the default STRING and mSigDB data. <br>
Other types of gene names are also available if the corresponding customized data (network and gene-set data) are used.
<hr>

#### 2. netGOVis (for visualization)
This function visualizes the result on the web browser (google chrome is recommended).<br> 
The result graphs and table are downloadable from the web browser.<br>

#### Input arguments
-	obj: A result data frame derived from ‘netGO’ function.<br>
It consists of 6 columns including 1) gene-set name and Q-values evaluated from 2) netGO , 3) netGO+, and 4) Fisher’s exact test.<br>
and additional 5) Overlap and 6) Network scores<br>

-	R (optional): Gene-set rank threshold, default is NULL (if 50, Top 50 gene-sets in either method will be shown).<br>

-	Q (optional): Gene-set Q-value threshold, default is 0.25. (Gene-sets with Q-value <= 0.25 will be used)<br>

-	genes, genesets, network: same as in the ‘netGO’ function.<br>
<hr>

<img src = 'https://user-images.githubusercontent.com/6457691/70369561-ee7ab380-18fe-11ea-9dcc-fe03d0ea37f0.png'></img>

*number after listening on may different*

<img src = 'https://user-images.githubusercontent.com/6457691/70369640-09015c80-1900-11ea-9eb3-f825e2cbf511.png'></img>


Question / Comment / Suggest : kjh0530@unist.ac.kr <br>
Thanks
