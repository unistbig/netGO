# netGO
R/Shiny package for network-integrated pathway enrichment analysis

### Install and Launch
#### * Prerequisites : Please install following R packages
- devtool, Rcpp, shiny, shinyjs, DT, doParallel, foreach, parallel, htmlwidgets, googleVis, V8, shinyCyJS
- Tip: <b>shinyCyJS</b> can be installed by  

```r
install_github(‘unistbig/shinyCyJS’)
```

#### R codes to run netGO

```r
library(devtools) # load devtools to use 'install_github' function
install_github('unistbig/shinyCyJS') # install shinyCyJS
install_github('unistbig/netGO') # install netGO
library(netGO) # load netGO
obj = netGO(genes, genesets, PPI, genesetV) # run netGO and save the result in 'obj' object
netGOVis(obj, genes, genesets, R, Q, PPI) # visualize the result
```

### Example (it may takes ~10 minutes to run)
#### To run toy example, please download following data from https://github.com/unistbig/netGO-Data/tree/master/Human .

* PPIString.RData (PPI data)

* brca.RData (DE gene list)

* c2gs.RData (gene-set list)

* genesetVString.z01~.zip (pre-calculated interaction data)

Tip: <b>genesetVString.z01~.zip</b> are split zipped files. Uncompressing <b> ‘genesetVString.zip’</b> will create ‘genesetVString.RData’. 

```r
> library(devtools) 
> install_github('unistbig/shinyCyJS') 
> install_github('unistbig/netGO')
> library(netGO) 
> load('c2gs.RData') # contains genesets
> load('brca.RData') # contains brca
> load('genesetVString.RData') # contains genesetV
> load('PPIString.RData') # contains PPI
> obj = netGO(genes = brca[1:10], genesets = genesets, PPI = PPI, genesetV = genesetV)
```
