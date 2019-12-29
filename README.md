netGO is an R/Shiny package for network-integrated pathway enrichment analysis.<br>
netGO provides user-interactive visualization of enrichment analysis results and related networks.<br>

Currently, netGO supports analysis for four species (*[Human](https://github.com/unistbig/netGO-Data/tree/master/Human), [Mouse](https://github.com/unistbig/netGO-Data/tree/master/Mouse), [Arabidopsis thaliana](https://github.com/unistbig/netGO-Data/tree/master/Arabidopsis),and [Yeast](https://github.com/unistbig/netGO-Data/tree/master/Yeast)*)<br>
These data are available from [netGO-Data](https://github.com/unistbig/netGO-Data) repository.<br>

## :clipboard: Prerequisites
The R packages listed below are required to be installed before running netGO.(Alphabetical order)<br>

*devtools, doParallel, doSNOW, DT, foreach, googleVis, htmlwidgets, shiny, shinyCyJS, shinyjs, V8*

* Most of the packages are avaiable from[CRAN](https://cran.r-project.org/), but [shinyCyJS](https://github.com/unistbig/shinyCyJS) needs to be installed from github.<br>

* Linux user has to install V8 after installing the other packages.<br>

* Note that netGO is not supported for centOS8, because V8 is not available in centoOS8.<br>

On Debian / Ubuntu : libv8-dev or libnode-dev. <br>
On Fedora : v8-devel <br>
[more information](https://cran.r-project.org/web/packages/V8/index.html)


The user may want to use the following codes to install the required packages.

``` R
install.packages('devtools') # 2.2.1
library(devtools)
install_github('unistbig/shinyCyJS')
install.packages('doParallel') # 1.0.15
install.packages('doSNOW') # 1.0.18
install.packages('DT') # 0.11
install.packages('foreach') # 1.4.7
install.packages('googleVis') # 0.6.4
install.packages('htmlwidgets') # 1.5.1
install.packages('shiny') # 1.4.0
install.packages('shinyjs') # 1.0
install.packages('V8') # 2.3
```

## :wrench: Running with an example data

Here are codes to run netGO for the breast tumor dataset (*GEO [GSE3744](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE3744).*)<br>

```r  
library(netGO) # load netGO library
DownloadExampleData() # Download and load the breast tumor data
obj = netGO(genes = brca[1:30], genesets, network, genesetV) 

# The user may also load the pre-calculated result using the following command
# load("brcaresult.RData")   
```
  
Running this exmaple takes 5 to 25 minutes depending on the system used. The analysis retuls of netGO is shown below.<br>
 
<img src ='https://user-images.githubusercontent.com/6457691/70370939-f5f68880-190f-11ea-9615-b11fb789fd0e.png'></img>

The analysis result can be visualized using the following codes:<br>  
```r  
netGOVis(obj, genes = brca[1:30], genesets, network, R = 50, Q = 0.25 ) # visualize netGO's result
```

<br>  
<img src = 'https://user-images.githubusercontent.com/6457691/70369640-09015c80-1900-11ea-9eb3-f825e2cbf511.png'></img>

If user wants to access result without shinyweb-application, the following functions can be used to export the result as text files<br>

```r 
# exportGraphTxt
table = exportGraphTxt(gene = brca[1:30], geneset = 
genesets[['SMID_BREAST_CANCER_NORMAL_LIKE_UP']], network) # table
head(table)

# exportGraph
graph = exportGraph(brca[1:30], geneset = 
genesets[['SMID_BREAST_CANCER_NORMAL_LIKE_UP']], network) # shinyCyJS graph object
shinyCyJS(graph)

# exportTable
table = exportTable(obj, R = 50, Q = 0.25) # table
head(table)

dtable = exportTable(obj, type='D', R = 50, Q = 0.25) # data.table
dtable
``` 

## :memo: Data

### <b>Example Datasets</b> [(netGO-Data repository)](https://github.com/unistbig/netGO-Data) <br>  
 
#### Human 

|Data|Genes|Genesets|Network|GenesetV|
|:---:|:---:|:---:|:---:|:---:|
|Breast Tumor|brca.RData|c2gs.RData|networkString.RData networkHumannet.RData|genesetVString1,2.RData genesetVHumannet1,2.RData|
|P53|p53.RData|c2gs.RData|networkString.RData networkHumannet.RData|genesetVString1,2.RData genesetVHumannet1,2.RData|
|Diabetes|dg.RData|cpGenesets.RData|networkString.RData networkHumannet.RData|cpgenesetV1,2.RData|

The user can download the  breast tumor data using *DownloadExampleData* function(Recommended)

#### Arabidopsis thaliana

|Data|Genes|Genesets|Network|GenesetV|
|:---:|:---:|:---:|:---:|:---:|
|ShadowResponse|Aragenes.RData|KEGGara.RData|networkAranet.RData|AragenesetV.RData|

#### Mouse & Yeast ( gene-set and networks available )

|Species|Genesets|Network|
|:----:|:----:|:----:|
|Mouse|KEGGmouse.Rdata|networkMousenet.Rdata|
|Yeast|KEGGyeast.Rdata|networkYeastnet.Rdata|

### <b>Data Formats</b>

netGO requires the follwoing four data types.<br>
- *genes* : a character vector of input genes (e.g., differentially expressed genes).<br>
- *genesets* : a list of gene-sets cto be tested.<br>
- *network* : a numeric matrix of network data. The network scores are normalized to the unit interval [0,1] by dividing each score by the maximum score<br>
- *genesetV* : A numeric matrix of pre-calculated interaction data between gene and gene-sets.<br>
  The dimension of matrix must be [{number of genes} , {number of gene-sets}]. <br>
  It can be built by using *BuildGenesetV* function with network and genesets objects as the input arguments.

  ```r
  genesetV = BuildGenesetV(network, genesets)
  ```  

## :white_circle: Functions
<hr>

### 1. netGO
netGO function tests the significance of the gene-sets for the input gene list<br>
and returns a data frame of gene-sets, their *p*-values, *q*-values derived from netGO+, Fisher’s exact test and netGO (optional) as well as the scores for the network interaction and overlap.<br>

<b>Input arguments</b>
 
* genes: a character vector of input genes (e.g., differentially expressed genes).<br>
* genesets: a list of gene-sets consisting of groups of genes.<br>
* network: A numeric matrix of network data. The network scores are normalized to the unit interval [0,1]. 1 represents strong interaction and 0 for no interaction<br>
  
  | |A|B|C|
  |:--:|:--:|:--:|:--:|
  |A|0|0.1|0.76|
  |B|0.1|0|0.324|
  |C|0.76|0.324|0|
  
* genesetV: a numeric matrix of pre-calculated interaction data between genes and gene-sets.<br>
  This object can be built with *BuildGenesetV* function.
  
  | |Gene-set1|Gene-set2|Gene-set3|
  |:--:|:--:|:--:|:--:|
  |A|0.837|1.647|0.074|
  |B|0|1.75|0.113|
  |C|0.464|0.486|2.442|
    
* alpha (optional): a numeric parameter ( ≥ 1; the default is 20) that weights the contribution of network connections in enrichment analysis.<br>
  
* beta (optional): a numeric parameter (∈[0,1]; the default is 0.5) that balances the weights between the relative and absolute network scores.<br>

<img src = 'https://user-images.githubusercontent.com/6457691/71439289-61c45800-273c-11ea-86db-a06bec993486.png' width = 500></img>

* nperm (optional): a numeric parameter to determine the bin size (number of genes) to be used during resampling. The default is NULL which assigns approximately 2000 genes to each bin<br>
  default is NULL. and it will divide ~ 2000 genes to each category.<br>
* pvalue (optional): a boolean parameter to determine whether to return Q-values only ( FALSE ) or both P-values and Q-values (TRUE)<br>
* plus (optinoal): a boolean parameter to determine whether to run both netGO and netGO+ (plus = FALSE) or netGO+ only ( plus = TRUE, default )<br>
* verbose (optional) : a boolean parameter<br>

After running the example, the user may see the following logs in R console.<br>

<img src='https://user-images.githubusercontent.com/6457691/71439716-0e530980-273e-11ea-8b27-3621c90416cc.png'></img>

**Notice** the input genes should be represented in **gene symbols** when using the default networks and gene-sets (STRING and MSigDB). <br>
Other types of gene names are also allowed if the corresponding customized data (networks and gene-set data) are used. <br>

<hr>

### 2. netGOVis 

netGOVis function visualizes the analysis results on the web browser (google chrome is recommended).<br> 
The resulting graphs (svg format) and table are downloadable from the web browser.<br>

<b>Input arguments</b>

* obj: the data frame of analysis results obtained by running netGO function.<br>
It consists of multiple columns including <br>

1.	gene-set name and p, q-values evaluated using netGO (optional), netGO+, and Fisher’s exact test as well as the scores for the overlap and networks.<Br>

* genes, genesets, network: the same as those in the *netGO* function.<br> 
* R (optional): gene-set rank threshold, The default is 50 (Top 50 gene-sets in either method will be shown).<br>
* Q (optional): Gene-set Q-value threshold, The default is 0.25. (gene-sets with Q-value ≤ 0.25 will be used)<br>

Afte running the netGO function, the user may see the following logs in the R console.

<img src = 'https://user-images.githubusercontent.com/6457691/71439835-6a1d9280-273e-11ea-922a-06bf35c45658.png'></img>
            
and user's default web browser (<b>netGO was built based on chrome environment</b>) will return the following interactive visualization:

<img src = 'https://user-images.githubusercontent.com/6457691/71439882-946f5000-273e-11ea-952f-6eda92b3f090.png' width = 700></img>
                        
<hr>

### 3. BuildGenesetV

BuildGenesetV function will build genesetV object using the given *network* and *genesets*. <br>
genesetV is pre-calculated interaction files used to reduce the running time of netGO.

<b>Input arguments</b>

* genesets, network: the same as those in the *netGO* function.<br> 
<hr>

### 4. DownloadExampleData()

This function will download example data in the user's working directory and load the data ( breast tumor, [GSE3744](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE3744) ) in user's R environment.<br>
Note that, if objects exist in the working directory, this function will not download the data again, so we recommand removing and downloading them again if netGO package is updated. <br>

<b>Input arguments</b>
* none
* R object named *brca, genesets, genesetV, network, obj* will be loaded. 
<hr>

### 5. exportGraph()

exportGraph function will export network data from the netGO analsysis result as graph object<br>

<b>Input arguments</b>
* genes, network : the same as those in the *netGO* function.<br> 

* geneset : a character vector of gene symbols (e.g., member of genesets object in *netGO*).<br>

for example, 

``` R
geneset = genesets[['SMID_BREAST_CANCER_NORMAL_LIKE_UP']]
graph = exportGraph(brca[1:30], geneset = 
genesets[['SMID_BREAST_CANCER_NORMAL_LIKE_UP']], network) # shinyCyJS graph object

shinyCyJS(graph)
```

However, the default viewer of R (not web browser) will not use the layout functions as shown below.<br>

<img src='https://user-images.githubusercontent.com/6457691/71440782-eb2a5900-2741-11ea-99ad-1e147e5e142b.gif' width = 500></img>
<hr>

### 6. exportGraphTxt()

exportGraphTxt function will export network data from the netGO analysis result as table format.<br>

<b>Input arguments</b>
* genes, network, geneset : the same as those in the *exportGraph* function.<br> 

For example, 

``` R
table = exportGraphTxt(brca[1:30], geneset, network)
head(table)
```

the exported data are shown as

|geneA|geneB|strength|type|
|:---:|:---:|:---:|:---:|
|A|B|0.1|Inter|
|C|D|0.82|Inner|

which means, there is interaction between genes A and B with the strength of 0.1 ( 0 means no interaction, 1 means strong interaction )<br>
and 'Inter' means A and B are not overlapped between *gene* and *geneset*.<br>

there is also interaction between genes C and D with the strength of 0.82<br>
while 'Inner' means C and D are overlapped between *gene* and *geneset*.

<hr>

### 7. exportTable()

exportTable will export the result object of netGO as table or data.table.

<b>Input arguments</b>
* obj, R, Q : the same as those in the *netGOVis* function.<br> 

for example, 

```R
table = exportTable(obj, R = 50, Q = 0.25) # table
head(table)

dtable = exportTable(obj, type='D', R = 50, Q = 0.25) # data.table
dtable
```

The exported data have the format as follows:<br>

|geneset name|netGO+ q-value| Fisher q-value |
|:---:|:---:|:---:|
|genesetA|0.11|0.2|

<hr>

## :blue_book: Visualization and exploration of netGO analysis results

The netGO analysis results are visualized through three panels: interaction networks, list of significant gene-sets, and the bubble chart.<br>

### <b>Interaction Network</b>

* The network panel displays the input genes, selected gene-set, and the network connections between the two.
* ![#48dbfb](https://placehold.it/15/48dbfb/000000?text=+) Sky blue nodes represent input genes (e.g., differentially expressed genes) 
* ![#feca57](https://placehold.it/15/feca57/000000?text=+) Yellow nodes represent genes in the selected gene-set 
* ![#1dd1a1](https://placehold.it/15/1dd1a1/000000?text=+) Green nodes represent the intersection of input genes and the gene-set. 
* The edge width represents the strength of interaction between two nodes. 
* Genes without edges will be not be displayed.
*	The gene-set can be selected by clicking on the gene-set name on the upper-right panel.
*	The user can download the graph image as SVG format.

 <img src = 'https://user-images.githubusercontent.com/6457691/70425850-4e1abf80-1ab5-11ea-96a3-84ac3e82d9f9.png' width = 900></img>

### <b>Significant gene-sets</b>

* This panel contains the list of significant gene-sets as well as their Q-values ( or P-values ) evaluated from netGO, netGO+ and Fisher’s exact test. It is downloadable by clicking the ‘Download Table’ button in the upper right corner of the table <br>

 <img src ='https://user-images.githubusercontent.com/6457691/70425051-c84a4480-1ab3-11ea-8eb4-4b45385943fe.png'></img>
 
### <b>Bubble chart</b>
 
* This module plots the bubble chart of significant gene-sets for the netGO+ results.
* The overlap (x-axis) and network (y-axis) scores of the significant gene-sets are represented.
*	The size of bubbles represents the significance level of each gene-set in -log10 scale (Qvalue).
* Hovering/Click on each bubble will show corresponding statistical values.
 
<img src ='https://user-images.githubusercontent.com/6457691/70425757-1c095d80-1ab5-11ea-99f6-4198fa48b384.png'></img>
 
## :blush: Contact

* Comments / suggestions and questions will be greatly appreciated,

* :octocat: Jinhwan Kim [@jhk0530](http://github.com/jhk0530) 

* prof. Dougu Nam *dougnam@unist.ac.kr*

## :memo: License

This project is [MIT](https://opensource.org/licenses/MIT) licensed

