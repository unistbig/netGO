# :yellow_heart: netGO <img src ='https://user-images.githubusercontent.com/6457691/71441958-1a8f9480-2747-11ea-8250-95f481180073.png' align = 'right' width = 150></img>


netGO is an R/Shiny package for network-integrated pathway enrichment analysis.<br>
netGO provides user-interactive visualization of enrichment analysis results and related networks.<br>

Currently, netGO provides network and annotation gene-set data for 4 species. (*[Arabidopsis](https://github.com/unistbig/netGO-Data/tree/master/Arabidopsis), [Human](https://github.com/unistbig/netGO-Data/tree/master/Human), [Mouse](https://github.com/unistbig/netGO-Data/tree/master/Mouse), and [Yeast](https://github.com/unistbig/netGO-Data/tree/master/Yeast)*)<br>
These data are available from [netGO-Data](https://github.com/unistbig/netGO-Data) repository.<br>

## :clipboard: Prerequisites

these R packages should be installed before run netGO. (Alphabetical ordered)<br>

*devtools, doParallel, doSNOW, DT, foreach, googleVis, htmlwidgets, shiny, shinyCyJS, shinyjs, V8*

most of packages are avaiable with [CRAN](https://cran.r-project.org/).<br>
* but [shinyCyJS](https://github.com/unistbig/shinyCyJS) needs to be install with github.<br>

```R
library(devtools)
install_github('unistbig/shinyCyJS')
library(shinyCyJS)
```

* linux user should install V8 after install additional packages.<br>

On Debian / Ubuntu : libv8-dev or libnode-dev. <br>
On Fedora : v8-devel <br>
[more information](https://cran.r-project.org/web/packages/V8/index.html)

if you need, use this code.

``` R
install.package('devtools') # 2.2.1
install.package('doParallel') # 1.0.15
install.package('doSNOW') # 1.0.18
install.package('DT') # 0.11
install.package('foreach') # 1.4.7
install.package('googleVis') # 0.6.4
install.package('htmlwidgets') # 1.5.1
install.package('shiny') # 1.4.0
install.package('shinyjs') # 1.0
install.package('V8') # 2.3
```

## :wrench: Example Run

We prepared example run with breast tumor dataset from *NCBI GEO [GSE3744](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE3744)<br>*

```r  
library(netGO) # load netGO library
DownloadExampleData() # Download and load example breast tumor datasets

# run netGO, this may takes a time.
# or user can load pre-calculated result using this command 
# load("brcaresult.RData")   

obj = netGO(genes = brca[1:30], genesets, network, genesetV) 

```
  
This example run takes around ~ 10 mins on Desktop,<br>
and 15 - 25 min on Laptop ( not included time for download data )<br>  

user can briefly see the netGO's analysis result.
  
```r
head(obj)
```

<img src ='https://user-images.githubusercontent.com/6457691/70370939-f5f68880-190f-11ea-9615-b11fb789fd0e.png'></img>

also user can visualize result with netGOVis function.
  
```r  
netGOVis(obj, genes = brca[1:30], genesets, network, R = 50, Q = 0.25 ) # visualize netGO's result
```

<br>  
<img src = 'https://user-images.githubusercontent.com/6457691/70369640-09015c80-1900-11ea-9eb3-f825e2cbf511.png'></img>

if user want to access result without shinyweb-application, <br>
user can export this result as text format using *exportGraphTxt, exportGraph, exportTable* functions.<br>
more description explained in API term.<br>

here are example codes.<br>

```r 
# exportGraphTxt
table = exportGraphTxt(gene = brca[1:30], geneset = genesets[['SMID_BREAST_CANCER_NORMAL_LIKE_UP']], network) # table
head(table)

# exportGraph
graph = exportGraph(brca[1:30], geneset = genesets[['SMID_BREAST_CANCER_NORMAL_LIKE_UP']], network) # shinyCyJS graph object
shinyCyJS(graph)

# exportTable
table = exportTable(obj, R = 50, Q = 0.25) # table
head(table)

dtable = exportTable(obj, type='D', R = 50, Q = 0.25) # data.table
dtable

``` 

## :memo: Datas

### <b>Example Datasets</b> [(netGO-Data repository)](https://github.com/unistbig/netGO-Data) <br>  
 
#### Human (in human directory)

|Data|Genes|Genesets|Network|GenesetV|
|:---:|:---:|:---:|:---:|:---:|
|Breast Tumor|brca.RData|c2gs.RData|networkString.RData networkHumannet.RData|genesetVString1,2.RData genesetVHumannet1,2.RData|
|P53|p53.RData|c2gs.RData|networkString.RData networkHumannet.RData|genesetVString1,2.RData genesetVHumannet1,2.RData|
|Diabete|dg.RData|cpGenesets.RData|networkString.RData networkHumannet.RData|cpgenesetV1,2.RData|

user can use Breast tumor data with *DownloadExampleData* function in netGO (Recommended)

#### Arabidopsis

|Data|Genes|Genesets|Network|GenesetV|
|:---:|:---:|:---:|:---:|:---:|
|ShadowResponse|Aragenes.RData|KEGGara.RData|networkAranet.RData|AragenesetV.RData|

#### Mouse & Yeast ( geneset, network available )

|Species|Genesets|Network|
|:----:|:----:|:----:|
|Mouse|KEGGmouse.Rdata|networkMousenet.Rdata|
|Yeast|KEGGyeast.Rdata|networkYeastnet.Rdata|

### <b>Data Formats</b>

netGO needs 4 data.<br>
- *genes* : A character vector of input genes (e.g., DE genes).<br>
- *genesets* : A list of gene-sets consisting of groups of genes.<br>
- *network* : A numeric matrix of network data. The network score range is normalized in scale [0,1].<br>
- *genesetV* : A numeric matrix of pre-calculated interaction data between gene and gene-sets.<br>
  The dimension of matrix must be [ {number of genes} , {number of gene-sets}]. <br>
  It can be built by using *BuildGenesetV* function with network and gene-set objects as input arguments.

  ```r
  genesetV = BuildGenesetV(network, genesets)
  ```  

## :white_circle: Functions
<hr>

### 1. netGO
netGO function calculates significant gene and genesets<br>
and returns a data frame of Gene-sets, and their P-value, Q-values derived from netGO+, Fisher’s exact test and *netGO (optional)*.<br>
and score for network interaction and overlap <br>

<b>Input arguments</b>
 
* genes: A character vector of input genes (e.g., DE genes).<br>
* genesets: A list of gene-sets consisting of groups of genes.<br>
* network: A numeric matrix of network data. The network score range is [0,1]. 1 for strong interaction and 0 for no interaction<br>
  
  | |A|B|C|
  |:--:|:--:|:--:|:--:|
  |A|0|0.1|0.76|
  |B|0.1|0|0.324|
  |C|0.76|0.324|0|
  
* genesetV: A numeric matrix of pre-calculated interaction data between gene and gene-sets.<br>
  this can be built with *BuildGenesetV* function.
  
  | |Gene-set1|Gene-set2|Gene-set3|
  |:--:|:--:|:--:|:--:|
  |A|0.837|1.647|0.074|
  |B|0|1.75|0.113|
  |C|0.464|0.486|2.442|
    
* alpha (optional): A numeric parameter weights how much network score will be effected. <br>
  The value is positive numeric value with > 1 and the default is 20.<br>
* beta (optional): A numeric parameter balancing the weights between the relative and absolute network scores.<br> 
  The value is between 0 and 1 and the default is 0.5.<br>

<img src = 'https://user-images.githubusercontent.com/6457691/71439289-61c45800-273c-11ea-86db-a06bec993486.png' width = 500></img>

* nperm (optional): A numeric parameter to tell netGO how many times permutation will be repeated, default is 10000<br>
* category (optional): A numeric parameter to decide number of genes will be categorized while resampling. <br>
  default is NULL. and it will divide ~ 2000 genes to each category.<br>
* pvalue (optional): A boolean parameter to tell whether return Q-values only ( FALSE ) or both P-values and Q-values (TRUE)
* plus (optinoal): A boolean parameter to tell whether netGO will calculate both netGO, netGO+ (plus = FALSE) or netGO+ only ( plus = TRUE, default )<br>
* verbose (optional) : A boolean parameter to show more process of netGO.<br>

for more information about alpha, beta and category, please refer our manuscript. 

After Example run, user may see these logs in R console.

<img src='https://user-images.githubusercontent.com/6457691/71439716-0e530980-273e-11ea-8b27-3621c90416cc.png'></img>

**Notice** member of genes should be gene symbols when using the default STRING and mSigDB data. <br>
Other types of gene names are also available if the corresponding customized data (network and gene-set data) are used.
<hr>

### 2. netGOVis 

netGOVis function visualizes netGo's result on the web browser (google chrome is recommended).<br> 
The result graphs *(svg format)* and table are downloadable from the web browser.<br>

<b>Input arguments</b>

* obj: A result data frame derived from <b>netGO</b> function.<br>
It consists of multiple columns including <br>
1) gene-set name and <br>
p, q-values evaluated from <br>
2) *netGO (optional)*, <br>
3) netGO+, and <br>
4) Fisher’s exact test.<br>
and additional<br>
5) Overlap and <br>
6) Network scores<br>

* genes, genesets, network: same as in the *netGO* function.<br> 
* R (optional): Gene-set rank threshold, default is 50 (Top 50 gene-sets in either method will be shown).<br>
* Q (optional): Gene-set Q-value threshold, default is 0.25. (Gene-sets with Q-value <= 0.25 will be used)<br>

After Example run, user may see these logs in R console.

<img src = 'https://user-images.githubusercontent.com/6457691/71439835-6a1d9280-273e-11ea-922a-06bf35c45658.png'></img>
            
and user's default web browser (<b>we built netGO based on chrome environment</b>) will return this interactive results.

<img src = 'https://user-images.githubusercontent.com/6457691/71439882-946f5000-273e-11ea-952f-6eda92b3f090.png' width = 700></img>
                        
<hr>

### 3. BuildGenesetV

BuildGenesetV function will build genesetV object with given network and geneset. <br>
genesetV is pre-calculated interaction files to reduce running time in netGO.

<b>Input arguments</b>

* genesets, network: same as in the *netGO* function.<br> 
<hr>

### 4. DownloadExampleData()

DownloadExampleData function will <br> 
1) Download example dataset object in user's working directory and <br>
2) load Example dataset ( breast tumor, [GSE3744](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE3744) ) in user's R environment.<br>
note that, if objects exist in working directory, it'll not download again, <br>
but we recommend that, remove and re-download them when netGO package updated.

<b>Input arguments</b>
* none
* R object named *brca, genesets, genesetV, network, obj* will be loaded. 
<hr>

### 5. exportGraph()

exportGraph function will export interaction in netGO's result object as graph object that can be accessed with *[shinyCyJS](https://github.com/unistbig/shinyCyJS)* function.<br>

<b>Input arguments</b>
* genes, network : same as in the *netGO* function.<br> 

* geneset : A character vector of gene symbols (e.g., one member of genesets object in *netGO*).<br>

for example, 

``` R
geneset = genesets[['SMID_BREAST_CANCER_NORMAL_LIKE_UP']]
graph = exportGraph(brca[1:30], geneset = genesets[['SMID_BREAST_CANCER_NORMAL_LIKE_UP']], network) # shinyCyJS graph object

shinyCyJS(graph)

```
note that, see graph with R's default viewer (not web browser), will not use layout functions.<br>
so result will have different shape. <br>

<img src='https://user-images.githubusercontent.com/6457691/71440782-eb2a5900-2741-11ea-99ad-1e147e5e142b.gif' width = 500></img>
<hr>

### 6. exportGraphTxt()

exportGraphTxt function will export interaction in netGO's result object as table format.<br>

<b>Input arguments</b>
* genes, network, geneset : same as in the *exportGraph* function.<br> 

for example, 

``` R
table = exportGraphTxt(brca[1:30], geneset, network)
head(table)
```

it will have format of 

|geneA|geneB|strength|type|
|:---:|:---:|:---:|:---:|
|A|B|0.1|Inter|
|C|D|0.82|Inner|

which means, there is Interaction in A & B gene with 0.1 strength ( 0 is weaker, 1 is stronger )<br>
and Inter means A & B gene are not overlapped between *gene* and *geneset*. <br>

also there is Interaction in C & D gene with 0.82 strength<br>
while Inner means C & D gene are overlapped between *gene* and *geneset*. <br>
<hr>

### 7. exportTable()

exportTable will export netGO's result object as table or data.table.

<b>Input arguments</b>
* obj, R, Q : same as in the *netGOVis* function.<br> 

for example, 

```R
table = exportTable(obj, R = 50, Q = 0.25) # table
head(table)

dtable = exportTable(obj, type='D', R = 50, Q = 0.25) # data.table
dtable
```

result will have format of 

|geneset name|netGO+ q-value| Fisher q-value |
|:---:|:---:|:---:|
|genesetA|0.11|0.2|

<hr>

## :blue_book: Explore with netGO

visualization page of netGO consists of 3 parts : Network, Table, Bubble.

### <b>Network</b>

* Network panel displays the input genes, selected gene-set, and the network connections between the two.  
* ![#48dbfb](https://placehold.it/15/48dbfb/000000?text=+) Sky blue nodes represent input genes (e.g., DE genes) 
* ![#feca57](https://placehold.it/15/feca57/000000?text=+) Yellow nodes represent genes in the selected gene-set 
* ![#1dd1a1](https://placehold.it/15/1dd1a1/000000?text=+) Green nodes represent the intersection of input genes and the gene-set. 
* The edge width represents the strength of interaction between two nodes. 
* ![#feca57](https://placehold.it/15/feca57/000000?text=+) Gene-set Nodes without edges will be discarded.
* The gene-set can be selected by clicking on the gene-set name in the Table on the right side. 
* The users can download the graph image as SVG format.

 <img src = 'https://user-images.githubusercontent.com/6457691/70425850-4e1abf80-1ab5-11ea-96a3-84ac3e82d9f9.png' width = 900></img>

### <b> Table </b>

* Table contains the names of gene-sets and their Q-values ( or P-values ) evaluated from netGO, netGO+ and Fisher’s exact test, respectively. It is downloadable by clicking the ‘Download Table’ button in the upper right of the table 
 <img src ='https://user-images.githubusercontent.com/6457691/70425051-c84a4480-1ab3-11ea-8eb4-4b45385943fe.png'></img>
 
### <b> Bubble </b>
 
* Bubble module plots the bubble chart of significant gene-sets. <br>
* The overlap (x-axis) and network scores (y-axis) of each significant gene-sets are represented. <br> 
* The size of bubbles represents the significance level of each gene-set in -log10 scale (Qvalue).
* Hovering/Click on each bubble will show statistcial values.
 
<img src ='https://user-images.githubusercontent.com/6457691/70425757-1c095d80-1ab5-11ea-99f6-4198fa48b384.png'></img>
 
## :blush: Authors
* :octocat: Jinhwan Kim [@jhk0530](http://github.com/jhk0530)

* comment / suggest / question will be really appreciated, *kjh0530@unist.ac.kr*

## :memo: License
Copyright :copyright: 2019 Jinhwan Kim
This project is [MIT](https://opensource.org/licenses/MIT) licensed

*This README was generated with :two_hearts: by [shinyReadme](http://github.com/jhk0530/shinyReadme)*
