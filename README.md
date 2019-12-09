# netGO User's Manual
<img width = 500 src = 'https://user-images.githubusercontent.com/6457691/70369320-d30da980-18fa-11ea-9b9e-d3eaf7400a4e.png'></img>
 @ 2019-12-09

### Introduction

netGO is an R/Shiny package for network-integrated pathway enrichment analysis.<br>
It also provides the conventional Fisher’s exact test.<br>
Specifically, it provides user-interactive visualization of enrichment analysis results and related networks.<br>
Currently, netGO provides network and annotation gene-set data for 4 species<br>
Including Human, Mouse, Yeast, and Arabidopsis.<br>
These data are available from this repository [netGO-Data](https://github.com/unistbig/netGO-Data) 
<br>


### Install and Example Run
<details>  
  <summary> 
    <b>Install Package</b>
  </summary>    
      
  - netGO uses other R packages, so please install following R packages before use netGO
  - devtools, doParallel, doSNOW, DT, foreach, googleVis, htmlwidgets, parallel, shiny, shinyCyJS, shinyJS, V8
  - User may want to use this codes.
  
  ```r
  install.packages(
    c('devtools','doParallel','doSNOW','DT','foreach','googleVis','htmlwidgets','parallel','shiny','shinyjs','V8')
  )
  library(devtools)
  install_github('unistbig/shinyCyJS') 
  ```  
</details>

<details>
  <summary>
    <b>Example Run</b>
  </summary>
  <br>

  We prepared example run with breast tumor dataset from NCBI GEO [GSE3744](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE3744)<br>

  ```r  
  library(netGO) # load netGO
  DownloadExampleData() # Download and load example breast tumor datasets
  obj = netGO(genes = brca[1:20], genesets, network, genesetV) # run netGO, this may takes a time.
  # or you can load pre calculated result using this command 
  # load("brcaresult.RData")   
  ```
  
  This example run takes around 10 min on Desktop, and 15 - 25 min on Laptop ( not including downloading data )
  <br>  
  user can briefly see the netGO's analysis result.
  
  ```r
  head(obj)
  ```

  <img src ='https://user-images.githubusercontent.com/6457691/70370939-f5f68880-190f-11ea-9615-b11fb789fd0e.png'></img>

  also user can visualize result with netGOVis function.
  
  ```r  
  netGOVis(obj, genes = brca[1:20], genesets, network) # visualize netGO's result
  ```
  
  <img src = 'https://user-images.githubusercontent.com/6457691/70369561-ee7ab380-18fe-11ea-9dcc-fe03d0ea37f0.png'></img>
  
  number after 'listening on' may different
  
  <br>  
  <img src = 'https://user-images.githubusercontent.com/6457691/70369640-09015c80-1900-11ea-9eb3-f825e2cbf511.png'></img>


</details>

### Datas
<details>
  <summary>     
    <b>Example Datasets</b> (https://github.com/unistbig/netGO-Data)    
  </summary>  
  <br>  
 
  |Species|Data Type|File name|Object Name|
  |:----:|:----:|:----:|:----:|
  |Human|Network|human/networkString.Rdata|network|
  |Human|Network|human/networkHumannet.Rdata|network|  
  |Human|Gene-set|human/c2gs.Rdata|genesets|
  |Human|Pre-calculated Interaction data|human/genesetVString1,2.Rdata|genesetV1,2|
  |Human|Pre-calculated Interaction data|human/genesetVHumannet1,2.Rdata|genesetV1,2|  
  |Arabidopsis|Network|networkMousenet.Rdata|network|
  |Arabidopsis|Gene-set|KEGGmouse.Rdata|genesets|
  |Arabidopsis|Pre-calculated Interaction data|human/AragenesetV.RData.Rdata|genesetV|
  |Mouse|Network|networkMousenet.Rdata|network|
  |Mouse|Gene-set|KEGGmouse.Rdata|genesets|
  |Mouse|Pre-calculated Interaction data|genesetVMousenet.Rdata|genesetV|    
  |Yeast|Network|networkYeastnet.Rdata|network|
  |Yeast|Gene-set|KEGGyeast.Rdata|genesets|    
  
</details>

<details>
  <summary>
    <b>Data Formats</b>
  </summary>  
  netGO needs 4 data.<br>
  -	genes: A character vector of input genes (e.g., DE genes).<br>
  -	genesets: A list of gene-sets consisting of groups of genes.<br>
  -	network: A numeric matrix of network data. The network score range is [0,1].<br>
  -	genesetV: A numeric matrix of pre-calculated interaction data between gene and gene-sets.<br>
    The matrix dimension must be [ {# of genes} X {# of gene-sets}]. <br>
    It can be built using BuildGenesetV function with network and gene-set objects as input arguments.

  ```r
  genesetV = BuildGenesetV(network, genesets)
  ```  
</details>

### Functions
<hr>

#### 1. netGO
This function returns a data frame of Gene-sets, and their Q-values derived from netGO+, netGO and Fisher’s exact test.<br>
also P-values available as option.<br>

<details>
 <summary>
  <b>Input arguments</b>
 </summary>
 
  -	genes: A character vector of input genes (e.g., DE genes).<br>
  -	genesets: A list of gene-sets consisting of groups of genes.<br>    
  -	network: A numeric matrix of network data. The network score range is [0,1]. 1 for strong interaction and 0 for no interaction<br>
  
  | |A|B|C|
  |:--:|:--:|:--:|:--:|
  |A|0|0.1|0.76|
  |B|0.1|0|0.324|
  |C|0.76|0.324|0|
  
  -	genesetV: A numeric matrix of pre-calculated interaction data between gene and gene-sets.<br>    
  
  | |Gene-set1|Gene-set2|Gene-set3|
  |:--:|:--:|:--:|:--:|
  |A|0.837|1.647|0.074|
  |B|0|1.75|0.113|
  |C|0.464|0.486|2.442|
    
  -	alpha (optional): A numeric parameter weights how much network score will be effected. <br>
  The value is positive numeric value with > 1 and the default is 20.<br>
  - beta (optional): A numeric parameter balancing the weights between the relative and absolute network scores.<br> 
  The value is between 0 and 1 and the default is 0.5.<br>
  -	nperm(optional): A numeric parameter to tell netGO how many times permutation will be repeated, default is 10000<br>
  - category(optional): A numeric parameter to decide how gene of networks will be categorized with their distribution, default is NULL. and it will divide ~ 2000 genes to each category.<br>
  - pvalue(optional): A boolean parameter whether return Q-values only ( FALSE ) or both P-values and Q-values (TRUE)

</details>
 
After Example run, user may see these logs in R console.

<img src = 'https://user-images.githubusercontent.com/6457691/70369534-7e6c2d80-18fe-11ea-877e-3aa6c79cd4f1.png'></img>

**Notice** member of genes should be gene symbols when using the default STRING and mSigDB data. <br>
Other types of gene names are also available if the corresponding customized data (network and gene-set data) are used.
<hr>

#### 2. netGOVis (for visualization)
This function visualizes the result on the web browser (google chrome is recommended).<br> 
The result graphs and table are downloadable from the web browser.<br>

<details>
 <summary>
  <b>Input arguments</b>
 </summary>

 -	obj: A result data frame derived from <b>netGO</b> function.<br>
 It consists of 6 columns including 1) gene-set name and Q-values evaluated from 2) netGO , 3) netGO+, and 4) Fisher’s exact test.<br>
 and additional 5) Overlap and 6) Network scores. if netGO runned with pvalue option as TRUE, obj will have 3 more columns : netGO, netGO+, FET's pvalue <br>
 
 -	R (optional): Gene-set rank threshold, default is NULL (if 50, Top 50 gene-sets in either method will be shown).<br>

 -	Q (optional): Gene-set Q-value threshold, default is 0.25. (Gene-sets with Q-value <= 0.25 will be used)<br>

 -	genes, genesets, network: same as in the ‘netGO’ function.<br>

</details>

After Example run, user may see these logs in R console.

<img src = 'https://user-images.githubusercontent.com/6457691/70423808-92a45c00-1ab1-11ea-94ac-9fea46a678ca.png'></img>
            
and user's default web browser (<b>we built netGO with chrome</b>) will return this interactive results.

<img src = 'https://user-images.githubusercontent.com/6457691/70423906-c8e1db80-1ab1-11ea-8f1d-6c57454b6062.png'></img>
                        
<hr>

#### 3. Explore with netGO

visualization page of netGO consists of 3 parts : Network, Table, Bubble.

<details>
 <summary>
  <b> Network </b>
 </summary>
 
 - Network panel displays the input genes, selected gene-set, and the network connections between the two.  
 - ![#48dbfb](https://placehold.it/15/48dbfb/000000?text=+) Sky blue nodes represent input genes (e.g., DE genes) 
 - ![#feca57](https://placehold.it/15/feca57/000000?text=+) Yellow nodes represent genes in the selected gene-set 
 - ![#1dd1a1](https://placehold.it/15/1dd1a1/000000?text=+) Green nodes represent the intersection of input genes and the gene-set. 
 - The edge width represents the strength of interaction between two nodes. 
 - ![#feca57](https://placehold.it/15/feca57/000000?text=+) Gene-set Nodes without edges will be discarded.
 - The gene-set can be selected by clicking on the gene-set name in the Table on the right side. 
 - The users can download the graph image as SVG format.

 <img src = 'https://user-images.githubusercontent.com/6457691/70425850-4e1abf80-1ab5-11ea-96a3-84ac3e82d9f9.png'></img>
</details>

<details>
 <summary>
  <b> Table </b>
 </summary>
 
 - Table contains the names of gene-sets and their Q-values ( or P-values ) evaluated from netGO, netGO+ and Fisher’s exact test, respectively. It is downloadable by clicking the ‘Download Table’ button in the upper right of the table 
 <img src ='https://user-images.githubusercontent.com/6457691/70425051-c84a4480-1ab3-11ea-8eb4-4b45385943fe.png'></img>
 
</details>

<details>
 <summary>
  <b> Bubble </b>
 </summary>
 
 - Bubble module plots the bubble chart of significant gene-sets. <br>
 - The overlap (x-axis) and network scores (y-axis) of each significant gene-sets are represented. <br> 
 - The size of bubbles represents the significance level of each gene-set in -log10 scale (Qvalue).
 - Hovering on each bubble will show statistcial values.
 
 <img src ='https://user-images.githubusercontent.com/6457691/70425757-1c095d80-1ab5-11ea-99f6-4198fa48b384.png'></img>
 
</details>

Question / Comment / Suggest : kjh0530@unist.ac.kr <br>
Thanks
