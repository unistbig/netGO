#' @import shiny
#' @import shinyjs
#' @import netGO
#' @import DT

# Data prepare


# netGOVis

library(shiny)
library(shinyjs)
library(netGO)
library(DT)
suppressPackageStartupMessages(library(googleVis))

obj = .GlobalEnv$.obj
genes = .GlobalEnv$.genes
network = .GlobalEnv$.network
genesets = .GlobalEnv$.genesets
Q = .GlobalEnv$.Q
R = .GlobalEnv$.R

buildIG = function(genes, color = 'sky'){
  if(length(genes)==0){return();}
  nodes = list()
  if(color == 'sky'){color = '#48DBFB'}
  if(color == 'yellow'){color = '#FCCE00'}
  if(color == 'green'){color = '#03CB5D'}
  for(i in 1:length(genes)){
    nodes[[i]] = buildNode(
      id = genes[i], bgColor = '#FFFFFF',
      borderColor = color, borderWidth = 2,
      fontSize = 10,width = 60, height = 20, opacity = 0)
  }
  nodes
}

nodetojs = function(genes){paste0("cy.nodes('", paste('#',genes,sep='',collapse = ','),"')")}

sigIdx = function(obj, R, Q){
  pv = obj$netGOP
  pvh = obj$FisherP
  names(pv) = names(pvh) = obj[,1]
  if(!is.null(Q)){
    idx = which(p.adjust(pv,'fdr') <= Q | p.adjust(pvh,'fdr') <= Q)
    return(idx)
  }
  else{
    idx = which(rank(pv, ties.method = 'first') <= R | rank(pvh, ties.method = 'first') <= R)
  }

  return(idx)
}

buildCol = function(obj, R, Q){
  if(!is.null(Q)){
    A = unname(which(p.adjust(obj$netGOP,'fdr')<=Q))
    B = unname(which(p.adjust(obj$FisherP,'fdr')<=Q))
  }
  else{
    A = unname(which(rank(obj$netGOP, ties.method = 'first')<=R))
    B = unname(which(rank(obj$FisherP, ties.method = 'first')<=R))
  }
  C = intersect(A,B)
  res = rep('NONE', length(obj$netGOP))
  if(length(A)){res[A] = 'netGO'}
  if(length(B)){res[B] = 'Fisher'}
  if(length(C)){res[C] = 'Both'}
  res
}

getIntersect = function(gene,geneset){
  elements = list()
  gs = intersect(geneset,rownames(network))
  g = intersect(gene,rownames(network))

  edges = list()
  nwe = c()

  for(i in 1:length(g)){
    E = network[g[i],names(which(network[g[i],gs]>0))]
    if(length(E)>0){
      if(length(E)==1){
        n = names(which(network[g[i],gs]>0))
        edges[[length(edges)+1]] =
          buildEdge(source = g[i], target = n, width = (E+.5)*3, lineColor = '#4B4B4B', opacity = 0)
        nwe = append(nwe, n, after = length(nwe))
        next
      }
      for(j in 1:length(E)){
        edges[[length(edges)+1]] =
          buildEdge(source = g[i], target = names(E)[j], width = (unname(E)[j]+.5)*3, lineColor = '#4B4B4B', opacity = 0)
        nwe = append(nwe, names(E)[j], after = length(nwe))
      }
    }
  }
  elements = append(elements, edges, after = length(elements))
  list(elements = elements, nwe = unique(nwe))
}

fit = function(genes, sGs){

  if(length(intersect(genes,sGs))){
    runjs(
      paste0('setTimeout(function(){',nodetojs(genes),
             ".layout({name:'random',boundingBox:{x1:0,x2:0.8*cy.width()/3,
             y1:0,y2:cy.height()},fit:false}).run()},3200);
             setTimeout(function(){",nodetojs(sGs),
             ".layout({name:'random',boundingBox:{x1:2.2*cy.width()/3,x2:cy.width(),
             y1:0,y2:cy.height()},fit:false}).run()},3300);
             setTimeout(function(){",nodetojs(intersect(genes,sGs)),
             ".layout({name:'random',boundingBox:{x1:1.2*cy.width()/3,x2:1.8*cy.width()/3,
             y1:0.2*cy.height(),y2:0.8*cy.height()},fit:false}).run()},3400);"
             )) # set butterfly layout
  }
  else{
    runjs(
      paste0('setTimeout(function(){',nodetojs(genes),
             ".layout({name:'random',boundingBox:{x1:0,x2:0.8*cy.width()/3,
             y1:0,y2:cy.height()},fit:false}).run()},3200);
             setTimeout(function(){",nodetojs(sGs),
             ".layout({name:'random',boundingBox:{x1:2.2*cy.width()/3,x2:cy.width(),
             y1:0,y2:cy.height()},fit:false}).run()},3300);"

      )) # set butterfly layout
  }

  runjs( # aftercall
    "setTimeout(function(){
    cy.nodes().style('text-valign','center');
    cy.fit();
    cy.nodes().style('opacity',1);
    cy.edges().style('opacity',0.3);
           },3500)"
  ) # show graphs
}
afterCall = function(){
  runjs("setTimeout(function(){cy.on('click','node',function(e){clicknode(e.target)})},3500)")
  runjs("setTimeout(function(){cy.on('click',function(e){if(e.target===cy){unclick()}})},3500)")
  runjs('setTimeout(function(){cy.panzoom()},3500)')
}

ui = function(){
  fluidPage(
    useShinyjs(),
    extendShinyjs(script = 'shinyjs.js'),
    tags$head(tags$style(type="text/css","html,body{width:100%;height:100%;overflow:hidden}")),
    tags$head(tags$style('.row,.col-sm-6,.container-fluid{height:100%}')),
    tags$head(tags$style('.well,#cy{height:95%;margin:1em}')),
    tags$head(tags$script(src="cytoscape-panzoom.js")),
    tags$head(tags$script(src="cytoscape-cy-svg-convertor.js")),
    tags$head(tags$script(src="cytoscape-svg-convertor.js")),
    tags$head(tags$script(src="svg.min.js")),
    tags$link(rel = "stylesheet", type = "text/css", href = "cytoscape.js-panzoom.css"),
    tags$head(tags$script(src="additional_script.js")),
    tags$head(tags$style('.dataTables_scrollHeadInner{width:100%;};')), # margin:1em;margin-top:4em

    tags$head(tags$style('#view{height:45%;left:5%;bottom:5%;position:absolute;width:90%; text-align:center};')), # margin:1em;margin-top:4em
    div(id='create', display='none'), # EMPTY DIV FOR DOWNLOAD SVG

    #tags$head(tags$script(src="cytoscape-svg-convertor.js")),
    sidebarLayout(
      position = 'right',
      sidebarPanel(
        downloadButton(outputId = "btn3", label = "Download Table",style='position:absolute;right:3.6em;'),
        div(
          DTOutput(outputId='table1',height = '25%', width ='100%'),
          style='height:30%; position:absolute;top:6%;width:90%;'),
        #actionButton(inputId='btn2',label='Plot Gene-set network',style='position:absolute;right:18em;margin-bottom:1em;'),

        htmlOutput("view"),

        width = 6
      ),
      mainPanel(
        #actionButton(
        #inputId = 'btn',label = 'Export',
        #style = 'position:absolute;top:1em;right:1em; z-index:9999;'
        #),
        ShinyCyJSOutput(outputId = 'cy', height = '95%'),
        actionButton(
          inputId = 'btn4',
          icon = icon('download'),
          label = 'Download Graph',
          style='position:absolute; z-index:9999;right:1em;top:2.5em;'
        ),
        span(textOutput(outputId = 'txt1'), style='font-weight: bold;position: absolute;z-index: 9999;left: 10%;top: 0.5em;'),
        span(tags$img(src = 'legend.png', style='width:20em'), style ='position: absolute;bottom: 1em;z-index: 9999;'),
        width = 6
      )
    )
  )
}

buildelement = function(genes, sGs){
  elements = list()
  elements = append(elements, buildIG(setdiff(genes,sGs) ), after = length(elements))
  elements = append(elements, buildIG(intersect(sGs, genes), '#03CB5D'), after = length(elements))
  isobj = getIntersect(genes,sGs)
  elements = append(elements, buildIG(isobj$nwe, '#FCCE00'), after = length(elements))
  elements = append(elements, isobj$elements, after = length(elements))
  elements
}

server = function(input,output,session){
  # build example network

  si = sigIdx(obj,R = R,Q = Q)

  myTab = cbind(names(si),round(cbind(p.adjust(obj$netGOP,'fdr'),p.adjust(obj$FisherP,'fdr'))[si,],4))

  myTab = data.frame(myTab, stringsAsFactors = FALSE)

  myTab[,2] = as.numeric(myTab[,2])
  myTab[,3] = as.numeric(myTab[,3])
  rownames(myTab) = myTab[,1]
  colnames(myTab) = c("Gene-set name","netGO<br>q-value","Fisher's exact test<br>q-value")
  myTab = myTab[order(myTab[,2]),]
  sGs = genesets[[myTab[1,1]]]

  elements = buildelement(genes, sGs)
  output$txt1 = renderText(paste0('Gene-set :',myTab[1,1]) )
  output$cy = renderShinyCyJS(shinyCyJS(elements))

  afterCall()
  fit(genes, sGs)

  output$table1 = renderDT(
    datatable(
    myTab,
    rownames = FALSE,
    extensions = c('Scroller', 'Buttons'),#,'Responsive'),
    options = list(
      processing = TRUE,
      order = list(list(1,'asc')),
      deferRender = TRUE,
      scrollY = "20vh",
      scroller = TRUE,
      scrollX = TRUE,
      dom = 'ltipr'
      ,autoWidth = FALSE
      ,columnDefs = list(
        #list(width ='10em', targets = 0),
        list(width ='200px', targets = 1),
        list(width ='200px', targets = 2)
      )
    ),
    selection = 'single',
    escape = FALSE
    )#,server = FALSE
  )

  myData = data.frame(
    name = names(si),
    network = unname(sapply(si, function(i){sum(network[intersect(rownames(network),genes),intersect(rownames(network), genesets[[i]])])/ length(genesets[[i]]) })),
    overlap = unname(sapply(si, function(i){length(intersect(genes,genesets[[i]]))/ length(genesets[[i]])})),
    pvalue_log10 = sapply(names(si), function(i){ as.numeric(-log10(as.numeric(myTab[i,2]) )) }),
    significant = buildCol(obj, R = R, Q = Q)[si]
  )

  myData = myData[order(myData[,'significant'], decreasing = TRUE),]


  output$view = renderGvis({
    gvisBubbleChart(
      myData,idvar = 'name',xvar = 'overlap',
      yvar = 'network', colorvar = 'significant',sizevar ='pvalue_log10',
      options = list(
        width = '100%', height = '100%',
        chartArea = "{left:'5%',top:'5%',width:'80%',height:'80%'}",
        # netGO Hyper Both
        # Red #ff7675 Blue #74b9ff Purple #a29bfe
        colors= "['#ff7675', '#74b9ff','#a29bfe']",
        hAxis = "{title : 'Overlap Score'}",
        vAxis = "{title : 'Network Score'}",
        colorAxis = "{legend:{position:'none'}}",
        bubble = '{textStyle:{color : "none" }}',# no label
        explorer ='{}'), chartid = 'BubbleChart'
    )
  })

  observeEvent(input$btn2, {
    runjs('cy.$().remove();')
    sIdx = input$table1_rows_selected
    sGs <<- genesets[[rownames(myTab)[sIdx]]]
    elements = buildelement(genes, sGs)

    output$cy = renderShinyCyJS(shinyCyJS(elements))
    fit(genes, sGs)
  })

  # download network svg form
  observeEvent(input$btn4,{ js$download() })

  observeEvent(input$table1_rows_selected, {
    runjs('cy.$().remove();')
    sIdx = input$table1_rows_selected
    sGs <<- genesets[[rownames(myTab)[sIdx]]]
    elements = buildelement(genes, sGs)
    output$txt1 = renderText(paste0('Gene-set :',rownames(myTab)[sIdx]))
    output$cy = renderShinyCyJS(shinyCyJS(elements))
    fit(genes, sGs)
  })

  output$btn3 = downloadHandler(
    filename = function(){paste0("netGO-", Sys.Date(), ".txt")},
    content = function(file) {
      set  = myTab[,1]
      size = sapply(1:nrow(myTab), function(i){length( genesets[[ myTab[i,1] ]] ) })
      pv1 = myTab[,2]
      pv2 = myTab[,3]
      g = sapply(1:nrow(myTab), function(i){genesets[[myTab[i,1]]] })
      for(i in 1:length(g)){
        for(j in 1:length(g[[i]])){
          if(g[[i]][j] %in% genes){g[[i]][j] = paste0(g[[i]][j],'(*)')}
        }
      }
      g = sapply(1:nrow(myTab), function(i){paste0(sort(g[[i]]), collapse = ',') })
      text =c('Gene-set\tSize\tnetGO-Qvalue\tFisher-Qvalue\tGenes')
      for(i in 1:(length(set))){
        text[i+1] = paste(set[i], size[i], pv1[i],pv2[i],g[[i]], sep = '\t')
      }
      writeLines(text, file)
    }
  )


}
shinyApp(ui,server)
