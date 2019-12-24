library(shiny)
library(shinyjs)
library(netGO)
library(DT)

suppressPackageStartupMessages(library(googleVis))

obj <- .GlobalEnv$.obj
genes <- .GlobalEnv$.genes
network <- .GlobalEnv$.network
genesets <- .GlobalEnv$.genesets
Q <- .GlobalEnv$.Q
R <- .GlobalEnv$.R

buildIG <- function(genes, color = "sky") {
  if (length(genes) == 0) {
    return()
  }
  nodes <- list()
  if (color == "sky") {
    color <- "#48DBFB"
  }
  if (color == "yellow") {
    color <- "#FCCE00"
  }
  if (color == "green") {
    color <- "#03CB5D"
  }
  for (i in 1:length(genes)) {
    nodes[[i]] <- buildNode(
      id = genes[i], bgColor = "#FFFFFF",
      borderColor = color, borderWidth = 2,
      fontSize = 10, width = 60, height = 20, opacity = 0
    )
  }
  nodes
}

nodetojs <- function(genes) {
  paste0("cy.nodes('", paste("#", genes, sep = "", collapse = ","), "')")
}

sigIdx <- function(obj, R, Q) {
  pv = obj$`netGO+P`
  pvh = obj$FisherP
  qv <- obj$`netGO+Q`
  qvh <- obj$FisherQ

  names(pv) <- names(pvh) <- obj$`gene-set`
  names(qv) <- names(qvh) <- obj$`gene-set`

  if (!is.null(Q)) {
    idx <- which(qv <= Q | qvh <= Q)
  }
  if (!is.null(R)) {
    idx2 <- which(rank(pv, ties.method = "first") <= R |
                    rank(pvh, ties.method = "first") <= R)
    if(!exists('idx')){return(idx2)}

    idx = union(names(idx),names(idx2))
    idx = sapply(idx, function(i){which(obj$`gene-set`==i)})


  }
  return(idx)
}

buildCol <- function(obj, R, Q) {
  A = B = c()
  if (!is.null(Q)) {
    A <- unname(which(obj$`netGO+Q` <= Q))
    B <- unname(which(obj$FisherQ <= Q))
  }

  if(!is.null(R)){
    if(length(A)){
      A = union(A, unname(which(rank(obj$`netGO+P`, ties.method = "first") <= R)))
      #A = intersect(A, unname(which(rank(obj$`netGO+P`, ties.method = "first") <= R)))
    }
    else{
      A = unname(which(rank(obj$`netGO+P`, ties.method = "first") <= R))
    }
    if(length(B)){
      B = union(B, unname(which(rank(obj$FisherP, ties.method = "first") <= R)))
      #B = intersect(B, unname(which(rank(obj$FisherP, ties.method = "first") <= R)))
    }
    else{
      B = unname(which(rank(obj$FisherP, ties.method = "first") <= R))
    }
  }

  C <- intersect(A, B)
  res <- rep("NONE", length(obj$`netGO+Q`))
  names(res) = obj$`gene-set`
  if (length(A)) {
    res[A] <- "netGO+"
  }
  if (length(B)) {
    res[B] <- "Fisher"
  }
  if (length(C)) {
    res[C] <- "Both"
  }

  res[which(res!='NONE')]
}

getIntersect <- function(gene, geneset) {
  elements <- list()
  gs <- intersect(geneset, rownames(network))
  g <- intersect(gene, rownames(network))

  edges <- list()
  nwe <- c()

  for (i in 1:length(g)) {
    E <- network[g[i], names(which(network[g[i], gs] > 0))]
    if (length(E) > 0) {
      if (length(E) == 1) {
        n <- names(which(network[g[i], gs] > 0))
        edges[[length(edges) + 1]] <-
          buildEdge(source = g[i], target = n, width = (E + .5) * 3, lineColor = "#4B4B4B", opacity = 0)
        nwe <- append(nwe, n, after = length(nwe))
        next
      }
      for (j in 1:length(E)) {
        edges[[length(edges) + 1]] <-
          buildEdge(source = g[i], target = names(E)[j], width = (unname(E)[j] + .5) * 3, lineColor = "#4B4B4B", opacity = 0)
        nwe <- append(nwe, names(E)[j], after = length(nwe))
      }
    }
  }
  elements <- append(elements, edges, after = length(elements))
  list(elements = elements, nwe = unique(nwe))
}

fit <- function(genes, sGs) {
  if (length(intersect(genes, sGs))) {
    runjs(
      paste0(
        "setTimeout(function(){", nodetojs(genes),
        ".layout({name:'random',boundingBox:{x1:0,x2:0.8*cy.width()/3,
             y1:0,y2:cy.height()},fit:false}).run()},3200);
             setTimeout(function(){", nodetojs(sGs),
        ".layout({name:'random',boundingBox:{x1:2.2*cy.width()/3,x2:cy.width(),
             y1:0,y2:cy.height()},fit:false}).run()},3300);
             setTimeout(function(){", nodetojs(intersect(genes, sGs)),
        ".layout({name:'random',boundingBox:{x1:1.2*cy.width()/3,x2:1.8*cy.width()/3,
             y1:0.2*cy.height(),y2:0.8*cy.height()},fit:false}).run()},3400);"
      )
    ) # set butterfly layout
  }
  else {
    runjs(
      paste0(
        "setTimeout(function(){", nodetojs(genes),
        ".layout({name:'random',boundingBox:{x1:0,x2:0.8*cy.width()/3,
             y1:0,y2:cy.height()},fit:false}).run()},3200);
             setTimeout(function(){", nodetojs(sGs),
        ".layout({name:'random',boundingBox:{x1:2.2*cy.width()/3,x2:cy.width(),
             y1:0,y2:cy.height()},fit:false}).run()},3300);"
      )
    ) # set butterfly layout
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
afterCall <- function() {
  runjs("setTimeout(function(){cy.on('click','node',function(e){clicknode(e.target)})},3500)")
  runjs("setTimeout(function(){cy.on('click',function(e){if(e.target===cy){unclick()}})},3500)")
  runjs("setTimeout(function(){cy.panzoom()},3500)")
}

ui <- function() {
  fluidPage(
    useShinyjs(),
    extendShinyjs(script = "shinyjs.js"),
    tags$head(tags$style(type = "text/css", "html,body{width:100%;height:100%;overflow:hidden}")),
    tags$head(tags$style(".row,.col-sm-6,.container-fluid{height:100%}")),
    tags$head(tags$style(".well,#cy{height:95%;margin:1em}")),
    tags$head(tags$script(src = "cytoscape-panzoom.js")),
    tags$head(tags$script(src = "cytoscape-cy-svg-convertor.js")),
    tags$head(tags$script(src = "cytoscape-svg-convertor.js")),
    tags$head(tags$script(src = "svg.min.js")),
    tags$link(rel = "stylesheet", type = "text/css", href = "cytoscape.js-panzoom.css"),
    tags$head(tags$script(src = "additional_script.js")),

    tags$head(tags$style(".dataTables_scrollHeadInner{width:100%;};")),

    div(id = "create", display = "none"), # EMPTY DIV FOR DOWNLOAD SVG

    sidebarLayout(
      position = "right",
      sidebarPanel(
        downloadButton(outputId = "btn3", label = "Download Table", style = "position:absolute;right:3.6em;"),
        div(
          DTOutput(outputId = "table1", height = "100%", width = "100%"),
          style = "height:40%; margin-top:3em;"
        ),
        span(
          htmlOutput("view", inline = TRUE),
          style = "height:40%;left:5%;bottom:5%;position:absolute;width:90%; text-align:center"
        ),
        width = 6
      ),
      mainPanel(

        ShinyCyJSOutput(outputId = "cy", height = "95%"),
        actionButton(
          inputId = "btn4",
          icon = icon("download"),
          label = "Download Graph",
          style = "position:absolute; z-index:9999;right:1em;top:2.5em;"
        ),
        span(textOutput(outputId = "txt1"), style = "font-weight: bold;position: absolute;z-index: 9999;left: 10%;top: 0.5em;"),
        span(
          tags$img(src = "legend.png", style = "width:20em"),
          style = "position: absolute;bottom: 1em;z-index: 9999;",
          span("×", onclick = "this.parentNode.style.display = 'none';", style = "cursor:pointer;width:1em;height:1em; position:absolute; top:1px;right:1px;")
        ),
        width = 6
      )
    )
  )
}

buildelement <- function(genes, sGs) {
  elements <- list()
  elements <- append(elements, buildIG(setdiff(genes, sGs)), after = length(elements)) # sky nodes
  elements <- append(elements, buildIG(intersect(sGs, genes), "#03CB5D"), after = length(elements)) # green nodes
  isobj <- getIntersect(genes, sGs)
  elements <- append(elements, buildIG(isobj$nwe, "#FCCE00"), after = length(elements)) # yellow nodes
  elements <- append(elements, isobj$elements, after = length(elements))
  elements
}

server <- function(input, output, session) {

  si <- sigIdx(obj, R = R, Q = Q)

  myTab <- cbind(names(si), round(cbind(obj$`netGO+Q`, obj$FisherQ)[si, ], 4))
  myTab <- data.frame(myTab, stringsAsFactors = FALSE)
  colnames(myTab) <- c("Gene-set name", "netGO+<br>q-value", "Fisher<br>q-value")

  if(!is.null(obj$netGOQ)){ # netGO and netGO+
    myTab <- cbind(names(si), round(cbind(obj$`netGO+Q`, obj$`netGOQ`, obj$FisherQ)[si, ], 4))
    myTab <- data.frame(myTab, stringsAsFactors = FALSE)
    colnames(myTab) <- c("Gene-set name", "netGO+<br>q-value", "netGO<br>q-value", "Fisher<br>q-value")
  }

  for(i in 2:ncol(myTab)){
    myTab[,i] = as.numeric(myTab[,i])
  }

  rownames(myTab) <- myTab[, 1]
  myTab <- myTab[order(myTab[, 2]), ] # sort by netGO+Q

  sGs <- genesets[[myTab[1, 1]]] # selected Geneset

  elements <- buildelement(genes, sGs)

  output$txt1 <- renderText(paste0("Gene-set :", myTab[1, 1]))
  output$cy <- renderShinyCyJS(shinyCyJS(elements))

  afterCall()
  fit(genes, sGs)

  output$table1 <- renderDT(
    datatable(
      myTab,
      rownames = FALSE,
      extensions = c("Scroller", "Buttons"),
      options = list(
        processing = TRUE,
        order = list(list(1, "asc")),
        deferRender = TRUE,
        scrollY = "34vh",
        scroller = TRUE,
        scrollX = TRUE,
        dom = "ltipr",
        autoWidth = TRUE,
        columnDefs = list(
          list(width = "60px", targets = c(1:(ncol(myTab)-1))),
          #list(width = "60px", targets = 2),
          #list(width = "60px", targets = 3),
          list(width = "100%", targets = 0)
        )
      ),
      selection = "single",
      escape = FALSE
    )
  )

  # network value

  myData <- data.frame(
    name = names(si),
    overlap = unname(sapply(names(si), function(i) {
      obj$OverlapScore[which(obj$`gene-set` == i)]
    })),
    network = unname(sapply(names(si), function(i) {
      obj$NetworkScore[which(obj$`gene-set` == i)]
    })),
    qvalue_log10 = sapply(names(si), function(i) {
      as.numeric(-log10(as.numeric(obj$`netGO+P`[which(obj$`gene-set` == i)])))
    }),
    significant = buildCol(obj, R = R, Q = Q)
  )

  myData <- myData[order(myData[, "significant"], decreasing = TRUE), ]

  output$view <- renderGvis({
    gvisBubbleChart(
      myData,
      idvar = "name", xvar = "overlap",
      yvar = "network", colorvar = "significant", sizevar = "qvalue_log10",
      options = list(
        width = "100%", height = "95%",
        chartArea = "{left:'5%',top:'5%',width:'80%',height:'80%'}",
        # netGO Hyper Both
        # Red #ff7675 Blue #74b9ff Purple #a29bfe
        colors = "['#ff7675', '#74b9ff','#a29bfe']",
        hAxis = "{title : 'Overlap Score'}",
        vAxis = "{title : 'Network Score'}",
        colorAxis = "{legend:{position:'none'}}",
        bubble = '{textStyle:{color : "none" }}', # no label
        explorer = "{}"
      ), chartid = "BubbleChart"
    )
  })

  observeEvent(input$btn2, {
    runjs("cy.$().remove();")
    sIdx <- input$table1_rows_selected
    sGs <<- genesets[[rownames(myTab)[sIdx]]]
    elements <- buildelement(genes, sGs)

    output$cy <- renderShinyCyJS(shinyCyJS(elements))
    fit(genes, sGs)
  })

  # download network svg form
  observeEvent(input$btn4, {
    js$download()
  })

  observeEvent(input$table1_rows_selected, {
    runjs("cy.$().remove();")
    sIdx <- input$table1_rows_selected
    sGs <<- genesets[[rownames(myTab)[sIdx]]]
    elements <- buildelement(genes, sGs)
    output$txt1 <- renderText(paste0("Gene-set :", rownames(myTab)[sIdx]))
    output$cy <- renderShinyCyJS(shinyCyJS(elements))
    fit(genes, sGs)
  })

  output$btn3 <- downloadHandler(
    filename = function() {
      paste0("netGO-", Sys.Date(), ".txt")
    },
    content = function(file) {
      set <- myTab[, 1]
      size <- sapply(1:nrow(myTab), function(i) {
        length(genesets[[ myTab[i, 1] ]])
      })
      pv1 <- myTab[, 2] # netGO+
      pv2 <- myTab[, 3] # Fisher

      # genes
      g <- sapply(1:nrow(myTab), function(i) {
        genesets[[myTab[i, 1]]]
      })
      for (i in 1:length(g)) {
        for (j in 1:length(g[[i]])) {
          if (g[[i]][j] %in% genes) {
            g[[i]][j] <- paste0(g[[i]][j], "(*)")
          }
        }
      }
      g <- sapply(1:nrow(myTab), function(i) {
        paste0(sort(g[[i]]), collapse = ",")
      })

      text <- c("Gene-set\tSize\tnetGO+-Qvalue\tFisher-Qvalue\tGenes")

      if(ncol(myTab)==4){# netGO exist
        pv3 = myTab[,4] # Fisher
        text <- c("Gene-set\tSize\tnetGO+-Qvalue\tnetGO-Qvalue\tFisher-Qvalue\tGenes")
        for (i in 1:(length(set))) {
          text[i + 1] <- paste(set[i], size[i], pv1[i], pv2[i], pv3[i], g[[i]], sep = "\t")
        }
      }
      else{
        for (i in 1:(length(set))) {
          text[i + 1] <- paste(set[i], size[i], pv1[i], pv2[i], g[[i]], sep = "\t")
        }
      }

      writeLines(text, file)
    }
  )
}
shinyApp(ui, server)
