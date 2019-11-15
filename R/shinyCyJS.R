library(htmlwidgets)

#' @export
shinyCyJS <- function(elements = list(), options = list(), layout = list(name = "cose"),
                      width = NULL, height = NULL, elementId = NULL) {
  input <- list(elements = elements, options = options, layout = layout)
  htmlwidgets::createWidget(
    name = "shinyCyJS",
    input,
    width = width,
    height = height,
    package = "netGO",
    elementId = elementId
  )
}

#' @export
buildIOptions <- function(...) {
  # touchTapThreshold, desktopTapthreshold not available
  options <- as.list(match.call())[-1] # remove function name
  options <- options[
    c(
      "minZoom", "maxZoom", "zoomingEnabled",
      "userZoomingEnabled", "panningEnabled", "userPanningEnabled",
      "boxSelectionEnabled", "selectionType", "autolock",
      "autoungrabify", "autounselectify"
    )
  ]
  options[which(!sapply(options, is.null))] # if option not given : remove
}

#' @export
buildROptions <- function(...) {
  options <- as.list(match.call())[-1] # remove function name
  options <- options[
    c(
      "headless", "styleEnabled", "hideEdgesOnViewport",
      "hideLabelsOnViewport", "textureOnViewport", "motionBlur",
      "motionBlurOpacity", "wheelSensitivity", "pixelRatio"
    )
  ]
  options[which(!sapply(options, is.null))] # if option not given : remove
}

#' @export
ShinyCyJSOutput <- function(outputId, width = "100%", height = "400px") {
  htmlwidgets::shinyWidgetOutput(outputId, "shinyCyJS", width, height, package = "shinyCyJS")
}

#' @export
renderShinyCyJS <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) {
    expr <- substitute(expr)
  } # force quoted
  htmlwidgets::shinyRenderWidget(expr, ShinyCyJSOutput, env, quoted = TRUE)
}

#' @export
buildNode <- function(
                      id, width = 15, height = 15, shape = "ellipse", bgColor = "#48DBFB",
                      bgOpacity = 1, bgFill = "solid", borderWidth = 0, borderStyle = "solid",
                      borderColor = "#8395a7", borderOpacity = 1, isParent = FALSE, labelColor = "#8395a7",
                      textOpacity = 1, fontSize = 16, textOutlineColor = "#222f3e", textOutlineOpacity = 1,
                      textOutlineWidth = 0, textbgColor = "#FFF", textbgOpacity = 0, textBorderColor = "#222f3e",
                      textBorderOpacity = 0, textBorderWidth = 0, parent = NULL, opacity = 1) {
  # megaman SKY #48DBFB
  # storm petrel #8395A7
  # impreial primer #222F3E

  l <- list(group = "nodes")

  options <- list(
    width = width, height = height, label = id, id = id, shape = shape,
    bgColor = bgColor, bgOpacity = bgOpacity, bgFill = bgFill,
    borderWidth = borderWidth, borderStyle = borderStyle, borderColor = borderColor,
    borderOpacity = borderOpacity, labelColor = labelColor, textOpacity = textOpacity,
    fontSize = fontSize, textOutlineColor = textOutlineColor, textOutlineOpacity = textOutlineOpacity,
    textOutlineWidth = textOutlineWidth, textbgColor = textbgColor, textbgOpacity = textbgOpacity,
    textBorderColor = textBorderColor, textBorderOpacity = textBorderOpacity, textBorderWidth = textBorderWidth, parent = parent,
    opacity = opacity
  )

  if (isParent) { # parent node
    options$bgColor <- "#c8d6e5" # light blue ballerina
    options$bgOpacity <- "0.5"
  }

  l$data <- options

  l
}

#' @export
buildEdge <- function(source, target, width = 3, curveStyle = "haystack",
                      lineColor = "#FECA57", lineStyle = "solid", sourceArrowColor = "#feca57",
                      targetArrowColor = "#feca57", sourceArrowShape = "none", targetArrowShape = "none",
                      opacity = 1) {
  # casandora YELLOW #FECA57
  l <- list()
  l$group <- "edges"
  options <- list(
    source = source, target = target, width = width,
    curveStyle = curveStyle, lineColor = lineColor, lineStyle = lineStyle,
    sourceArrowColor = sourceArrowColor, targetArrowColor = targetArrowColor,
    sourceArrowShape = sourceArrowShape, targetArrowShape = targetArrowShape,
    opacity = opacity
  )
  l$data <- options
  l
}
