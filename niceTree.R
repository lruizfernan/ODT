#' niceTree function
#' 
#' A graphical display of the tree. It can also be saved as an image in the selected directory.
#' 
#' @param tree A party of the trained tree with the treatments assigned to each node.
#' @param folder Directory to save the image
#' @param colors Selection of colors for the boxes
#' @param fontname The font name
#' @param fontstyle The font style
#' @param shape format of the boxes of the different genes
#' @param output_format image format
#'
#' @details 
#' \itemize{
#'  \item The user has already a defined style for the plot, the parameters
#'        are previously set if the parameters are not modified when calling
#'        niceTree
#' }
#'  
#'  
#'   
#'
#' @return (Invisibly) returns a list. The representation of the tree in the command window and the plot of the tree
#'
#' @examples 
#' 
#' \dontrun{
#'   #This first example is a basic example of how to perform niceTree: 
#'   
#'   data(DataODT.rda)
#'   ODTmut <- trainTree(PatientData = mutations_w12,
#'   PatientSensitivity=drug_response_w12, minbucket =10)
#'   niceTree(ODTMut,folder="")
#' 
#'   
#'
#'   #The next example, is the same as the first
#'   # one but, plotting the tree trained for gene expressions. 
#'   
#'   data("DATA_LOADED.RData")
#'   ODTExp <- trainTree(PatientData=expression_w34,
#'   PatientSensitivity=drug_response_w34, minbucket = 20)
#'   niceTree(ODTExp,folder="")

#'  }
#'  
#' @importFrom data.tree as.Node
#' @export


niceTree <- function(tree, folder = NULL,
                     colors = c("", "#367592", "#39A7AE", "#96D6B6",
                     "#FDE5B0", "#F3908B", "#E36192",
                                "#8E4884", "#A83333"),
                     fontname = "Roboto", fontstyle = "plain",
                     shape = "diamond", output_format = "png") {
  # library(DiagrammeR)
  # library(DiagrammeRsvg)
  # library(svglite)
  # library(ggparty)
  # library(data.tree)
  
  
  treeNode <- as.Node(tree)
  
  levels <- max(treeNode$Get(function(x) c(level = x$level)))
  # bucle
  for (i in 1:levels) {
    color_node <- colors[(i %% length(colors)) + 1]
    nodes_level <- Traverse(treeNode, filterFun = function(x) x$level == i)
    if (i > 1) {
      for (j in 1:length(nodes_level)) {
        if (nodes_level[[j]][["splitLevel"]] == "<= 1") {
          nodes_level[[j]][["splitLevel"]] <- 'WT'
        } else if (nodes_level[[j]][["splitLevel"]] == "> 1") {
          nodes_level[[j]][["splitLevel"]] <- 'Mut'}}}
    
    Do(nodes_level, function(node) {
      SetNodeStyle(node,
                   label = function(node) paste0(node$splitname),
                   tooltip = function(node) paste0(nrow(node$data), " observations"),
                   shape = shape,
                   style = "filled",
                   color = color_node,
                   fillcolor = paste(color_node, "88", sep = ""),
                   fontcolor = "black",
                   fontname = fontname,
                   fontstyle = fontstyle)})}
  
  SetEdgeStyle(treeNode,
               arrowhead = "none",
               label = function(node) node$splitLevel,
               fontname = fontname,
               fontstyle = fontstyle,
               penwidth = function(node) 12 * nrow(node$data) / nrow(node$root$data),
  )
  
  Do(treeNode$leaves, function(node) SetNodeStyle(node, shape = "box", fontname = fontname, fontstyle = fontstyle))
  
  if (!is.null(folder)) {
    # Save plot
    tmp <- DiagrammeRsvg::export_svg(plot(treeNode))
    
    if (output_format %in% c("png", "jpg")) {
      # PatientSensitivityonvert SVG to PNG using magick
      magick::image_write(magick::image_read_svg(rawToChar(charToRaw(tmp))), 
                          path = paste(folder, "/tree_plot.", output_format, sep = ""))
    } else if (output_format == "svg") {
      # Save SVG directly
      writeBin(charToRaw(tmp), paste(folder, "/tree_plot.svg", sep = ""))
    } else if (output_format == "pdf") {
      # Save PDF directly
      pdf(paste(folder, "/tree_plot.pdf", sep = ""), width = 12, height = 12)
      plot(treeNode)
      dev.off()
    } else {
      warning("Unsupported output format. Defaulting to PNG.")
      # PatientSensitivityonvert SVG to PNG using magick
      magick::image_write(magick::image_read_svg(rawToChar(charToRaw(tmp))), 
                          path = paste(folder, "/tree_plot.png", sep = ""))
    }
  }
  return(list(Tree = tree, Plot = plot(treeNode)))
}

