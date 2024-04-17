#' trainTree function
#' 
#' Given a matrix containing gene PatientData levels or a binary matrix indicating mutations, it trains the model.
#'
#' @param PatientData A matrix indicating mutations or a gene expression. Depending on the data, this matrix will
#' be binary(mutational data) or not (expression data).
#' @param PatientSensitivity Drug Response. Higher values of response are considered greater resistance and thus less sensitivity.
#' (Depending on the interpretation of the response the user may have to change the sign of the data)
#' @param minbucket The minimum number of patients on a node in order to split the patients
#'
#'
#' @return  A party of the trained tree with the treatments assigned to each node.
#' 
#' 
#' 
#'
#' @examples 
#' 
#'   \dontrun{
#'   #This first example is a basic example of how to perform trainTree: 
#'   
#'   data(DataODT.rda)
#'   ODTmut <- trainTree(PatientData = mutations_w12,PatientSensitivity=drug_response_w12, minbucket =10)
#'   plot(ODTmut)
#'   
#'   #The next example, is the same as the first one but, 
#'   #using a matrix with gen PatientData data: 
#'   
#'   data(DataODT.rda)
#'   ODTExp <- trainTree(PatientData=expression_w34,PatientSensitivity=drug_response_w34, minbucket = 20)
#'   plot(ODTExp)
#'   
#'   
#'   }
#'
#'
#' @import matrixStats
#' @import partykit
#' @export

trainTree<-function(PatientData, PatientSensitivity, minbucket = 20) {
  exp<-PatientData
  
  if (length(unique(c(unlist(exp)))) == 2){
    PatientData <- PatientData - min(PatientData) + 1L
    mode(PatientData) <- "integer"
    nodes <- growtreeMut(id = 1L, PatientSensitivity, PatientData, minbucket = minbucket)}
  else {
    PatientData<-t(PatientData)
    nodes <- growtreeExp(id = 1L, PatientSensitivity, PatientData, minbucket = minbucket)}
  
  tree <- party(nodes, data = as.data.frame(PatientData))
  tree$node$info <- NULL
  return(tree) }

