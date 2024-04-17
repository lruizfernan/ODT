#' predictTree function
#' 
#' Prediction of test data using ODT
#' 
#'
#' @param tree Trained tree
#' @param PatientSensitivity Drug Response. Higher values of response are considered greater resistance and thus less sensitivity.
#' (Depending on the interpretation of the response the user may have to change the sign of the data)
#' @param PatientData A matrix indicating mutations or a gene expression. Depending on the data, this matrix will
#' be binary(mutational data) or not (expression data).
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' @return A party of the predicted tree with the treatments assigned to each node.
#'
#' @examples 
#' 
#'    \dontrun{
#'   #This first example is a basic 
#'   #example of how to perform predictTree. 
#'   
#'   data(DataODT.rda)
#'   ODTmut <- trainTree(PatientData = mutations_w12,PatientSensitivity=drug_response_w12, minbucket =10)
#'   
#'   
#'   ODT_mutpred<-predictTree(tree=ODTmut, PatientSensitivity=drug_response_w34, PatientData=mutations_w34)
#'   
#'   #The next example, is the same as the first one but,
#'   # using a matrix with genomic expression data: 
#'   
#'   data(DataODT.rda)
#'   ODTExp <- trainTree(PatientData=expression_w34,PatientSensitivity=drug_response_w34, minbucket = 20)
#'   
#'   
#'   ODT_EXPpred<-predictTree(tree=ODTExp, PatientSensitivity=drug_response_w12, PatientData=expression_w12)

#'   }
#'
#' @export

predictTree <- function(tree,PatientSensitivity, PatientData) {
  exp<-PatientData
  if (length(unique(c(unlist(exp)))) == 2){
    PatientData <- PatientData - min(PatientData) + 1L
    mode(PatientData) <- "integer"}
  else{
    PatientData<-t(PatientData)}
  
  treatments <- unlist(nodeapply(tree,  predict.party(tree, as.data.frame(PatientData)), info_node))
  TratamientoTree <- match(treatments, colnames(PatientSensitivity))
  TratamientoTree <- factor(TratamientoTree, levels = 1:ncol(PatientSensitivity))
}