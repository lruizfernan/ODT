#' ODT Internal Functions
#'
#' Internal functions used by trainTree in the different
#' steps of the algorithm
#'
#' @keywords internal
#' @name InternalFunctions
#' @return Internal outputs
#'
#'
NULL


#' @rdname InternalFunctions
findsplitExp <- function(PatientSensitivity, PatientData, minimum = 1, weights=NULL, verbose = F) {
  # PatientSensitivity: IC50
  # X: mutation matrix
  if (verbose) print("Split!")
  if (!is.null(weights)) {
    PatientSensitivity <- PatientSensitivity[weights == 1,]
    PatientData <- PatientData[weights == 1,]
  }
  if (is.matrix(PatientSensitivity)==TRUE){
    if (nrow(PatientSensitivity)<=minimum) return(NULL)}
  else {return(NULL)}
  
  Gene <- which.min(apply(PatientData, 2, getsumic50v3, PatientSensitivity)) # Find optimal split
  names(Gene) <- colnames(PatientData[,Gene])
  Output <- getSplit(Gene, PatientData, PatientSensitivity)
  return(partysplit(varid = as.integer(Gene), 
                    breaks = Output$expressionSplit,
                    index = c(1L,2L),
                    info = list(treatments = c(names(Output$T1), names(Output$T2)))))
}

getTreatment <- function(PatientSensitivity, weights) {
  # PatientSensitivity: IC50
  PatientSensitivity <- PatientSensitivity[weights == 1,,drop=F]
  sumtreat <- colSums2(PatientSensitivity)
  biomk <- which.min(sumtreat)
  names(biomk) <- colnames(PatientSensitivity)[biomk]
  return(biomk)
}


growtreeExp <- function(id = 1L, PatientSensitivity, PatientData, minbucket = 10, weights = NULL) {
  if (is.null(weights))  weights <- rep(1L, nrow(PatientData))
  
  ## for less than minbucket observations stop here
  if (sum(weights) < minbucket) return(partynode(id = id, info = names(getTreatment(PatientSensitivity, weights))))
  
  ## find best split
  sp <- findsplitExp(PatientSensitivity, PatientData, minimum = minbucket, weights)
  datadf <- as.data.frame(PatientData)
  ## no split found, stop here
  if (is.null(sp)) return(partynode(id = id, info = names(getTreatment(PatientSensitivity, weights))))
  
  ## actually split the data
  kidids <- kidids_split(sp, data = datadf)
  if(length(unique(kidids[weights==1]))==1) 
    return(partynode(id = id, info = names(getTreatment(PatientSensitivity, weights))))
  if(min(table(kidids[weights==1])) < minbucket) 
    return(partynode(id = id, info = names(getTreatment(PatientSensitivity, weights))))
  
  ## set up all daugther nodes
  kids <- vector(mode = "list", length = max(kidids, na.rm = TRUE))
  for (kidid in 1:length(kids)) {
    ## select observations for current node
    w <- weights
    w[kidids != kidid] <- 0
    ## get next node id
    if (kidid > 1) {
      myid <- max(nodeids(kids[[kidid - 1]]))
    } else {
      myid <- id
    }
    ## start recursion on this daugther node
    kids[[kidid]] <- growtreeExp(id = as.integer(myid + 1), PatientSensitivity, 
                                 PatientData, weights =w, minbucket = minbucket)
  }
  
  ## return nodes
  return(partynode(id = as.integer(id), split = sp, kids = kids,
                   info = ""))
}

getSplit <- function(gene, PatientData, PatientSensitivity) {
  I <- order(PatientData[,gene])
  Id <- order(PatientData[,gene], decreasing = T)
  
  M1 <- rbind(0,colCumsums(PatientSensitivity[I,]))
  M2 <- rbind(0,colCumsums(PatientSensitivity[Id,]))
  M2 <- M2[c(nrow(M1):1),,drop=F]
  
  
  # plot(rowMins(M1) + rowMins(M2))
  th <- which.min(rowMins(M1) + rowMins(M2))
  sumic50 <- min(rowMins(M1) + rowMins(M2))
  T1 = max.col(-M1)[th]
  names(T1) = colnames(PatientSensitivity)[T1]
  T2 = max.col(-M2)[th]
  names(T2) = colnames(PatientSensitivity)[T2]
  
  sumic50
  return(list(sumic50 = sumic50, 
              T1 = T1,
              T2 = T2,
              split = th,
              expressionSplit = PatientData[I[th],gene]-1e-6))
}

getsumic50v3 <- function(genePatientResponse, PatientSensitivity) {
  I <- order(genePatientResponse)
  Id <- I[length(I):1]
  
  M1 <- colCumsums(PatientSensitivity, rows=I)
  M2 <- colCumsums(PatientSensitivity, rows=Id)
  dummy <- c(0,rowMins(M1)) + c(rowMins(M2, rows=(nrow(M2):1)),0)
  sumic50 <- min(dummy)
  return(sumic50)
}

# library(partykit)

trainTreeExp <- function(PatientData, PatientSensitivity, minbucket = 20) {
  nodes <- growtreeExp(id = 1L, PatientSensitivity, PatientData, minbucket = minbucket)
  tree <- party(nodes, data = as.data.frame(PatientData))
  tree$node$info <- NULL
  return(tree)
}

predictTreeExp <- function(tree,PatientSensitivity, PatientData) {
  treatments <- unlist(nodeapply(tree, 
                                 predict.party(tree, as.data.frame(PatientData))
                                 , info_node))
  TratamientoTree <- match(treatments, colnames(PatientSensitivity))
  TratamientoTree <- factor(TratamientoTree, levels = 1:ncol(PatientSensitivity))
}

findsplitMut <- function(PatientSensitivity, X, minimum = 1, weights) {
  # PatientSensitivity: IC50
  # X: mutation matrix
  PatientSensitivity <- PatientSensitivity[weights == 1,]
  X <- X[weights == 1,]
  
  if (is.matrix(PatientSensitivity)==TRUE){
    if (nrow(PatientSensitivity)<=minimum) return(NULL)}
  else {return(NULL)}
  
  mut <- max(X) # Mutation is the largest value of X
  WT <- min(X) # Wildtype is the smalles value of X
  
  tA <- t(PatientSensitivity) %*% (X==mut) # treatment with mutation
  tB <- t(PatientSensitivity) %*% (X==WT)  # Treatment for WT
  
  #   tA Sum of the IC50 of the mutated samples treated with each treatment
  #   tB Sum of the IC50 of the wild type samples treated with each treatment
  
  
  # The biomarker mutation
  optimal <- colMins(tA) + colMins(tB) # Optimal selection of treatments for each mutation
  biomk <- which.min(optimal) # Optimal mutation
  names(biomk) <- colnames(tA)[biomk]
  
  # If both mutated and WT samples return the same treatment, don't split
  if (length(unique(X[,biomk])) == 1) return(NULL)  ##Parece que sirve para asegurarnos
  #que para esa mutacion hay pacientes que la tienen y otros que no. Si todos la tuviesen
  # o todos no la tuviesen no tendria sentido separar en esa mutacion.
  
  # The treatments
  T1Mut <- which.min(tA[,biomk])
  T1WT <- which.min(tB[,biomk])
  
  # If both mutated and WT samples return the same treatment, don't split. Probably already done.
  #Comprueba que no es el mismo farmaco para los que tienen y no tienen mutacion.
  if (T1Mut == T1WT) return(NULL)
  
  index = c(1L,2L)
  names(index) <- c(paste0(names(biomk),"mut"), paste0(names(biomk),"WT"))
  return(partysplit(varid = as.integer(biomk), 
                    index = index,
                    info = list(treatments = c(names(T1Mut), names(T1WT)))))
}

growtreeMut <- function(id = 1L, PatientSensitivity, PatientData, minbucket = 10, weights = NULL, findsplit = findsplitMut) { # nolint: line_length_linter.
  if (is.null(weights))  weights <- rep(1L, nrow(PatientData))
  
  # Parameters:
  # id: unique Id of the node
  # response: IC50, or in general, the funciton to minminze.
  # data: PatientData, mutations, or in general, the data to do the split
  # minbucket: minimum number of samples in child node.
  # weights: don't set it. Variable to state which are the samples under study.
  
  ## for less than "minbucket"  observations stop here
  if (sum(weights) < minbucket) return(partynode(id = id, info = names(getTreatment(PatientSensitivity, weights)))) # nolint: line_length_linter.
  
  ## find best split
  sp <- findsplitMut(PatientSensitivity, PatientData, minimum = minbucket, weights)
  
  ## no split found, stop here
  if (is.null(sp)) return(partynode(id = id, info = names(getTreatment(PatientSensitivity, weights))))
  
  
  ## actually split the data
  datadf <- as.data.frame(PatientData)
  kidids <- kidids_split(sp, data = datadf)
  
  ## If only one kid, return. Solo un tipo, return.
  if(length(unique(kidids))==1) 
    return(partynode(id = id, info = names(getTreatment(PatientSensitivity, weights))))
  
  ## If any of the splits smaller than minbucket, return
  if(min(table(kidids[weights==1,drop=F])) < minbucket) 
    return(partynode(id = id, info = names(getTreatment(PatientSensitivity, weights))))
  
  ## set up all daugther nodes
  kids <- vector(mode = "list", length = max(kidids, na.rm = TRUE)) #na.rm=TRUE quita los missing values
  #para hacer el calculo del max
  for (kidid in 1:length(kids)) {
    ## select observations for current node
    w <- weights
    w[kidids != kidid] <- 0
    ## get next node id
    if (kidid > 1) {
      myid <- max(nodeids(kids[[kidid - 1]]))
    } else {
      myid <- id
    }
    ## start recursion on this daugther node
    kids[[kidid]] <- growtreeMut(id = as.integer(myid + 1), PatientSensitivity, 
                                 PatientData, weights =w, minbucket = minbucket, findsplit = findsplitMut)
  }
  
  ## return nodes
  return(partynode(id = as.integer(id), split = sp, kids = kids,
                   info = ""))
}

trainTreeMut <- function(PatientData, PatientSensitivity, minbucket = 20) {
  PatientData <- PatientData - min(PatientData) + 1L
  mode(PatientData) <- "integer"
  nodes <- growtreeMut(id = 1L, PatientSensitivity, PatientData, minbucket = minbucket)
  tree <- party(nodes, data = as.data.frame(PatientData))
  tree$node$info <- NULL
  return(tree)
}

predictTreeMut <- function(tree,PatientSensitivity, PatientData) {
  PatientData <- PatientData - min(PatientData) + 1L
  mode(PatientData) <- "integer"
  predictTree(tree,PatientSensitivity, PatientData)
}


