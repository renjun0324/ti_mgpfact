
#' %nin%
#'
#' @export
'%nin%' <- function(x,y)!('%in%'(x,y))

#' LogNormalize
#'
#'
LogNormalize <- function(
    data,
    scale.factor = 1e4,
    margin = 2L,
    verbose = TRUE,
    ...
) {
  ncells <- dim(x = data)[margin]
  if (isTRUE(x = verbose)) {
    pb <- txtProgressBar(file = stderr(), style = 3)
  }
  for (i in seq_len(length.out = ncells)) {
    x <- if (margin == 1L) {
      data[i, ]
    } else {
      data[, i]
    }
    xnorm <- log1p(x = x / sum(x) * scale.factor)
    if (margin == 1L) {
      data[i, ] <- xnorm
    } else {
      data[, i] <- xnorm
    }
    if (isTRUE(x = verbose)) {
      setTxtProgressBar(pb = pb, value = i / ncells)
    }
  }
  return(data)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                             
#                                 object
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The DimReduc-class
#'
#' @slot PCA pca result
#' @slot tSNE tsne result
#' @slot DM diffusion map result
#' @slot UMAP umap result
#'
#' @name DimReduc-class
#' @rdname DimReduc-class
#' @exportClass DimReduc
#'·
DimReducs <- setClass(
  Class = 'DimReducs',
  slots = c(
    PCA = 'ANY',
    tSNE = 'ANY',
    DM = 'ANY',
    UMAP = 'ANY'
  )
)

#' The MGPfact-class
#'
#' @slot assay
#' @slot OptimResult
#' @slot MURP
#' @slot DimReducs
#' @slot MetaData
#' @slot Date
#'
#' @name MGPfact-class
#' @rdname MGPfact-class
#' @exportClass MGPfact
#'
MGPfact <- setClass(
  Class = 'MGPfact',
  slots = c(
    assay = 'list',
    OptimResult = 'ANY',
    MURP = 'ANY',
    DimReducs = 'DimReducs',
    MetaData = 'data.frame',
    Settings = 'ANY',
    Date = 'ANY',
    BeautGene = 'list',
    Tree = 'list',
    GPR = 'list',
    Command = 'list'
  )
)

#' The PseOptimResult-class
#'
#' @slot sample_value sampling results for each iteration
#' @slot dimnames names of sampling value
#' @slot burnin burning period
#' @slot chains chain number
#' @slot iter iteration number
#' @slot hasinits is there an initial value
#' @slot inits inits
#' @slot iter_chain_cor
#' @slot t_pred
#' @slot sortindex_cor
#' @slot diagnostic
#' @slot logpdf
#'
#'
PseOptimResult <- setClass(
  Class = 'PseOptimResult',
  slots = c(
    sample_value = 'array',
    dimnames = 'list',
    burnin = 'numeric',
    chains = 'vector',
    iter = 'numeric',
    hasinits = 'logical',
    inits = 'list',
    iter_chain_cor = 'list',
    t_pred = 'list',
    sortindex_cor = 'list',
    diagnostic = 'list',
    logpdf = 'list'
  )
)

#' The TrackOptimResult-class
#'
#' @slot sample_value sampling results for each iteration
#' @slot dimnames names of sampling value
#' @slot burnin burning period
#' @slot chains chain number
#' @slot iter iteration number
#' @slot hasinits is there an initial value
#' @slot inits inits
#'
TrackOptimResult <- setClass(
  Class = 'TrackOptimResult',
  slots = c(
    sample_value = 'array',
    dimnames = 'list',
    burnin = 'numeric',
    chains = 'vector',
    iter = 'numeric',
    hasinits = 'logical',
    inits = 'list',
    diagnostic = 'list',
    logpdf = 'list',
    gene_weight = 'list'
  )
)

#' The Settings-class
#'
#' @slot sample_value sampling results for each iteration
#' @slot dimnames names of sampling value
#' @slot burnin burning period
#' @slot chains chain number
#' @slot iter iteration number
#' @slot hasinits is there an initial value
#' @slot inits inits
#'
Settings <- setClass(
  Class = 'Settings',
  slots = c(
    datasetTitle = "ANY",
    settings="list"
  )
)

#' Create PseOptimResult Object
#'
#' @param sample_value sampling results for each iteration
#' @param dimnames names of sampling value
#' @param burnin burning period
#' @param chains chain number,vector
#' @param iter iteration number
#' @param hasinits is there an initial value
#' @param inits inits
#'
#' @export
#'
CreatePseOptimResultObject <- function(sample_value = NULL,
                                       dimnames = NULL,
                                       burnin = NULL,
                                       chains = NULL,
                                       iter = NULL,
                                       hasinits = FALSE,
                                       inits = list(),
                                       iter_chain_cor = list(),
                                       sortindex_cor = list(),
                                       t_pred = list(),
                                       diagnostic = list(),
                                       logpdf = list()
                                       
){
  
  # check sample_value and dimnames
  if(missing(x = sample_value) && missing(x = dimnames)){
    stop("Must provide either 'sample_value' or 'dimnames'")
  }
  
  # check hasinits and inits
  if(hasinits){
    if(missing(inits)){
      stop("has inits? if has, please provide it in 'inits = '")
    }
    if(!is.list(inits)){
      stop("inits must be list")
    }
  }
  
  # check inits
  object.list <- list(
    'Class' = 'PseOptimResult',
    'sample_value' = sample_value,
    'dimnames' = dimnames,
    'chains' = chains,
    'burnin' = burnin,
    'iter' = iter,
    'hasinits' = hasinits,
    'inits' = inits,
    'iter_chain_cor' = iter_chain_cor,
    't_pred' = t_pred,
    'sortindex_cor' = sortindex_cor,
    'diagnostic' = diagnostic,
    'logpdf' = logpdf
  )
  
  object <- do.call(what = 'new', args = object.list)
  return(object)
}

#' Create TrackOptimResult Object
#'
#' @param sample_value sampling results for each iteration
#' @param dimnames names of sampling value
#' @param burnin burning period
#' @param chains chain number,vector
#' @param iter iteration number
#' @param hasinits is there an initial value
#'
#' @export
#'
CreateTrackOptimResultObject <- function(sample_value = NULL,
                                         dimnames = NULL,
                                         burnin = NULL,
                                         chains = NULL,
                                         iter = NULL,
                                         hasinits = FALSE,
                                         inits = list(),
                                         diagnostic = list(),
                                         logpdf = list(),
                                         gene_weight = list()
){
  
  # check sample_value and dimnames
  if(missing(x = sample_value) && missing(x = dimnames)){
    stop("Must provide either 'sample_value' or 'dimnames'")
  }
  
  # check hasinits and inits
  if(hasinits){
    if(missing(inits)){
      stop("has inits? if has, please provide it in 'inits = '")
    }
    if(!is.list(inits)){
      stop("inits must be list")
    }
  }
  
  # check inits
  object.list <- list(
    'Class' = 'TrackOptimResult',
    'sample_value' = sample_value,
    'dimnames' = dimnames,
    'chains' = chains,
    'burnin' = burnin,
    'iter' = iter,
    'hasinits' = hasinits,
    'inits' = inits,
    'diagnostic' = diagnostic,
    'logpdf' = logpdf,
    'gene_weight' = gene_weight
  )
  
  object <- do.call(what = 'new', args = object.list)
  return(object)
}

#' Create DimReducs Object
#'
#' @export
#'
CreateDimReducObject <- function(pca = NULL,
                                 tsne = NULL,
                                 umap = NULL,
                                 dm = NULL){
  
  # check sample_value and dimnames
  # if(missing(x = reduc_result)){
  #   stop("Must provide 'reduc_result'")
  # }
  
  # pca = new('Class' = 'PCAResult')
  # tsne = new('Class' = 'tSNEResult')
  # umap = new('Class' = 'UMAPResult')
  # dm = new('Class' = 'DMResult')
  
  # pca = NULL
  # tsne = NULL
  # umap = NULL
  # dm = NULL
  
  dimreduc = new(Class = 'DimReducs',
                 PCA = NULL,
                 tSNE = NULL,
                 DM = NULL,
                 UMAP = NULL)
  
  return(dimreduc)
  
}

#' Create DimReducs Object
#'
#' @param datasetTitle
#' @param settings
#' @export
#'
CreateSettings <- function(datasetTitle = NULL,
                           settings = NULL){
  
  Settings = new(Class = 'Settings',
                 datasetTitle = datasetTitle,
                 settings = list())
  return(Settings)
}

#' Create MGPfact Object
#'
#' @param data_matrix gene expression matrix, maybe the expression matrix after normalization or centralization
#' @param count_matrix
#' @param datasetTitle project name
#' @param OptimResult optim result from julia package MGPfact
#' @param MURP MURP result from R package MURP
#' @param DimReducs
#' @param MetaData cell information
#'
#' @export
#'
CreateMGPfactObject <- function(data_matrix = NULL,
                                count_matirx = NULL,
                                datasetTitle = "project",
                                dir = NULL,
                                OptimResult = NULL,
                                MURP = NULL,
                                DimReducs = CreateDimReducObject(),
                                BeautGene = list(),
                                MetaData = NULL,
                                Tree = NULL,
                                GPR = NULL,
                                Command = NULL){
  
  if(is.null(dir)){
    dir = "MGPfact_result"
    dir.create(dir)
  }else{
    if(!dir.exists(dir)) dir.create(dir)
  }
  setwd(dir)
  Initialize()
  
  # check count_matirx and data_matrix
  if(missing(x = data_matrix)){
    stop("Must provide either 'data_matrix'")
  }else if (!missing(x = count_matirx)) {
    if (!inherits(x = count_matirx, what = 'dgCMatrix')) {
      count_matirx <- as(object = as.matrix(x = count_matirx), Class = 'dgCMatrix')
    }
  }
  
  # save to 0_input
  save(data_matrix, file = "0_input/data_matrix.rda")
  save(MetaData, file = "0_input/MetaData.rda")
  
  # check metadata
  if(missing(MetaData)){
    MetaData <- data.frame(row.names = rownames(data_matrix),
                           cellname = rownames(data_matrix),
                           stringsAsFactors = FALSE)
  }else {
    if (is.null(x = rownames(x = MetaData))) {
      stop("Row names not set in metadata. Please ensure that rownames of metadata match column names of data matrix")
    }
    if (length(x = setdiff(x = rownames(x = MetaData), y = rownames(x = data_matrix)))) {
      warning("Some cells in meta.data not present in provided data_matrix matrix.")
      MetaData <- MetaData[intersect(x = rownames(x = MetaData), y = rownames(x = data_matrix)), , drop = FALSE]
    }
    if (is.data.frame(x = MetaData)) {
      save(MetaData, file = "~/x.rda")
      new.MetaData <- data.frame(row.names = rownames(x = data_matrix))
      for (ii in 1:ncol(MetaData)) {
        # cat(ii, "\n")
        # cat(colnames(MetaData)[ii],"\n")
        new.MetaData[rownames(MetaData), colnames(x = MetaData)[ii]] <- MetaData[, ii, drop = FALSE]
      }
      MetaData <- new.MetaData
    }
  }
  
  # check count
  if (!missing(x = count_matirx)){
    # check colnames
    if (anyDuplicated(colnames(x = count_matirx)) | anyDuplicated(rownames(x = data_matrix))) {
      stop( "Non-unique cell names (colnames) present in the input matrix, making unique" )
    }
    # check rownames
    if (anyDuplicated(rownames(x = count_matirx)) | anyDuplicated(colnames(x = data_matrix))) {
      stop( "Non-unique cell names (colnames) present in the input matrix, making unique" )
    }
    if (nrow(x = count_matirx) > 0 && is.null(x = rownames(x = count_matirx))) {
      stop("No feature names (rownames) names present in the input matrix")
    }
    if (any(rownames(x = count_matirx) == '')) {
      stop("Feature names of count_matirx matrix cannot be empty", call. = FALSE)
    }
    if (is.null(x = colnames(x = count_matirx))) {
      stop("No cell names (colnames) names present in the input matrix")
    }
    
    # Ensure row- and column-names are vectors, not arrays
    if (!is.vector(x = rownames(x = count_matirx))) {
      rownames(x = count_matirx) <- as.vector(x = rownames(x = count_matirx))
    }
    if (!is.vector(x = colnames(x = count_matirx))) {
      colnames(x = count_matirx) <- as.vector(x = colnames(x = count_matirx))
    }
  }
  
  # create settings
  settings <- CreateSettings(datasetTitle = datasetTitle,
                             settings = NULL)
  
  # create MGPfact object
  MGPfact <- new(
    Class = 'MGPfact',
    assay = list(count_matrix = count_matirx, data_matrix = data_matrix),
    OptimResult = OptimResult,
    MURP = MURP,
    DimReducs = DimReducs,
    MetaData = MetaData,
    Date = date(),
    Settings = settings,
    BeautGene = list(),
    Tree = list(),
    GPR = list(),
    Command = list())
  if(!is.null(BeautGene)){
    MGPfact@BeautGene = BeautGene
  }
  if(!is.null(Tree)){
    MGPfact@Tree = Tree
  }
  if(!is.null(GPR)){
    MGPfact@GPR = GPR
  }
  return(MGPfact)
}

#' set parameters of settings
#'
#' @description
#' A folder about MGPfact will be created in the current path
#'
Initialize <- function(){
  # create directory in the local
  # dir.create("1_input", showWarnings=FALSE)
  dir.create("0_input", showWarnings=FALSE)
  dir.create("1_murp", showWarnings=FALSE)
  
  dir.create("2_pseudotime", showWarnings=FALSE)
  dir.create("2_pseudotime/2.1_julia_result", showWarnings=FALSE)
  dir.create("2_pseudotime/2.2_binarytree", showWarnings=FALSE)
  dir.create("2_pseudotime/2.3_tbtree", showWarnings=FALSE)
  dir.create("2_pseudotime/2.4_cov_matrix", showWarnings=FALSE)
  
  dir.create("3_tracking", showWarnings=FALSE)
  dir.create("3_tracking/3.1_optim_result", showWarnings=FALSE)
  dir.create("3_tracking/3.2_trajectory", showWarnings=FALSE)
  dir.create("3_tracking/3.3_gene_weight", showWarnings=FALSE)
  
  dir.create("4_differential_genes", showWarnings=FALSE)
  dir.create("5_combine_plot", showWarnings=FALSE)
}

#' set parameters of settings
#'
#' @description
#' Set parameters necessary in optimization
#'
#' @param object MGPfact object
#' @param murp_pc_number the number of principal components of the expression matrix after MURP downsampling
#' @param trajectory_number the number of trajectories to be deconstructed
#' @param pse_optim_iterations the number of optimization iterations to predict pseudotime
#' @param start_murp root point, If not required, set it to a number larger than the number of murps
#' @param chains_number the number of Markov chains
#'
#' @export
#'
SetSettings <- function(object,
                        murp_pc_number = 3,
                        trajectory_number = 3,
                        pse_optim_iterations = 10,
                        start_murp = 1,
                        chains_number = 10){
  
  setting = list(cell_number = nrow(object@assay$data_matrix),
                 gene_number = ncol(object@assay$data_matrix),
                 murp_number = object@MURP$Recommended_K,
                 murp_pc_number = murp_pc_number,
                 trajectory_number = trajectory_number,
                 pse_optim_iterations = pse_optim_iterations,
                 start_murp = start_murp,
                 chains_number = chains_number)
  new_c = setdiff(names(object@Settings@settings),names(setting))
  x = c(object@Settings@settings[new_c],setting)
  object@Settings@settings = x
  return(object)
}

#' assignSettings
#'
#' @description
#' add a single parameter
#'
#' @export
#'
assignSettings <- function(object,
                           slotName = NULL,
                           values = NULL){
  if(slotName %in% names(object@Settings@settings)){
    # object@Settings@settings[[slotName]] = c(object@Settings@settings[[slotName]],values) %>% unique
    object@Settings@settings[[slotName]] = values
  }else{
    len = length(object@Settings@settings)
    object@Settings@settings[[len+1]] = values
    names(object@Settings@settings)[len+1] = slotName
  }
  return(object)
}

#' assignSettingsMulti
#'
#' @description
#' add multiple parameters
#'
#' @param object MGPfact object
#' @param slot a list containing parameters, must be named with the parameter name
#'
#' @export
#'
assignSettingsMulti <- function(object,
                                slot){
  for(i in 1:length(slot) ){
    slotName = names(slot)[i]
    values = slot[[i]]
    object = assignSettings(object, slotName, values)
  }
  return(object)
}

#' writeSettings
#'
#' @description
#' Write out the parameter settings for easy viewing
#'
#' @param object MGPfact object
#' @export
#'
writeSettings <- function(object){
  x = object@Settings@settings
  write.table("\n", file = paste0(getwd(),"/settings.txt"),
              quote = FALSE, append = FALSE, row.names = FALSE, col.names = FALSE, sep = "\n")
  for(i in 1:length(x)){
    write.table(paste0("*", names(x)[i]), file = paste0(getwd(),"/settings.txt"),
                quote = FALSE, append = TRUE, row.names = FALSE, col.names = FALSE, sep = "\n")
    write.table(paste0(x[[i]], collapse = ", "), file = paste0(getwd(),"/settings.txt"),
                quote = FALSE, append = TRUE, row.names = FALSE, col.names = FALSE, sep = "\n\n")
    write.table("\n", file = paste0(getwd(),"/settings.txt"),
                quote = FALSE, append = TRUE, row.names = FALSE, col.names = FALSE, sep = "\n")
  }
}


#' GetCommand
#'
#' @description
#' @noRd
GetCommand <- function(){
  time.stamp <- Sys.time()
  command.name = sys.calls()[[1]]
  command.name = strsplit(as.character(command.name[[1]]),"\\(")[[1]]
  
  argnames <- argnames <- names(x = formals(fun = sys.function(which = sys.parent(n = 1))))
  params <- list()
  p.env <- parent.frame(n = 1)
  argnames <- intersect(x = argnames, y = ls(name = p.env))
  argnames <- setdiff(argnames, c("mamba","murp","metadata","object"))
  for (arg in argnames) {
    param_value <- get(x = arg, envir = p.env)
    # if (inherits(x = param_value)) {
    #   next
    # }
    params[[arg]] <- param_value
  }
  
  # p = list(name = command.name, time.stamp = time.stamp, argnames = argnames, params = params)
  p = list(name = command.name, time.stamp = time.stamp, params = params)
  return(p)
}

#' @export
setMethod("show", "MGPfact",
          function(object) {
            cat("An object of class ", class(object), "\n", sep = "")
            cat("Default data: data_matrix \n")
            cat(" ", ncol(object@assay$data_matrix), " genes across ",
                nrow(object@assay$data_matrix), " samples.\n", sep = "")
            invisible(NULL)
          }
)

#' @export
setMethod("dim", "MGPfact",
          function(x) {
            d1 = dim(x@assay$data_matrix)
            d2 = dim(x@MURP$Recommended_K_cl$centers)
            cat("expression matrix: ", d1[1], " cells, ", d1[2], " genes ",  "\n", sep = "")
            cat("murp matrix: ", d2[1], " murps, ", d2[2], " genes ",  "\n", sep = "")
            return(d1)
          }
)

#' @export
setMethod("$", "MGPfact",
          function(x, name) {
            x = x@MetaData[,name]
            return(x)
            # invisible(NULL)
          }
)

#' @export
setReplaceMethod("$", signature = "MGPfact",
                 function(x, name, value) {
                   x@MetaData <- cbind(x@MetaData, value)
                   colnames(x@MetaData)[ncol(x@MetaData)] = name
                   return(x)
                 })

##### Get Params
#' @name getParams
#' @rdname Settings-class
#' @export getParams
setGeneric(name="getParams",
           def=function(object, ...) standardGeneric("getParams"))

#' @export getParams
setMethod("getParams",
          signature="MGPfact",
          definition = function(object, slotName=NULL)
          {
            if(is.null(slotName))
            {
              print(object@Settings@settings)
            }else{
              x = slotName
              return(object@Settings@settings[[x]])
              # if(slotName=="cell_number") return(object@Settings@settings$cell_number)
              # if(slotName=="gene_number") return(object@Settings@settings$gene_number)
              # if(slotName=="murp_number") return(object@Settings@settings$murp_number)
              # if(slotName=="murp_pc_number") return(object@Settings@settings$murp_pc_number)
              # if(slotName=="trajectory_number") return(object@Settings@settings$trajectory_number)
              # if(slotName=="pse_optim_iterations") return(object@Settings@settings$pse_optim_iterations)
              # if(slotName=="root") return(object@Settings@settings$root)
              # if(slotName=="track_optim_iterations") return(object@Settings@settings$track_optim_iterations)
              # if(slotName=="label") return(object@Settings@settings$label)
            }
            invisible(NULL)
          })

##### Get Pse Object
#' @name getPse
#' @rdname Settings-class
#' @export getPse
setGeneric(name="getPse",
           def=function(object, ...) standardGeneric("getPse"))

#' @export getPse
setMethod("getPse",
          signature="MGPfact",
          definition = function(object)
          {
            return(object@OptimResult$pse)
            invisible(NULL)
          })

##### Get Track Object
#' @name getTrack
#' @rdname Settings-class
#' @export getTrack
setGeneric(name="getTrack",
           def=function(object, ...) standardGeneric("getTrack"))

#' @export getTrack
setMethod("getTrack",
          signature="MGPfact",
          definition = function(object)
          {
            return(object@OptimResult$track)
            invisible(NULL)
          })

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                                    MURP
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' om
#'
#' @Description:
#' omega compute
#'
#' @param x: expression matrix, row is cell, col is gene
#' @export
#'
om <- function(x){
  r = -0.036 + 78.34 * (1/nrow(x)) + 1.85 * 1e-5 * ncol(x)
  return(r)
}

##### MURPDownsampling
#' @param name description
MURPDownsampling <- function(object,
                             omega = 0.5,
                             iter = 10,
                             seed = 723,
                             fast = T,
                             max_murp = 500,
                             cores = 1,
                             pca.center = FALSE,
                             pca.scale = FALSE,
                             plot = F){
  
  object@MURP <- MURP(Data = object@assay$data_matrix,
                      cores = cores,
                      omega = omega,
                      iter = iter,
                      seed = seed,
                      fast = fast,
                      max_murp = max_murp)
  object@MetaData$murp_cluster <- object@MURP$Recommended_K_cl$cluster
  # object@MURP$centersPCA <- prcomp(object@MURP$Recommended_K_cl$centers, center = pca.center, scale. = pca.scale)
  object <- MURPPCA(object, pca.center = pca.center, pca.scale = pca.scale)
  
  if(plot){
    ggsave("1_murp/murp_bic_k.pdf", MURPNestedGridPlot(object@MURP), width = 3.5, height = 3.5)
    
    pl = PCANestedGridPlot(pca_result = object@MURP$centersPCA, sd_cutoff = 1, max_pc = 100)
    p = wrap_plots(pl, ncol = 3, guides = "collect")
    ggsave("1_murp/pca_scree_30pc.pdf",p, width = 10, height = 3.7)
  }
  
  command = GetCommand()
  object@Command$murp$MURPDownsampling = command
  return(object)
}

##### MURPDownsampling
#' @param name description
MURPPCA <- function(object,
                    pca.center = FALSE,
                    pca.scale = FALSE){
  
  object@MURP$centersPCA <- prcomp(object@MURP$Recommended_K_cl$centers, center = pca.center, scale. = pca.scale)
  
  command = GetCommand()
  object@Command$murp$MURPPCA = command
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                                   convert
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' SaveMURPDatToJulia
#'
#' @Description:
#' extract murp data to julia environment
#' @param object MGPfact object
#' @param murp_pc_number the number of principal components of the expression matrix after MURP downsampling
#' @export
#'
SaveMURPDatToJulia <- function(object,
                               murp_pc_number = 5){

  # murp_matrix = ct@MURP$Recommended_K_cl$centers
  # save(murp_matrix, file = "1_murp/murp_matrix.rda")
  Q = murp_pc_number
  pca = object@MURP$centersPCA$x
  murp_matrix_pca = pca[,1:Q,drop=FALSE]
  save(murp_matrix_pca, file = "1_murp/murp_matrix_pca.rda")
}

#' RunningmodMGPpseudoT
#'
#' @Description:
#' running modMGPpseudoT model in julia
#'
#' @param object MGPfact object
#' @param julia_home Julia's bin path
#' @export
#'
RunningmodMGPpseudoT <- function(object,
                                 julia_home,
                                 seed = 723,
                                 cores = 1){
  require(JuliaCall)
  julia_setup(JULIA_HOME = julia_home)

  cat("/// parameters setting \n")
  # list(getParams(object, "trajectory_number"),
  #      getParams(object, "murp_pc_number"),
  #      getParams(object, "pse_optim_iterations"),
  #      getParams(object, "start_murp"),
  #      getParams(object, "chains_number")) -> x
  # save(x, file = "~/x.rda")

  cmd = list(paste0("trajectory_number=", getParams(object, "trajectory_number")),
             paste0("pc_number=", getParams(object, "murp_pc_number")),
             paste0("iterations=", getParams(object, "pse_optim_iterations")),
             paste0("start_id=[", paste0(getParams(object, "start_murp"),collapse = ";"),"]"),
             paste0("chains_number=", getParams(object, "chains_number")))
  invisible(sapply(cmd, julia_command))
  # cmd = paste0("trajectory_number=", getParams(ct, "trajectory_number"),
  #              "; pc_number=", getParams(ct, "murp_pc_number"),
  #              "; iterations=", getParams(ct, "pse_optim_iterations"),
  #              "; start_id=", getParams(ct, "start_murp"),
  #              "; chains_number=", getParams(ct, "chains_number"))
  # cmd = paste0("trajectory_number=", getParams(ct, "trajectory_number"),
  #              "\n pc_number=", getParams(ct, "murp_pc_number"),
  #              "\n iterations=", getParams(ct, "pse_optim_iterations"),
  #              "\n start_id=", getParams(ct, "start_murp"),
  #              "\n chains_number=", getParams(ct, "chains_number"))
  # julia_command(paste0("trajectory_number=", getParams(ct, "trajectory_number"),
  #                      "; pc_number=", getParams(ct, "murp_pc_number"),
  #                      "; iterations=", getParams(ct, "pse_optim_iterations"),
  #                      "; start_id=", getParams(ct, "start_murp"),
  #                      "; chains_number=", getParams(ct, "chains_number")))

  julia_command("using Distributed")
  julia_command(paste0("addprocs(",cores,")"))

  cat("/// load data \n")

  cmd = '@everywhere using MGPfact, Mamba, RData, Distributions
data_path="1_murp/murp_matrix_pca.rda"
yx = RData.load(data_path)
yx = yx["murp_matrix_pca"][:,1:pc_number]'
#   cmd = '@everywhere using MGPfact, Mamba, RData, Distributions
# data_path="1_murp/murp_matrix_pca.rda"
# yx = RData.load(data_path)
# yx = yx["murp_matrix_pca"][:,1:pc_number]'
  cmd = gsub("\n", ";", cmd)
  julia_command(cmd)

  cat("/// running model \n")

  julia_command("using Random")
  julia_command(paste0("Random.seed!(",seed,")"))

  cat("/// running model \n")

  cmd = 'model = MGPfact.modMGPpseudoT()
SC,inits,scheme = MGPfact.Initialize(yx, pc_number, trajectory_number, start_id, iterations, chains_number)
setinputs!(model, SC)
setinits!(model, inits)
setsamplers!(model, scheme)
@time sim = mcmc(model, SC, inits, iterations, burnin = 0, chains = chains_number)'
  cmd = gsub("\n", ";", cmd)
  julia_command(cmd)

  cat("/// saving \n")

  cmd = 'using JLD2
write(string("2_pseudotime/2.1_julia_result/iter",sim.model.iter,"_bi",sim.model.burnin,".jls"), sim)
@save string("2_pseudotime/2.1_julia_result/inits.jld2") inits'
  cmd = gsub("\n", ";", cmd)
  julia_command(cmd)

  object = ImportPseResult(object = object, init = TRUE)
  return(object)
}

#' ReadPseSim
#'
#' @Description
#' read and convert the object optimized by julia into R
#'
#' @param sim julia result
#' @param init
#' Logical value, whether to import the initial value.
#' it can be viewed in the directory: 2_pseudotime/2.1_julia_result/inits.jld2
#'
#' @export
#'
ReadPseSim <- function(sim = NULL,
                       julia_home = NULL,
                       init = FALSE){
  library(JuliaCall)
  julia_setup(JULIA_HOME = julia_home)

  ## load sim.jls
  cmd = 'using Distributed, RData
nc = 10
addprocs(nc)
@everywhere using Pkg
@everywhere using MGPfact, Mamba, RData, Distributions, KernelFunctions'
  cmd=paste0(cmd,'\nsim1=read("2_pseudotime/2.1_julia_result/',sim,'",ModelChains)')
  cmd = gsub("\n", ";", cmd)
  julia_command(cmd)

  ## load inits
  if(init){
    julia_command("using JLD2")
    julia_command('@load "2_pseudotime/2.1_julia_result/inits.jld2"')
  }
}

#' ImportPseSimToR
#'
#' @Description
#' Import the pseudotime object optimized by julia into R
#'
#' @param sim julia result
#' @param init
#' Logical value, whether to import the initial value.
#' it can be viewed in the directory: 2_pseudotime/2.1_julia_result/inits.jld2
#'
#' @export
#'
ImportPseResult <- function(object,
                            init = FALSE,
                            julia_home = FALSE){

  library(JuliaCall)
  julia_setup(JULIA_HOME = julia_home)

  if(init){
    cmd = 'chains = sim.chains
names2 = sim.names
values1 = sim.value
burnin = sim.model.burnin
iter =  sim.model.iter
hasinits = sim.model.hasinits
@rput chains burnin iter inits names2 values1'
    cmd = gsub("\n", ";", cmd)
    julia_command(cmd)
    OptimResult <- CreatePseOptimResultObject(sample_value = values1,
                                              dimnames = list(paste0("iter",1:dim(values1)[1]),
                                                              names=names2,
                                                              paste0("chain",1:dim(values1)[3])),
                                              burnin = burnin,
                                              chains = chains,
                                              iter = iter,
                                              inits = as.list(inits),
                                              hasinits = TRUE)
  }else{
    cmd = 'chains = sim.chains
names2 = sim.names
values1 = sim.value
burnin = sim.model.burnin
iter =  sim.model.iter
hasinits = sim.model.hasinits
@rput chains burnin iter names2 values1'
    cmd = gsub("\n", ";", cmd)
    julia_command(cmd)
    OptimResult <- CreatePseOptimResultObject(sample_value = values1,
                                              dimnames = list(paste0("iter",1:dim(values1)[1]),
                                                              names=names2,
                                                              paste0("chain",1:dim(values1)[3])),
                                              burnin = burnin,
                                              chains = chains,
                                              iter = iter,
                                              # inits = as.list(inits),
                                              hasinits = FALSE)
  }

  object@OptimResult$pse = OptimResult
  return(object)
}

#' GetSigmaFromJulia
#'
#' @Description:
#' get the covariance matrix
#'
#' @param object MGPfact object
#' @param pse_sdf tThe data frame obtained by the function GetPseSdf
#' @param load logical value, whether to load julia's environment
#'
#' @export
#'
GetSigmaFromJulia <- function(object,
                              load = FALSE,
                              julia_home = NULL,
                              init = T){
  ## 1. load
  if(load){
    ReadPseSim(sim = paste0("iter",getParams(object, "pse_optim_iterations"),"_bi0.jls"),
               julia_home = julia_home,
               init = init)
#    cmd = 'using Distributed, RData
# nc = 10
# addprocs(nc)
# @everywhere using Pkg
# @everywhere using MGPfact, Mamba, RData, Distributions, KernelFunctions'
# cmd = gsub("\n", ";", cmd)
# julia_command(cmd)
  }

  ## 2. prepare in julia
  sdf = GetMURPInfo(object=object)
  # cat(dim(sdf))
  julia_assign("sdf",sdf)
  # julia_command("using RCall; @rget sdf")

  julia_command("using KernelFunctions")
  cmd = 'P = size(sdf)[1]
L = size(filter(x->occursin("Tb", string(x)), names(sdf)))[1]
Q = size(filter(x->occursin("PC", string(x)), names(sdf)))[1]
Tb = Vector(sdf[1,filter(x->occursin("Tb", string(x)), names(sdf))])
lambda1 = Vector(sdf[1,filter(x->occursin("lambda1", string(x)), names(sdf))])
lambda2 = Vector(sdf[1,filter(x->occursin("lambda2", string(x)), names(sdf))])
lambda3 = Vector(sdf[1,filter(x->occursin("lambda3", string(x)), names(sdf))])
lambda = hcat(lambda1,lambda2,lambda3)
T = sdf.T
X = sdf.X
m_t = sdf.mt[1]
s2_t = sdf.s2t[1]
s2_x = sdf.s2x[1]
rho = sdf.rho[1]'
  cmd = gsub("\n", ";", cmd)
  julia_command(cmd)

  julia_command('tmp = sdf[:,filter(x->occursin(r"^C_", string(x)), names(sdf))]')
  julia_command("C = convert(Array{Int},Matrix(tmp)')")

  ## 3. get cov matrix
  cmd = 'sigma = 1.0
delta = 1.0
P = Int(P)
L = Int(L)
S_list = [( tb = Tb[i];
            cl = C[i,:];
            le = lambda[i,:];
            k = MGPfact.track.bifurcateKernel(tb, le, sigma, 0.0, delta);
            x = [MGPfact.track.Trejactory(T[j], cl[j], X[j]) for j in 1:P];
            ss = KernelFunctions.kernelmatrix(k, x);
            Matrix(Hermitian(ss .+ 1E-6 .* Matrix(I, P, P)))) for i in 1:L ]
S = sum(S_list)'
  cmd = gsub("\n", ";", cmd)
  cmd = gsub(";;", ";", cmd)
  julia_command(cmd)

  ## 4. get 100 points sigma
  # cat("2 \n")
  cmd = 'P_ = 100
T1_ = collect(1:(P_) )/(P_)
T2_ = collect(1:(P_) )/(P_)
T_ = vec(vcat(T1_,T2_))
C1_ = fill(1, L, P_)
C2_ = fill(2, L, P_)
C_ = hcat(C1_,C2_)
sort_index = sortperm(T_)
T_ = T_[sort_index]
C_ = C_[:,sort_index]
P_ = P_*2
L = Int(L)
X_ = rand(truncated(Normal(0, sqrt(s2_x[1])), 0.0, 1.0), P_)
S_list_ = [( tb = Tb[i];
             le = lambda[i,:];
             k = MGPfact.track.bifurcateKernel(tb, le, sigma, 0.0, delta);
             x = [MGPfact.track.Trejactory(vcat(T,T_)[j], hcat(C,C_)[i,j], vcat(X, X_)[j]) for j in 1:(P+P_)];
             ss = KernelFunctions.kernelmatrix(k, x);
             Matrix(Hermitian(ss .+ 1E-6 .* Matrix(I, P+P_, P+P_)))) for i in 1:L ]
S_ = sum(S_list_)'
  cmd = gsub("\n", ";", cmd)
  cmd = gsub(";;", ";", cmd)
  julia_command(cmd)

  ## 5. save in julia
  cat("3 \n")
  cmd = 'using JLD2
  @save "2_pseudotime/2.4_cov_matrix/s_.jld2" S_ S_list_ C_ T_ X_ P_ rho m_t s2_t s2_x Tb lambda L Q P
  @save "2_pseudotime/2.4_cov_matrix/s.jld2" S S_list C T X rho m_t s2_t s2_x Tb lambda L Q P
  @rput S_ S_list_ C_ T_ X_ P_ rho m_t s2_t s2_x Tb lambda L Q P
  @rput S S_list C T X rho m_t s2_t s2_x Tb lambda L Q P'
  cmd = gsub("\n", ";", cmd)
  julia_command(cmd)

  ## 6. save in R
  S_list_ = as.list(S_list_)
  S_list = as.list(S_list)
  save(S_,S_list_,C_,T_,X_,P_,rho,m_t,s2_t,s2_x,Tb,lambda,L,Q,P, file = "2_pseudotime/2.4_cov_matrix/s_.rda")
  save(S,S_list,C,T,X,rho,m_t,s2_t,s2_x,Tb,lambda,L,Q,P, file = "2_pseudotime/2.4_cov_matrix/s.rda")

  object@GPR$reg_cov = list(cov = S_, cov_l = S_list_,
                            C = C_, T = T_, X = X_, P = P_,
                            rho = rho, m_t = m_t, s2_t = s2_t, s2_x = s2_x, Tb = Tb, lambda = lambda)
  object@GPR$murp_cov = list(cov = S, cov_l = S_list,
                             C = C, T = T, X = X,
                             rho = rho, m_t = m_t, s2_t = s2_t, s2_x = s2_x, Tb = Tb, lambda = lambda)
  return(object)
}

#' ImportLogpdf
#'
#' @Description:
#' load probability density values for different parameters
#'
ComputeLogpdf <- function(object,
                          import_sim = FALSE,
                          julia_home=NULL ){

  if(import_sim){
    ReadPseSim(sim =  paste0("iter",getParams(object, "pse_optim_iterations"),"_bi0.jls"),
               julia_home = julia_home, init = TRUE)
  }

  cmd = 'using JLD2
chain = size(sim1)[3]
iterations = size(sim1.value)[1]
burnin = size(sim1)[1] - size(sim1.value)[1]
sigma = sim1.model.nodes[:sigma]
delta = sim1.model.nodes[:delta]
P = sim1.model.nodes[:P]
L = sim1.model.nodes[:L]
Q = sim1.model.nodes[:Q]'
  cmd = gsub("\n", ";", cmd)
  julia_command(cmd)

  # cat("all \n")
  cmd='value = logpdf(sim1)
l_all = zeros(iterations, chain)
for i in 1:iterations l_all[i, :] .= value.value[i,1,:] end'
  cmd = gsub("\n", ";", cmd)
  julia_command(cmd)

  # cat("rho & pi \n")
  cmd = 'rho = hcat([sim1[burnin:(burnin+iterations), "rho", :].value ]...)
a = hcat([[sim1[burnin:(burnin+iterations), "a[$i, $j]", :].value for i in 1:L] for j in 1:P]...)
p = hcat([[sim1[burnin:(burnin+iterations), "p[$i, $j]", :].value for i in 1:L] for j in 1:P]...)
C = hcat([[sim1[burnin:(burnin+iterations), "C[$i, $j]", :].value for i in 1:L] for j in 1:P]...)
l_rho = zeros(iterations, chain)
for i in 1:iterations l_rho[i, :] .= logpdf.(truncated(Beta(1.0, 1.0), 1E-8, 0.99999999), vcat(rho[i,:,1:chain]...)) end
l_p = zeros(iterations, chain)
for i in 1:iterations l_p[i, :] .= [sum(hcat([vcat([logpdf(Beta(a[l,j][i,1,ch]*(1-rho[i, 1, ch]), a[l,j][i,1,ch]*rho[i, 1, ch]), p[l,j][i,1,ch]) for l in 1:L]...) for j in 1:P]...)) for ch in 1:chain] end'
  cmd = gsub("\n", ";", cmd)
  julia_command(cmd)

  # cat("Tb \n")
  cmd = 'Tb = hcat([sim1[burnin:(burnin+iterations), "Tb[$i]", :].value for i in 1:L ]...)
l_tb = zeros(iterations, chain)
for i in 1:iterations l_tb[i, :] .= [sum([logpdf(truncated(Gamma(l, 1/L), 0.0, 1.0),Tb[i, l, ch]) for l in 1:L]) for ch in 1:chain] end'
  cmd = gsub("\n", ";", cmd)
  julia_command(cmd)

  # cat("lambda \n")
  cmd = 'lambda = hcat([[sim1[burnin:(burnin+iterations), "lambda[$i, $j]", :].value for i in 1:L] for j in 1:3]...)
lambda1 = hcat(lambda[:, 1]...)
lambda2 = hcat(lambda[:, 2]...)
lambda3 = hcat(lambda[:, 3]...)
l_l1 = zeros(iterations, chain)
for i in 1:iterations l_l1[i, :] .= [sum(logpdf(sim1.model.nodes[:lambda1], lambda1[i,:,ch])) for ch in 1:chain] end
l_l2 = zeros(iterations, chain)
for i in 1:iterations l_l2[i, :] .= [sum(logpdf(sim1.model.nodes[:lambda2], lambda2[i,:,ch])) for ch in 1:chain] end
l_l3 = zeros(iterations, chain)
for i in 1:iterations l_l3[i, :] .= [sum(logpdf(sim1.model.nodes[:lambda3], lambda3[i,:,ch])) for ch in 1:chain] end'
  cmd = gsub("\n", ";", cmd)
  julia_command(cmd)

  # cat("T \n")
  cmd = 'm_t = hcat([sim1[burnin:(burnin+iterations), "m_t", :].value ]...)
l_mt = zeros(iterations, chain)
for i in 1:iterations l_mt[i, :] .= logpdf.(sim1.model.nodes[:m_t], vcat(m_t[i,:,1:chain]...)) end
s2_t = hcat([sim1[burnin:(burnin+iterations), "s2_t", :].value ]...)
l_s2t = zeros(iterations, chain)
for i in 1:iterations l_s2t[i, :] .= logpdf.(sim1.model.nodes[:s2_t], vcat(s2_t[i,:,1:chain]...)) end'
  cmd = gsub("\n", ";", cmd)
  julia_command(cmd)

  cmd = 'T = hcat([sim1[burnin:(burnin+iterations), "T[$i]", :].value for i in 1:P]...)
l_t = zeros(iterations, chain)
rt = sim1.model.nodes[:rt]'
  cmd = gsub("\n", ";", cmd)
  julia_command(cmd)

  cmd = 'for i in 1:iterations
    tmp = zeros(P,chain)
    for ch in 1:chain
      for i2 in 1:P
        if i2 in rt
          tmp[i2,ch]=logpdf(truncated(Normal(0, 0.01), 0.0, 1.0), T[i,i2,ch])
        else
          m_t2 = m_t[i,1,ch]
          s2_t2 = s2_t[i,1,ch]
              tmp[i2,ch]=logpdf(truncated(Normal(m_t2, sqrt(s2_t2)), 0.0, 1.0), T[i,i2,ch])
          end
      end
  end
  l_t[i,:] .=  vcat(sum(tmp, dims=1)...)
end'
  julia_command(cmd)

  cat("X \n")
  cmd = 's2_x = hcat([sim1[burnin:(burnin+iterations), "s2_x", :].value ]...)
l_s2x = zeros(iterations, chain)
for i in 1:iterations l_s2x[i, :] .= logpdf.(sim1.model.nodes[:s2_x], vcat(s2_x[i,:,1:chain]...)) end
X = hcat([sim1[burnin:(burnin+iterations), "X[$i]", :].value for i in 1:P]...)
l_x = zeros(iterations, chain)
for i in 1:iterations l_x[i, :] .= [sum(logpdf(truncated(Normal(0, sqrt(s2_x[i,1,ch])), -1.0, 1.0), X[i, :, ch])) for ch in 1:chain] end'
  cmd = gsub("\n", ";", cmd)
  julia_command(cmd)

  # cat("9 \n")
  # cmd='S = [[
  #        Matrix(Hermitian(
  #         sum([(
  #             tb = Tb[i, l, ch];
  #             cl = hcat(C[l,:]...)[i,:,ch];
  #             le = [lambda1[i,l,ch], lambda2[i,l,ch], lambda3[i,l,ch]];
  #             k = MGPfact.pseudot.bifurcateKernel(tb, le, sigma, 0.0, delta);
  #             gp = [MGPfact.pseudot.Trejactory(T[i, p, ch], cl[p], X[i, p, ch]) for p in 1:P];
  #             KernelFunctions.kernelmatrix(k, gp)
  #             )
  #         for l in 1:L]) .+ 1E-6 .* Matrix(I, P, P) ))
  # for i in 1:5] for ch in 1:10]'
  #     julia_command(cmd)

  cmd='@rput l_all l_rho l_tb l_mt l_s2t l_s2x l_l1 l_l2 l_l3 l_t l_x l_p chain iterations
@save "2_pseudotime/2.1_julia_result/logpdf.jld2" l_all l_rho l_tb l_mt l_s2t l_s2x l_l1 l_l2 l_l3 l_t l_x l_p chain iterations'
  cmd = gsub("\n", ";", cmd)
  julia_command(cmd)

  logpdf = list(l_all = l_all,
                l_p = l_p,
                l_rho = l_rho,
                l_tb = l_tb,
                l_l1 = l_l1, l_l2 = l_l2, l_l3 = l_l3,
                l_x = l_x, l_s2x = l_s2x,
                l_t = l_t, l_mt = l_mt, l_s2t = l_s2t)

  mamba=object@OptimResult$pse
  mamba@logpdf$logpdf = logpdf

  chain = length(mamba@chains)
  iterations = mamba@iter
  logpdf_cor = lapply(c("spearman","pearson"), function(m){

    lapply(names(logpdf), function(x){
      lp = logpdf[[x]]
      apply(lp, 2, function(y) cor(y, 1:iterations, method = m) )
    }) %>% do.call(rbind,.) -> cor_df
    dimnames(cor_df) = list(names(logpdf),paste0("chain",1:chain))
    return(cor_df)
  })
  names(logpdf_cor) = c("spearman","pearson")

  corr_spearman = logpdf_cor[[1]]
  corr_pearson = logpdf_cor[[1]]

  write.csv(corr_spearman, file = "2_pseudotime/2.1_julia_result/logpdf_iter_cor_spearman.csv")
  write.csv(corr_pearson, file = "2_pseudotime/2.1_julia_result/logpdf_iter_cor_pearson.csv")

  save(corr_spearman, corr_pearson, file = "2_pseudotime/2.1_julia_result/logpdf_iter_cor.rda")
  mamba@logpdf$cor = list(pearson = corr_pearson,
                          spearman = corr_spearman)
  object@OptimResult$pse=mamba
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                                    cor
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' GetIterCor
#'
#' @description
#' use this function, obtain correlations under different
#' iterations of different chains
#'
#'
#' @param mamba mamba result object
#' @param murp murp object
#' @param metadata attribute of cell
#' @param iteration_list a list, list(c(1, 2000),c(1,4000))
#' @param adjust TRUE
#'
GetIterCor <- function(object,
                       iteration_list = NULL,
                       chains = NULL){
  
  require(dplyr)
  mamba = getPse(object)
  murp = object@MURP
  metadata = object@MetaData
  
  if(is.null(chains)){
    chains <- mamba@chains
  }
  if(is.null(iteration_list)){
    seq_tmp <- c(seq(0,mamba@iter,2000))[-1]
    iteration_list <- lapply(seq_tmp, function(x) c(1, x))
  }
  
  svalue <- mamba@sample_value
  parnames <-  mamba@dimnames[[2]]
  iter_chain_cor <- mamba@iter_chain_cor
  
  cellnum <- nrow(metadata)
  cluster_label <- murp$Recommended_K_cl$cluster
  
  # get T_pred
  T_pred <- lapply(iteration_list, function(iter){
    # cat(iter)
    c = lapply(chains, function(ch){
      if(iter[1]==iter[2]){
        svalue[iter[1]:iter[2], grep(pattern="T\\[", parnames),ch]
      }else{
        apply(svalue[iter[1]:iter[2], grep(pattern="T\\[", parnames),ch], 2, mean)
      }
      
    })
    names(c) = paste0("chain",chains)
    c
  })
  names(T_pred) <- paste0("iter_", lapply(iteration_list, paste0, collapse="_"))
  
  # get pseudotime
  pseudotime <- lapply(1:length(iteration_list), function(i){
    tmp = T_pred[[i]]
    pt = lapply(1:length(chains), function(ch){
      # GetPseudoTt(CellNum = cellnum, pred_t = tmp[[ch]], Cluster_Label = cluster_label)
      MapMURPLabelToAll(vecc =  tmp[[ch]], orig = cluster_label)
    })
    names(pt) = paste0("chain",chains)
    pt
  })
  names(pseudotime) <- paste0("iter_",iteration_list)
  
  # get cor
  if("truetime" %in% colnames(metadata)){
    truetime <- as.numeric(as.character(metadata$truetime))
    ct_spearman <- matrix(0, nc = length(iteration_list), nr = length(chains))
    colnames(ct_spearman) <- sapply(iteration_list, function(c){paste0("iter_",c[1],"_",c[2])})
    rownames(ct_spearman) <- paste0('chain', chains)
    ct_pearson <- ct_spearman
    for (k in 1:length(iteration_list)){
      for (i in 1:length(chains) ){
        ct_spearman[i, k] <- cor(truetime, pseudotime[[k]][[i]], method = 'spearman')
        ct_pearson[i, k] <- cor(truetime, pseudotime[[k]][[i]], method = 'pearson')
      }
    }
    
  }else{
    ct_spearman = ct_pearson = NULL
  }
  
  iter_chain_cor = list(t_pred = T_pred,
                        pseudot = pseudotime,
                        spearman = ct_spearman,
                        pearson = ct_pearson)
  
  # 不管矫正不矫正都加,get adjust t_pred
  chain_t = do.call(rbind, tail(T_pred,1)[[1]])
  chain_cor = cor(t(chain_t))
  c = which.max(apply(chain_cor,1,mean))
  
  # 根据所有的chain的均值得到要反转的chain是哪几个
  if(length(which(chain_cor[c,] > 0)) < (nrow(chain_t)/2) ){
    ch_rev <- which(chain_cor[c,] > 0)
  } else{
    ch_rev <- which(chain_cor[c,] < 0)
  }
  # cat(ch_rev,"\n")
  
  # reverse
  adjust_T_pred <- lapply(1:length(T_pred), function(i){
    chain_t = do.call(rbind, T_pred[[i]])
    chain_t[ch_rev,] <- 1 - chain_t[ch_rev,]
    chain_t
  })
  names(adjust_T_pred) = paste0("iter_", lapply(iteration_list, paste0, collapse="_"))
  
  # get adjust pseudotime
  adjust_pseudotime <- lapply(1:length(adjust_T_pred), function(i){
    tmp = adjust_T_pred[[i]]
    pt = lapply(1:length(chains), function(ch){
      # GetPseudoTt(CellNum = cellnum, pred_t = tmp[ch,], Cluster_Label = cluster_label)
      MapMURPLabelToAll(vecc =  tmp[ch,], orig = cluster_label)
    })
    names(pt) = paste0("chain",chains)
    pt
  })
  names(adjust_pseudotime) <- paste0("iter_", lapply(iteration_list, paste0, collapse="_"))
  
  # get adjust cor
  adjust_ct_spearman <- matrix(0, nc = length(iteration_list), nr = length(chains))
  colnames(adjust_ct_spearman) <- sapply(iteration_list, function(c){paste0("iter_",c[1],"_",c[2])})
  rownames(adjust_ct_spearman) <- paste0('chain', chains)
  adjust_ct_pearson <- adjust_ct_spearman
  
  if("truetime" %in% colnames(metadata)){
    for (k in 1:length(iteration_list)){
      for (i in 1:length(chains) ){
        adjust_ct_spearman[i, k] <- cor(truetime, adjust_pseudotime[[k]][[i]], method = 'spearman')
        adjust_ct_pearson[i, k] <- cor(truetime, adjust_pseudotime[[k]][[i]], method = 'pearson')
      }
    }
  }else{
    adjust_ct_spearman = adjust_ct_pearson = NULL
  }
  
  # prepare
  if(length(ch_rev)==0){
    adjust_chain = NULL
  }else{
    adjust_chain = paste0("chain",ch_rev)
  }
  
  iter_chain_cor = list(t_pred = T_pred,
                        pseudot = pseudotime,
                        spearman = ct_spearman,
                        pearson = ct_pearson,
                        adjust_chain = adjust_chain,
                        adjust_t_pred = adjust_T_pred,
                        adjust_pseudot = adjust_pseudotime,
                        adjust_spearman = adjust_ct_spearman,
                        adjust_pearson = adjust_ct_pearson )
  
  # save to object
  mamba@iter_chain_cor = iter_chain_cor
  object@OptimResult$pse = mamba
  command = GetCommand()
  object@Command$pse$GetIterCor = command
  
  return(object)
  
}

#' GetPredT
#'
#' @description
#' use this function, obtain correlations under different
#' iterations of different chains
#'
#' @param mamba mamba result object
#' @param murp murp object
#' @param metadata attribute of cell
#' @param iteration_list a list, list(c(1, 2000),c(1,4000))
#' @param chains c(1,2,3,4)
#' @param adjust TRUE
#' @param iter_chain_cor_index index of iter_chain_cor from function GetIterCor
#'
GetPredT <- function(object,
                     chains = c(1:10),
                     filter_chain = TRUE,
                     mean_th = 0.45,
                     adjust = TRUE){
  
  # object=ct
  # chains = c(1:10)
  # filter_chain = TRUE
  # mean_th = 0.45
  # adjust = FALSE
  
  require(dplyr)
  mamba = getPse(object)
  murp = object@MURP
  metadata = object@MetaData
  
  t_pred_list = mamba@t_pred
  iter_chain_cor = mamba@iter_chain_cor
  
  # get data from iter_chain_cor
  if(object@Settings@settings$chains_number>1 & adjust){
    t_pred = iter_chain_cor$adjust_t_pred
    pseudot = iter_chain_cor$adjust_pseudot
    ch_rev = iter_chain_cor$adjust_chain
  }else{
    t_pred = iter_chain_cor$t_pred
    pseudot = iter_chain_cor$pseudot
    ch_rev = NULL
  }
  pseudot_r = lapply(pseudot, function(pse){ do.call(rbind,pse) })
  names(pseudot_r) = names(pseudot)
  
  # filter chains
  if(object@Settings@settings$chains_number > 1 & filter_chain){
    # cat(1)
    chain_cor = cor(t( tail(pseudot_r,1)[[1]][chains,,drop=FALSE] ))
    chain_cor_mean = apply(chain_cor, 2, mean)
    filter_chain = rownames(chain_cor)[which(chain_cor_mean > mean_th)]
  }else{
    filter_chain = paste0("chain", chains)
  }
  
  # get mean pseudot
  mean_pseudot_r = lapply(pseudot_r, function(df) apply(df[filter_chain, ,drop=FALSE], 2, mean) )
  
  # get cor
  if("truetime" %in% colnames(metadata)){
    truetime = as.numeric(as.character(metadata$truetime))
    ct_cor <- matrix(0, nc = length(names(pseudot_r)), nr = 2)
    colnames(ct_cor) <- names(pseudot_r)
    rownames(ct_cor) <- c("spearman", "pearson")
    for (k in 1:length(pseudot_r)){
      ct_cor["spearman", k] <- cor(truetime, mean_pseudot_r[[k]], method = 'spearman')
      ct_cor["pearson", k] <- cor(truetime, mean_pseudot_r[[k]], method = 'pearson')
    }
  }else{
    ct_cor = NULL
  }
  
  t_pred_list = list(t_pred = t_pred,
                     pseudot = pseudot,
                     mean_pseudot = mean_pseudot_r,
                     cor = ct_cor,
                     keep_chain = filter_chain,
                     adjust = adjust,
                     adjust_chain = ch_rev)
  
  mamba@t_pred = t_pred_list
  object@OptimResult$pse = mamba
  
  command = GetCommand()
  object@Command$pse$GetPredT = command
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                                    param
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' GetParamMeanChain
#'
#' @description
#' Get the average value of each chain for any parameter
#'
#' @param mamba mamba object
#' @param murp murp result
#' @param param parameter name
#' @param whichPos 'all' / c(1,2,3)
#' @param iter_range Iteration range, default is c(1, 3000)
#' @param cores threads
#'
#' @export
GetParamMeanChain <- function(mamba = NULL,
                              murp = NULL,
                              param = NULL,
                              whichPos = 'all',
                              iter_range = c(1,3000),
                              cores = 1
){
  
  require(reshape)
  require(stringr)
  require(pbmcapply)
  require(doParallel)
  
  # if(is.null(dimnames(mamba_object$forcast_time@sample_value))){
  #   dimnames(mamba_object$forcast_time@sample_value) <- mamba_object$forcast_time@dimnames
  # }
  chains <- mamba@chains
  K <- length(grep("T",mamba@dimnames[[2]]))
  all_params <- mamba@dimnames[[2]]
  sample_value <- mamba@sample_value
  param_number <- as.data.frame.numeric(table(str_split_fixed(all_params,"\\[",2)[,1]))
  
  dimnames(sample_value) <- list(paste0("iter",1:dim(sample_value)[1]),
                                 all_params,
                                 paste0("chain",1:dim(sample_value)[3]))
  
  # check param number
  if(param %in% rownames(param_number) ){
    
    if(param_number[param,]==1){
      
      # t = sample_value[,which(all_params==param),]
      t = sample_value[,grep(param,all_params),,drop=FALSE]
      # if(is.vector(t)){
      #   t = as.data.frame(t)
      # }
      df = mclapply(chains, function(ch){
        if(is.numeric(iter_range) & length(iter_range)==2){
          if(iter_range[2]>=iter_range[1] ){
            tmp = t[iter_range[1]:iter_range[2],,ch,drop=F]
          }
        }else{
          tmp = t[,ch,drop=F]
        }
        mean(tmp)
      }, mc.cores = cores)
      df = do.call(rbind, df)
      rownames(df) = paste0("chain", chains)
      
    }else{
      if(whichPos == 'all'){
        t = sample_value[,grep(pattern=paste0(param, "\\["),all_params),,drop=FALSE]
        # if(is.vector(t)){
        #   t = as.data.frame(t)
        # }
        dimnames(t)[[2]] = all_params[grep(pattern=paste0(param, "\\["),all_params)]
      }else{
        if(param == "rho" | param == "Y"){
          index_name = unlist(lapply(whichPos, function(i) paste0(param,"[",i[1],", ",i[2],"]")) )
          t = sample_value[, index_name, drop=F]
          dimnames(t)[[2]] = index_name
        }else if(param == "pi"){
          index_name = unlist(lapply(whichPos, function(i) paste0(param,"[",i[1],", ",i[2],", ",i[3],"]")) )
          t = sample_value[, index_name, drop=F]
          dimnames(t)[[2]] = index_name
        }else{
          t = sample_value[, paste0(param, "[",whichPos,"]"), ]
          dimnames(t)[[2]] = paste0(param, "[",whichPos,"]")
        }
      }
      
      # compute mean value
      df <- mclapply(chains, function(ch){
        if(is.numeric(iter_range) & length(iter_range)==2){
          if(iter_range[2]>=iter_range[1] ){
            if(length(chains)==1 & length(dim(t))==2){
              tmp = t[iter_range[1]:iter_range[2],,drop=F]
            }else{
              tmp = t[iter_range[1]:iter_range[2],,ch,drop=F]
            }
          }
        }else{
          tmp = t[,,ch,drop=F]
        }
        if(is.vector(tmp)){
          tmp
        }else{
          apply(tmp, 2, mean)
        }
        
      }, mc.cores = cores)
      df = do.call(rbind, df)
      rownames(df) = paste0("chain", chains)
      
    }
  }
  return(df)
}

#' GetAllParamMeanChain
#'
#' @description
#' Get the average value of all parameters on each chain
#'
#' @param mamba mamba object
#' @param murp murp result
#' @param param parameter name
#' @param iter_range Iteration range, default is c(1, 3000)
#' @param cores threads
#' @export
#'
GetAllParamMeanChain <- function(object,
                                 iter_range = NULL,
                                 cores = 1,
                                 aspect = "pse",
                                 save = F){
  
  # aspect = "pse"
  # iter_range = rep(getParams(ct, "pse_optim_iterations"),2)
  # whichPos = 'all'
  # iter_range = iter_range
  # cores = 1
  # param = "m_t"
  
  mamba = object@OptimResult[[aspect]]
  murp = object@MURP
  all_params <- mamba@dimnames[[2]]
  param_names <- rownames(as.data.frame.numeric(table(str_split_fixed(all_params,"\\[",2)[,1])))
  param_number <- as.data.frame.numeric(table(str_split_fixed(all_params,"\\[",2)[,1]))
  
  param_meanvalue <- lapply(param_names, function(param){
    print(param)
    GetParamMeanChain(mamba = mamba,
                      murp = murp,
                      param = param,
                      whichPos = 'all',
                      iter_range = iter_range,
                      cores = 1)
  })
  names(param_meanvalue) <- param_names
  
  if(save){
    if(aspect=="pse"){
      save(param_meanvalue,
           file = paste0("2_pseudotime/param_meanvalue_", iter_range[1],"_",iter_range[2], ".rda"))
    }
    if(aspect=="track"){
      save(param_meanvalue,
           file = paste0("3_tracking/param_meanvalue_", iter_range[1],"_",iter_range[2], ".rda"))
    }
  }
  return(param_meanvalue)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                                    pse
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' MapMURPLabelToAll
#'
#' @description
#' mapping MURP label to all cells
#'
#' @param vecc
#' @param orig
#'
MapMURPLabelToAll <- function(vecc,
                              orig){
  
  vecc_names = paste0("T[",1:length(vecc),"]")
  murp_t = data.frame(row.names = vecc_names,
                      names = vecc_names,
                      value = vecc)
  
  if(is.null(names(orig))){
    orig_names = paste0("c",1:length(orig))
  }else{
    orig_names = names(orig)
  }
  cell_t = data.frame(row.names = orig_names,
                      cell_id = orig_names,
                      names =  paste0("T[",orig,"]"))
  cell_t = merge(murp_t, cell_t, by = "names")
  rownames(cell_t) = cell_t$cell_id
  cell_t = cell_t[orig_names,]
  
  orig_value = cell_t$value
  names(orig_value) = cell_t$cell_id
  return(orig_value)
}

#' GetMURPMapLabel
#' @description
#' map the existing attribute labels about cells to murp
#'
#' Arguments:
#' @param object MGPfact object
#' @param labels
#'
#' @export
#'
GetMURPMapLabel <- function(object, labels = NULL){
  
  ind = which(labels %in% colnames(object@MetaData)) %>% length
  if(ind==0) labels = NULL
  len = length(labels)
  
  ######### 1. get murp_cellinfo
  if(!is.null(labels)){
    mcl = paste0("T[", object@MetaData$murp_cluster, "]")
    tmp <- lapply(labels, function(lab){
      # cat(lab, "\n")
      tab = table(object@MetaData[,lab],mcl) %>% as.matrix
      tab_prop = apply(tab, 2, function(x){ x/sum(x) })
      if(is.null(dim(tab_prop))) {
        tmp = rep(rownames(tab),length(tab_prop))
        names(tmp) = names(tab_prop)
        tmp
      }else{
        apply(tab_prop, 2, function(x){rownames(tab_prop)[which.max(x)]})
      }
      
    })
    tmp <- do.call(cbind, tmp) %>% as.data.frame
    colnames(tmp) <- d2_label
    object = AddMURPMetadata(object, tmp)
  }
  object = assignSettings(object,"label",labels)
  return(object)
}

#' GetMURPInfo
#'
#' @param object
#'
GetMURPInfo <- function(object){
  r = object@MURP$murp_cellinfo
  return(r)
}

#' AddMURPMetadata
#'
#' @param object
#' @param df
#'
AddMURPMetadata <- function(object, df){
  
  meta = GetMURPInfo(object)
  if(is.null(meta)){
    meta = data.frame(row.names = paste0("T[",1:nrow(object@MURP$Recommended_K_cl$centers),"]"),
                      murp_cluster = 1:object@MURP$Recommended_K)
  }
  i = which(colnames(meta) %nin% colnames(df))
  cols = colnames(meta)[i]
  meta = cbind(meta[,cols,drop=F], df[rownames(meta),,drop=F])
  object@MURP$murp_cellinfo = meta
  return(object)
}

#' AddMURPMetadata
#'
#' @param object
#' @param df
#'
AddMetadata <- function(object, df){
  
  meta = object@MetaData
  i = which(colnames(meta) %nin% colnames(df))
  cols = colnames(meta)[i]
  meta = cbind(meta[,cols,drop=F], df)
  
  object@MetaData = meta
  return(object)
}

#' GetPseSdf
#'
#' @description
#' get all param info from optim result
#'
#' @param object
#' @param param_meanvalue
#' @param unified_direction
#' @param rm_adjust_chain
#' @param rev
#'
GetPseSdf <- function(object,
                      param_meanvalue = NULL,
                      unified_direction = FALSE,
                      rm_adjust_chain = FALSE,
                      rev = FALSE){
  
  # param_meanvalue = GetAllParamMeanChain(object = ct, aspect = "pse",
  #                                        iter_range = rep(getParams(ct, "pse_optim_iterations"),2))
  # object = ct
  # unified_direction = FALSE
  # rm_adjust_chain = TRUE
  # rev = FALSE
  
  if(is.null(param_meanvalue)){
    iter = getParams(ct,"pse_optim_iterations")
    # load(paste0("2_pseudotime/param_meanvalue_",iter,"_",iter,".rda"))
    param_meanvalue <- GetAllParamMeanChain(object = object, aspect = "pse",
                                            iter_range = rep(getParams(object, "pse_optim_iterations"),2))
  }
  if(!is.null(labels)){
    labels = getParams(object,"label")
  }
  
  mamba = object@OptimResult$pse
  murp = object@MURP
  metadata = object@MetaData
  chains = mamba@chains
  
  P <- getParams(object, "murp_number")
  L <- getParams(object, "trajectory_number")
  Q <- getParams(object, "murp_pc_number")
  
  ## 哪些chain是最终保留下来的
  t_pred =mamba@t_pred
  if(object@Settings@settings$chains_number>1){
    if(t_pred$adjust){
      rev_ch <- t_pred$adjust_chain
    }else{
      rev_ch <- NULL
    }
    if(rm_adjust_chain){
      filter_ch <- setdiff(t_pred$keep_chain,rev_ch)
    }else{
      filter_ch <- t_pred$keep_chain
    }
  }else{
    filter_ch <- t_pred$keep_chain
  }
  cat("filter_chain", filter_ch, "\n")
  
  param_tb <- param_meanvalue$Tb[filter_ch,,drop=FALSE]
  param_T <- param_meanvalue$T[filter_ch,,drop=FALSE]
  param_X <- param_meanvalue$X[filter_ch,,drop=FALSE]
  param_lambda1 <- param_meanvalue$lambda1[filter_ch,, drop=FALSE]
  param_lambda2 <- param_meanvalue$lambda2[filter_ch,, drop=FALSE]
  param_lambda3 <- param_meanvalue$lambda3[filter_ch,, drop=FALSE]
  param_mt <- param_meanvalue$m_t[filter_ch,, drop=FALSE] %>% mean
  param_s2t <- param_meanvalue$s2_t[filter_ch,, drop=FALSE] %>% mean
  param_s2x <- param_meanvalue$s2_x[filter_ch,, drop=FALSE] %>% mean
  param_rho <- param_meanvalue$rho[filter_ch,, drop=FALSE] %>% mean
  
  ### 用P确定C
  x <- param_meanvalue$C[filter_ch, , drop=FALSE]
  param_p <- param_meanvalue$p[filter_ch, , drop=FALSE]
  param_C <- ifelse(param_p>0.5, 1, 2)
  colnames(param_C) <- colnames(x)
  
  ### 方向的统一性，需要同时翻转T和Tb
  if(unified_direction & !rm_adjust_chain){
    if(length(which(rev_ch %in% filter_ch))!=0 ){
      inter = intersect(rev_ch, filter_ch)
      param_tb[inter, ] = 1 - param_tb[inter, ]
      param_T[inter, ] = 1 - param_T[inter, ]
    }
  }
  
  ### 方向翻转
  if(rev){
    param_T <- 1- param_T
    param_tb <- 1 - param_tb
  }
  
  ### 取均值
  T <- apply(param_T, 2, mean)
  X <- apply(param_X, 2, mean)
  
  ## 判断是不是单链
  if(nrow(param_tb)==1){
    
    tb_list = as.list( param_tb )
    l1_list = param_lambda1[1,,drop=F]
    l2_list = param_lambda2[1,,drop=F]
    l3_list = param_lambda3[1,,drop=F]
    tmp = matrix(param_C[1,], nr=L)
    
    c_combine = tmp
    rownames(c_combine) = paste0("L_", 1:L)
    
    c0_combine <- c_combine
    for(i in 1:L){
      ind = which(T < mean(tb_list[[i]]))
      if(length(ind)!=0)
        c0_combine[i,ind] = 0
    }
    
  }else{
    
    ### 1. 对10个chain的Tb向量化并聚类，得到均值Tb
    param_tb_2 <- param_tb
    colnames(param_tb_2) <- paste0("L",1:L)
    tb_chain <- as.vector(param_tb_2)
    names(tb_chain) <- unlist(lapply(colnames(param_tb_2),
                                     function(x) paste0(rownames(param_tb_2),"_",x)))
    tb_pam <- cluster::pam(tb_chain, L, metric = "manhattan")
    pdf("2_pseudotime/2.2_binarytree/tb_cluster.pdf")
    plot(x = tb_chain, y = tb_chain, col = tb_pam$clustering)
    dev.off()
    # Tb <- tapply(tb_chain, tb_pam$clustering, mean)
    # names(Tb) <- paste0("Tb", 1:L)
    
    ### 2. 准备一个漂亮的C
    C_list <- lapply(filter_ch, function(ch){
      tmp = matrix(param_C[ch, ], nr=L)
      rownames(tmp) = c(paste0(ch,"_L",1:L))
      tmp
    })
    names(C_list) <- c(filter_ch)
    c <- do.call(rbind, C_list)
    
    ### 3. 根据聚类对C | Tb | lambda1 | lambda2 | lambda3分组
    cluster <- tb_pam$clustering
    split_c <- split(names(cluster), as.factor(cluster))
    group_list <- lapply(split_c, function(x){
      if(length(x)==1){
        tmp = matrix(c[x,],nr=1)
        dimnames(tmp) = list(x,colnames(c))
        tmp
      }else{
        c[x,]
      }
    })
    
    param_tb_2 <- param_tb
    colnames(param_tb_2) <- paste0("L",1:L)
    tb <- melt(param_tb_2)
    tb <- data.frame(row.names = paste0(tb[,1],"_",tb[,2]), tb = tb[,3])
    tb_list <- split(tb[names(cluster),], as.factor(cluster))
    
    param_l1_2 <- param_lambda1
    colnames(param_l1_2) <- paste0("L",1:L)
    l1 <- melt(param_l1_2)
    l1 <- data.frame(row.names = paste0(l1[,1],"_",l1[,2]), l1 = l1[,3])
    l1_list <- split(l1[names(cluster),], as.factor(cluster))
    
    param_l2_2 <- param_lambda2
    colnames(param_l2_2) <- paste0("L",1:L)
    l2 <- melt(param_l2_2)
    l2 <- data.frame(row.names = paste0(l2[,1],"_",l2[,2]), l2 = l2[,3])
    l2_list <- split(l2[names(cluster),], as.factor(cluster))
    
    param_l3_2 <- param_lambda3
    colnames(param_l3_2) <- paste0("L",1:L)
    l3 <- melt(param_l3_2)
    l3 <- data.frame(row.names = paste0(l3[,1],"_",l3[,2]), l3 = l3[,3])
    l3_list <- split(l3[names(cluster),], as.factor(cluster))
    
    ### 5. 翻转
    group_rev_list <- lapply(group_list, function(xt){
      xi = as.matrix(xt)
      x_ind = which(apply(xi, 1, sd) != 0)
      x = xi[x_ind, ]
      if(is.null(dim(x))){
        x = as.matrix(x) %>% t
        max_i = which.max(apply(cor(t(x)),1,mean, na.rm = TRUE))
        x1 = x[max_i,]
      }else if(nrow(x)==0){
        x1 = NULL
      }else{
        max_i = which.max(apply(cor(t(x)),1,mean, na.rm = TRUE))
        x1 = x[max_i,]
      }
      if(is.null(x1)){
        x_rev = xi
      }else{
        x_rev = apply(x, 1, function(y){
          y1 = y;
          y2 = y; y2[which(y==1)]=2; y2[which(y==2)]=1
          c1 = cor(x1, y1); c2 = cor(x1, y2)
          if(c2>c1) y2 else y1
        }) %>% t
      }
      return(x_rev)
    })
    
    ### 6. 合并
    c_combine = lapply(group_rev_list, function(tb_C){
      apply(tb_C, 2, function(x){
        t = table(x)
        ind = which.max(t)
        as.numeric(names(t)[ind])
      })
    })
    c_combine = do.call(rbind, c_combine)
    rownames(c_combine) = paste0("L_", 1:L)
    
    ### 7. 生成合并的c0_combine
    c0_combine <- c_combine
    for(i in 1:L){
      ind = which(T < mean(tb_list[[i]]))
      if(length(ind)!=0)
        c0_combine[i,ind] = 0
    }
    
  }
  
  
  ###############################################
  #####             create sdf              #####
  ###############################################
  ### 1. 准备sdf
  ord <- 1:L
  Tb <- sapply(tb_list, mean)
  lambda1 <- sapply(l1_list, mean)
  lambda2 <- sapply(l2_list, mean)
  lambda3 <- sapply(l3_list, mean)
  mt <- param_mt
  s2t <- param_s2t
  s2x <- param_s2x
  rho <- param_rho
  
  C <- t(c_combine)
  C0 <- t(c0_combine)
  
  ord <- order(Tb, decreasing = FALSE)
  Tb <- Tb[ord]
  C <- C[,ord]
  C0 <- C0[,ord]
  lambda1 <- lambda1[ord]
  lambda2 <- lambda2[ord]
  lambda3 <- lambda3[ord]
  
  sdf <- data.frame(T = T,
                    C0,
                    C,
                    matrix(rep(Tb, each=P),nr=P),
                    object@MURP$centersPCA$x[,1:getParams(ct, "murp_pc_number")],
                    X = X,
                    matrix(rep(lambda1, each=P),nr=P),
                    matrix(rep(lambda2, each=P),nr=P),
                    matrix(rep(lambda3, each=P),nr=P),
                    mt = rep(mt, P),
                    s2t = rep(s2t, P),
                    s2x = rep(s2x, P),
                    rho = rep(rho, P))
  colnames(sdf) = c("T",
                    paste0("C0_", 1:L),
                    paste0("C_", 1:L),
                    paste0("Tb_", 1:L),
                    paste0("PC_", 1:Q),
                    "X",
                    paste0("lambda1_", 1:L),
                    paste0("lambda2_", 1:L),
                    paste0("lambda3_", 1:L),
                    "mt","s2t","s2x","rho")
  
  sdf = sdf[1:nrow(sdf),1:ncol(sdf)]
  sdf$unified_direction = unified_direction
  
  ### 2. 准备sdf_orig
  orig_sdf <- data.frame(row.names = rownames(object@assay$data_matrix),
                         names = rownames(object@assay$data_matrix),
                         T = MapMURPLabelToAll(vecc = sdf$T, orig = murp$Recommended_K_cl$cluster) )
  for(l in 1:L){
    c0 = paste0("C0_",l)
    c = paste0("C_",l)
    orig_sdf$a = MapMURPLabelToAll(vecc = sdf[,c0], orig = murp$Recommended_K_cl$cluster)
    orig_sdf$b = MapMURPLabelToAll(vecc = sdf[,c], orig = murp$Recommended_K_cl$cluster)
    colnames(orig_sdf)[(ncol(orig_sdf)-1):ncol(orig_sdf)] = c(c0,c)
  }
  
  ### 3. save
  # save(sdf, file = "2_pseudotime/sdf.rda")
  # save(orig_sdf, file = "2_pseudotime/orig_sdf.rda")
  
  ### 4. save to object
  # save(sdf, file = "~/sdf.rda")
  object = AddMURPMetadata(object, sdf)
  object = AddMetadata(object, orig_sdf)
  
  command = GetCommand()
  object@Command$pse$GetPseSdf = command
  
  return(object)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                                    tree
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' split_tb
#'
#' cut T (pseudotime) according to tb (bifurcation time)
#'
#' @param t pseudotime of murps
#' @param tb bifurcation time of different trajectory
#'
#' @return a matrix with Gene Number * Cell Number
#'
split_tb <- function(t, tb){
  pos = lapply(1:(length(tb)+1), function(i){
    if(i==1){
      which(t < tb[i])
    }else if(i==(length(tb)+1)){
      which(t >= tb[i-1])
    }else {
      which((t >= tb[i-1]) & (t < tb[i]) )
    }
  })
  return(pos)
}


#' split_binary_tree
#'
#' 对拆分之后的前后分组再根据分叉情况分组
#'
#' @param df
#' @param sep trajectory name
#'
#' @return list
#'
split_binary_tree <- function(df, sep = "L1"){
  
  require(dplyr)
  require(tibble)
  
  final_result <- lapply(1:length(df), function(i){
    print(i)
    if(i == 1){
      x = df[[1]] %>%
        rownames_to_column("c_names") %>%
        dplyr::select(c_names, T) %>%
        dplyr::group_split()
      names(x) = "g1"
    }else{
      c = sep
      x = df[[i]] %>% rownames_to_column("c_names") %>%
        dplyr::select(c_names, T, c) %>%
        dplyr::group_by(across(all_of(c)) ) %>%
        dplyr::group_split()
      tmp = do.call(rbind,lapply(x, function(y)y[1,3:ncol(y)]))
      if(!is.null(tmp)){
        names(x) = apply(tmp,1,function(x)paste0("g",i,"_",paste0(x,collapse = "_")))
      }
    }
    x
  })
  
  return(final_result)
}

#' split_tb_tree
#'
#' 根据Tb进行分段之后，对每一个分段再进行不同分支的分组
#' 每一个分段的分组依据是当前tb以及之前tb的L的排列组合
#'
#' @param df_list
#'
#' @return list
#'
#' @export
#'
split_tb_tree <- function(df_list){
  
  require(dplyr)
  require(tibble)
  
  tree <- lapply(1:length(df_list), function(i){
    print(i)
    df = df_list[[i]]
    
    if(nrow(df)==0){
      x = NULL
    }else if(i == 1){
      x = df_list[[1]] %>%
        rownames_to_column("c_names") %>%
        dplyr::select(c_names, T) %>%
        dplyr::group_split()
      names(x) = "g1"
    }else{
      # c = paste0("L",1:(i-1))
      c = tail(colnames(df), ncol(df)-1)
      x = df_list[[i]] %>% rownames_to_column("c_names") %>%
        dplyr::select(c_names, T, c) %>%
        dplyr::group_by(dplyr::across(all_of(c)) ) %>%
        dplyr::group_split()
      tmp = do.call(rbind,lapply(x, function(y)y[1,3:ncol(y)]))
      names(x) = apply(tmp,1,function(x) paste0("g",i,"_",paste0(x,collapse = "_")))
    }
    x
  })
  
  return(tree)
}


#' adj_Tb_tree / adj_binary_tree
#'
#' @description
#' adjacency matrix of binary tree
#'
#' @param tb_group_list
#' @param P node number
#' @param point_names
#'
#' @return matrix
#'
#' @export
#'
adj_tb_tree <- function(tb_group_list, P, point_names = NULL){
  
  require(stringr)
  require(tibble)
  require(dplyr)
  
  Group = tb_group_list
  
  # construct adj matrix
  adj_matrix <- matrix(0, nr = P, nc = P)
  if(is.null(point_names)){
    rownames(adj_matrix) <- paste0("c",1:P)
    colnames(adj_matrix) <- paste0("c",1:P)
  }else{
    rownames(adj_matrix) <- point_names
    colnames(adj_matrix) <- point_names
  }
  
  
  # 在按照tb分割后的每一段分别进行循环
  for (g in 1:length(Group)){
    
    G = Group[[g]]
    
    # cat(paste0("G",g,"\n"))
    if(is.null(G)){
      # cat("G",g,"is null \n")
      Group[[g]] = NA
      next
    }else if(length(G)==0){
      # cat("G",g,"is null \n")
      Group[[g]] = NA
      next
    }else if(nrow(G[[1]])==0){
      # cat("G",g,"is null \n")
      Group[[g]] = NA
      next
    }else {
      # cat("G",g,"is not null \n")
    }
    # if(nrow(as.data.frame(G))==0){
    #   cat("G",g,"is null \n")
    #   G = NULL
    #   next
    # }
    
    # 1. 在这一段中按照不同的簇分组，组内按照pseudotime连接
    # cat("Loop in every layer \n")
    for(d in G){
      if(nrow(d)==0) next;
      n = d$c_names
      # cat(n, "\n")
      for(i in 1:nrow(d) ){
        #print(i)
        if(i==nrow(d)) break;
        adj_matrix[n[i], n[i+1]] = 1
        #adj_matrix[n[i+1], n[i]] = 1
      }
    }
    
    # 2. 为这一段中的每一个小组找到上家
    # g=1; g=2; g=...
    # cat("Connect between layers \n")
    if(g==1){
      next;
    }else if(g==2){
      
      G2 = Group[[1]]
      
      # check if all tb < T
      # if(nrow(G2[[1]])==0){
      if(is.na(G2)|is.null(G2)){
        tmp = do.call(rbind,G)
        tmp = tmp[order(tmp$T),]
        min.name = tmp[1,]$c_names
        for(n in names(G)){
          x_c = G[[n]][1,]$"c_names"
          if(min.name!=x_c){
            adj_matrix[min.name, x_c] = 1
          }
        }
      }else{
        y_c = tail(G2[[1]],1)$c_names
        for(n in names(G)){
          x_c = G[[n]][1,]$"c_names"
          adj_matrix[y_c, x_c] = 1
        }
      }
      
    }else{
      
      G_n = stringr::str_split_fixed(names(G),"_",2)
      G_n[,2] = paste0("_", G_n[,2])
      
      for(n in G_n[,2]){
        # cat("Group name: ", n, "\n")
        for(g2 in (g-1):1){
          # cat("----Last Group: ", g2, "\n")
          G2 = Group[[g2]]
          aim = str_sub(n, 1, 2*(g2-1))
          ind = grep(aim, names(G2))
          if(length(ind)!=0){
            # cat("----Find it!! \n")
            x_n = paste0(G_n[1,1],n)
            y_n = names(G2)[ind]
            x_c = G[[x_n]][1,]$"c_names"
            y_c = tail(G2[[y_n]],1)$c_names
            
            adj_matrix[y_c, x_c]=1
            break;
          }
          
          # 没找到上家，直接连接到第一个点
          if(g2==1){
            # cat("----No find. \n ")
            # 生成原始的df
            # dd = lapply(tb_group_list, function(x) {
            #   if(!is.null(x)){
            #     do.call(rbind,x) %>% select(c_names, T) %>% data.frame
            #   } })
            # orig_df = do.call(rbind,dd)
            # start = orig_df[order(orig_df$T),][1,"c_names"]
            # x_n = paste0(G_n[1,1],n)
            # x_c = G[[x_n]][1,]$"c_names"
            # adj_matrix[start, x_c] = 1
          }
          
        }
      }
      ##
    }
  }
  return(adj_matrix)
}

#' GetBinTree
#'
#' @description
#' Get the deconstructed multiple trajectories,
#' and display them through a single-fork binary tree
#'
#' @param object MGPfact object
#' @param pse_sdf the data frame obtained by the function GetPseSdf
#'
#' @export
#'
GetBinTree <- function(object,
                       save = TRUE){
  
  sdf = GetMURPInfo(object)
  sdf = sdf[order(sdf$T, decreasing = FALSE), ]
  P = getParams(object, "murp_number")
  L = getParams(object, "trajectory_number")
  Q = getParams(object, "murp_pc_number")
  
  mamba = object@OptimResult$pse
  murp = object@MURP
  metadata = object@MetaData
  chains = mamba@chains
  
  ## 0. 获取bin_f和Tb
  bin_f <- sdf[, c("T", paste0("C0_",1:L))]
  Tb <- sdf[1,paste0("Tb_",1:L)] %>% as.matrix %>% as.vector
  names(Tb) <- paste0("Tb_", 1:L)
  
  ## 1. Segment bin_f according to different differentiation situations (tb)
  ## length(df_list) = length(Tb)
  tb_cut_list <- lapply(1:L, function(i){
    tb = paste0("Tb_", i)
    # cat("tb: ", tb, "\n")
    pos = split_tb(bin_f$T, Tb[tb])
    lapply(pos, function(x) bin_f[x,])
  })
  names(tb_cut_list) <- paste0("Tb_", 1:L)
  
  ## 2. group the (bin_f) of each trajectory according to the fork (T>tb)
  ## length(group_list) = length(Tb)
  tb_group_list <- lapply(1:L, function(i){
    tb = paste0("Tb_", i)
    l = paste0("C0_",i)
    # cat("tb: ", tb, " C: ", l, "\n")
    split_binary_tree(tb_cut_list[[tb]], sep = l)
  })
  
  ## 3. get the adjacency matrix and get the binary tree
  binary_tree_list <- lapply(tb_group_list, function(tb_group){
    # print(paste0("tree:",i))
    graph_from_adjacency_matrix(adj_tb_tree(tb_group, P, point_names =  rownames(sdf)) )
  })
  
  ## 4. get layout of each tree
  layout_list <- lapply(1:L, function(i){
    layout_as_tree(binary_tree_list[[i]],
                   root = rownames(bin_f)[1],
                   circular = FALSE,
                   mode = 'out',
                   flip.y = TRUE)
  })
  
  ## 5. get each bin tree
  centers <- murp$Recommended_K_cl$centers
  rownames(centers) <- paste0("T[",1:nrow(centers),"]")
  centers <- centers[rownames(bin_f),]
  bintree_all <- lapply(1:L, function(i){
    tree <- binary_tree_list[[i]]
    layout <- layout_list[[i]]
    edge <- as.data.frame(get.edgelist(tree), stringsAsFactors = FALSE)
    vertex <- data.frame( name = rownames(bin_f),
                          label = str_sub(rownames(bin_f), start = 3, end = -2),
                          c_label = bin_f[ ,paste0("C0_",i)],
                          pse = bin_f[ ,"T"],
                          centers,
                          sdf)
    newGraph <- graph_from_data_frame(edge, directed = TRUE, vertices = vertex)
    list(edge = edge, vertex = vertex, layout = layout, graph = newGraph)
  })
  names(bintree_all) <- paste0("trajectory",1:L)
  if(save){
    save(bintree_all, file = "2_pseudotime/bintree_all.rda")
  }
  object@Tree$bintree = bintree_all
  
  return(object)
}


#' GetTbTree
#'
#' @description
#' Get a multi-forked binary-tree by merging different trajectories
#'
#' @param object MGPfact object
#' @param pse_sdf the data frame obtained by the function GetPseSdf
#' @param save Logical value, whether to save the result
#'
#' @export
#'
GetTbTreeDrq <- function(object,
                         save = TRUE){
  
  ### prepare
  sdf = GetMURPInfo(object)
  sdf = sdf[order(sdf$T, decreasing = FALSE), ]
  P = getParams(object, "murp_number")
  L = getParams(object, "trajectory_number")
  Q = getParams(object, "murp_pc_number")
  mamba = object@OptimResult$pse
  murp = object@MURP
  metadata = object@MetaData
  chains = mamba@chains
  
  ### function
  make_adj <- function(names){
    tmp = matrix(0, length(names), length(names)-1)
    diag(tmp) = 1
    tmp = cbind(matrix(0,length(names),1), tmp)
    dimnames(tmp) = list(names, names)
    return(tmp)
  }
  
  ###########################################
  ### 0. 获取tb_f和Tb
  tb_f <- sdf[, c("T", paste0("C0_",1:L))]
  Tb <- sdf[1,paste0("Tb_",1:L)] %>% as.matrix %>% as.vector
  names(Tb) <- colnames(sdf)[(2+L*2):(1+3*L)]
  
  ### 1. make a null matrix
  adj = matrix(0, P, P)
  dimnames(adj) = list(rownames(sdf), rownames(sdf))
  
  ## 2. 寻找跟节点在哪个TB处
  # for(i in 1:L){ # 选择第一个大于>min(T)的tb，相当于前面的都要剔除
  #   cat("root:", i, "\n")
  #   tb = Tb[i]
  #   if(tb<min(tb_f$T)){
  #     next
  #   } else{
  #     names = rownames(tb_f)[which(tb_f$T<=tb)] # root 节点前面的连接
  #     adj[names, names] = make_adj(names)
  #     end = tail(names, 1) ### end point
  #     root = i
  #     break
  #   }
  # }
  
  index = which(Tb > min(tb_f$T)) # 所有Tb都要保留
  if(length(index)==0){ # 没有tb > min(T), 一条直线
    adj <- matrix(0, nr = P, nc = P)
    rownames(adj) <- paste0("T[",1:P,"]")
    colnames(adj) <- paste0("T[",1:P,"]")
    for(i in 1:(nrow(tb_f)-1) ){
      adj[rownames(tb_f)[i],rownames(tb_f)[i+1]] = 1
    }
  }else if(length(index)==1){ # 只有1个tb > min(T), 一个二叉树
    pos = split_tb(tb_f$T, Tb[index])
    tb_cut = lapply(pos, function(x) tb_f[x,])
    # tb_group = split_binary_tree(tb_cut, sep = paste0("L",index))
    tb_group = split_binary_tree(tb_cut, sep = paste0("C0", str_sub(names(index), start= 3) ))
    adj = adj_tb_tree(tb_group, P, point_names =  rownames(sdf))
  }else{ # 两个及以上，无论如何以第一个tb作为第一个root
    tb = Tb[1]
    names = rownames(tb_f)[which(tb_f$T<=tb)] # root 节点前面的连接
    if(length(names)==0){
      names = rownames(tb_f)[1]
      first_t =  tb_f$T[1]
      end = rownames(tb_f)[1]
    }else{
      adj[names, names] = make_adj(names)
      end = tail(names, 1) ### end point
      first_t =  Tb[1]
    }
    # if(length(names)!=0){
    #   adj[names, names] = make_adj(names)
    #   end = tail(names, 1) ### end point
    #   first_t =  Tb[1]
    # }else{
    #   end = rownames(tb_f)[1]
    #   first_t =  tb_f$T[1]
    # }
    root = 1
    
    ## 3. 在root之后的分成c1, c2
    root_c = paste0("C0_",root)
    root_ind1 = rownames(tb_f)[which(tb_f$T > first_t & tb_f[,root_c]==1)]
    root_ind2 = rownames(tb_f)[which(tb_f$T > first_t & tb_f[,root_c]==2)]
    ind1_list = list( list(ind = root_ind1, end = end) ) # 把应该连接的上一个根节点也存起来
    ind2_list = list( list(ind = root_ind2, end = end) )
    adj[root_ind1, root_ind1] = make_adj(root_ind1) ## 第一个二叉树建立好
    adj[root_ind2, root_ind2] = make_adj(root_ind2)
    adj[end, root_ind1[1]] = 1
    adj[end, root_ind2[1]] = 1
    
    for(i in (root+1):L){
      
      # cat("adj", i, "\n")
      tb = Tb[i]
      c = paste0("C0_",i)
      new_ind1_list = list() # 创建空list准备放入
      new_ind2_list = list()
      
      # 上一个L中c1的后面
      for(j in 1:length(ind1_list)){
        
        # cat(j, "\n")
        ind = ind1_list[[j]]$ind
        if(is.null(ind)) next
        tmp = tb_f[ind,]
        median = rownames(tmp)[which(tmp$T<=tb)] # 中间的
        if(length(median)==0){
          new_end = end
        }else{
          new_end = tail(median,1) # 更新本次tb的end
        }
        
        upper = rownames(tmp)[which(tmp$T > tb)]
        adj[upper, ] = 0
        adj[ ,upper] = 0
        # adj[new_end, ] = 0
        
        # if(length(new_end)==0) new_end = end
        tmp2 = tmp[upper,,drop=FALSE]
        if(nrow(tmp2)==0) next
        new_id1 = rownames(tmp2)[which(tmp2[, c]==1)]
        new_id2 = rownames(tmp2)[which(tmp2[, c]==2)]
        # adj[new_end, new_id1[1]] = 1
        # adj[new_end, new_id2[1]] = 1
        
        # 把对应的新生成的c1,c2放在新的list里面
        if(length(new_id1)!=0){
          adj[new_end, new_id1[1]] = 1
          adj[new_id1,new_id1] = make_adj(new_id1)
          new_ind1_list = append(new_ind1_list,  list(list(ind = new_id1, end = new_end)) )
        }
        if(length(new_id2)!=0){
          adj[new_end, new_id2[1]] = 1
          adj[new_id2,new_id2] = make_adj(new_id2)
          new_ind2_list = append(new_ind2_list, list(list(ind = new_id2, end = new_end)) )
        }
      }
      
      # 上一个L中c2的后面
      for(j in 1:length(ind2_list)){
        ind = ind2_list[[j]]$ind
        if(is.null(ind)) next
        tmp = tb_f[ind,]
        median = rownames(tmp)[which(tmp$T<=tb)] # 中间的
        if(length(median)==0){
          new_end = end
        }else{
          new_end = tail(median,1) # 更新本次tb的end
        }
        
        upper = rownames(tmp)[which(tmp$T > tb)]
        adj[upper, ] = 0; adj[ ,upper] = 0;
        # adj[new_end, ] = 0
        
        # if(length(new_end)==0) new_end = end
        tmp2 = tmp[upper,]
        if(nrow(tmp2)==0) next
        new_id1 = rownames(tmp2)[which(tmp2[, c]==1)]
        new_id2 = rownames(tmp2)[which(tmp2[, c]==2)]
        # adj[new_end, new_id1[1]] = 1
        # adj[new_end, new_id2[1]] = 1
        
        # 把对应的新生成的c1,c2放在新的list里面
        if(length(new_id1)!=0){
          adj[new_end, new_id1[1]] = 1
          adj[new_id1,new_id1] = make_adj(new_id1)
          new_ind1_list = append(new_ind1_list,  list(list(ind = new_id1, end = new_end)) )
        }
        if(length(new_id2)!=0){
          adj[new_end, new_id2[1]] = 1
          adj[new_id2,new_id2] = make_adj(new_id2)
          new_ind2_list = append(new_ind2_list, list(list(ind = new_id2, end = new_end)) )
        }
      }
      ind1_list = new_ind1_list
      ind2_list = new_ind2_list
    }
  }
  
  y = adj
  # x = c(names, median, new_id1, new_id2)
  # y = adj[x,x]
  tmp <- graph_from_adjacency_matrix(y)
  edge <- as.data.frame(get.edgelist(tmp), stringsAsFactors = FALSE)
  vertex = data.frame(row.names = rownames(y),
                      name = rownames(y))
  g <- graph_from_data_frame(edge, directed = TRUE, vertices = vertex)
  layout <- layout_as_tree(g,
                           root = rownames(y)[1],
                           circular = FALSE,
                           mode = 'out',
                           flip.y = TRUE)
  ggraph(g, layout = layout) +
    geom_edge_diagonal(alpha = 0.7, width = 0.5, check_overlap = FALSE) +
    geom_node_point(color = "blue",size = 6, alpha = 0.6) +
    geom_node_text(aes(label = name), size = 4) +
    coord_flip() +
    scale_y_reverse() +
    labs(title = paste0("L",i), color = "Bif") +
    rj.graph.ftheme
  
  tb_adj_matrix = adj
  
  ### 2. 画图准备, get tbtree list
  tmp <- graph_from_adjacency_matrix(tb_adj_matrix)
  edge <- as.data.frame(get.edgelist(tmp), stringsAsFactors = FALSE)
  rat <- table(murp$Recommended_K_cl$cluster)/length(murp$Recommended_K_cl$cluster)
  names(rat) <- paste0("T[", names(rat), "]")
  
  centers <- murp$Recommended_K_cl$centers
  rownames(centers) <- paste0("T[",1:nrow(centers),"]")
  centers <- centers[rownames(tb_f),]
  
  vertex = data.frame(row.names = rownames(tb_f),
                      name = rownames(tb_f),
                      label = str_sub(rownames(tb_f), start = 3, end = -2),
                      pse = tb_f$T,
                      grp = rownames(tb_f),
                      rat = as.vector(rat[rownames(tb_f)]),
                      sdf,
                      centers)
  
  ### 做一个C_all
  # vertex = tbtree$vertex
  # edge = tbtree$edge
  # layout = tbtree$layout
  tbx = unlist(vertex[1,paste0("Tb_",1:L),drop=TRUE])
  tbx = c(0,tbx,1)
  tmp = cut(vertex$T,tbx,labels = FALSE)
  tmp[which(tmp==1)] = "start"
  for(i in 2:length(tbx)){
    ind = which(tmp==i)
    l = i-1
    tmp[ind] = paste0("L",l,"_",vertex[ind,paste0("C_",l)])
  }
  vertex$C_all = tmp
  
  for(i in 1:L){ vertex[,paste0("C0_",i)] = as.factor(vertex[,paste0("C0_",i)]) }
  graph_tbtree <- graph_from_data_frame(edge, directed = TRUE, vertices = vertex)
  layout <- layout_as_tree(graph_tbtree,
                           root = rownames(tb_f)[1],
                           circular = FALSE,
                           mode = 'out',
                           flip.y = TRUE)
  tbtree <- list(edge = edge,
                 vertex = vertex,
                 graph = graph_tbtree,
                 layout = layout,
                 adj = tb_adj_matrix)
  
  if(save){
    save(tbtree, file = "2_pseudotime/tbtree.rda")
  }
  object@Tree$tbtree = tbtree
  
  return(object)
}

#' GetTbTree
#'
#' @description
#' get tbtree result of murp
#'
#' @param object MGPfact object
#' @param pse_sdf the data frame obtained by the function GetPseSdf
#'
#' @export
#'
GetTbTree <- function(object,
                      save = TRUE){
  
  ### judge sdf
  sdf = GetMURPInfo(object)
  sdf = sdf[order(sdf$T, decreasing = FALSE), ]
  P = getParams(object, "murp_number")
  L = getParams(object, "trajectory_number")
  Q = getParams(object, "murp_pc_number")
  
  mamba = object@OptimResult$pse
  murp = object@MURP
  metadata = object@MetaData
  chains = mamba@chains
  
  ### 0. 获取tb_f和Tb
  tb_f <- sdf[, c("T", paste0("C0_",1:L))]
  Tb <- sdf[1,paste0("Tb_",1:L)] %>% as.matrix %>% as.vector
  names(Tb) <- paste0("Tb_",1:L)
  # names(Tb) <- colnames(sdf)[(2+L*2):(1+3*L)]
  
  
  ### 1. 获得Tb-tree的邻接矩阵tb_adj_matrix
  index = which(Tb > min(tb_f$T))
  if(length(index)==0){
    itb = length(Tb)
    pos = split_tb(tb_f$T, Tb[itb])
    tb_cut = lapply(pos, function(x) tb_f[x,])
    tb_group = split_binary_tree(tb_cut, sep = paste0("C0_", itb))
    adj_matrix = adj_tb_tree(tb_group, P, point_names =  rownames(sdf))
    
  }else if(length(index)==1){
    pos = split_tb(tb_f$T, Tb[index])
    tb_cut = lapply(pos, function(x) tb_f[x,])
    # tb_group = split_binary_tree(tb_cut, sep = paste0("L",index))
    tb_group = split_binary_tree(tb_cut, sep = paste0("C0", str_sub(names(index), start= 3) ))
    adj_matrix = adj_tb_tree(tb_group, P, point_names =  rownames(sdf))
  }else{
    pos <- split_tb(tb_f$T, Tb)
    tb_cut_list <- lapply(pos, function(x) tb_f[x,])
    tb_group_list <- split_tb_tree(tb_cut_list)
    adj_matrix <- adj_tb_tree(tb_group_list, P, point_names =  rownames(sdf))
  }
  tb_adj_matrix = adj_matrix
  
  ### 2. 画图准备, get tbtree list
  tmp <- graph_from_adjacency_matrix(tb_adj_matrix)
  edge <- as.data.frame(get.edgelist(tmp), stringsAsFactors = FALSE)
  rat <- table(murp$Recommended_K_cl$cluster)/length(murp$Recommended_K_cl$cluster)
  names(rat) <- paste0("T[", names(rat), "]")
  
  centers <- murp$Recommended_K_cl$centers
  rownames(centers) <- paste0("T[",1:nrow(centers),"]")
  centers <- centers[rownames(tb_f),]
  
  vertex = data.frame(row.names = rownames(tb_f),
                      name = rownames(tb_f),
                      label = str_sub(rownames(tb_f), start = 3, end = -2),
                      pse = tb_f$T,
                      grp = rownames(tb_f),
                      rat = as.vector(rat[rownames(tb_f)]),
                      sdf,
                      centers)
  
  ### 做一个C_all
  # vertex = tbtree$vertex
  # edge = tbtree$edge
  # layout = tbtree$layout
  tbx = unlist(vertex[1,paste0("Tb_",1:L),drop=TRUE])
  
  if(1 %in% tbx){ tbx = c(0,tbx) }else{ tbx = c(0,tbx,1) }
  tmp = cut(vertex$T,tbx,labels = FALSE)
  tmp[which(tmp==1)] = "start"
  for(i in 2:length(tbx)){
    ind = which(tmp==i)
    l = i-1
    tmp[ind] = paste0("L",l,"_",vertex[ind,paste0("C_",l)])
  }
  vertex$C_all = tmp
  
  for(i in 1:L){ vertex[,paste0("C0_",i)] = as.factor(vertex[,paste0("C0_",i)]) }
  graph_tbtree <- graph_from_data_frame(edge, directed = TRUE, vertices = vertex)
  layout <- layout_as_tree(graph_tbtree,
                           root = rownames(tb_f)[1],
                           circular = FALSE,
                           mode = 'out',
                           flip.y = TRUE)
  tbtree <- list(edge = edge,
                 vertex = vertex,
                 graph = graph_tbtree,
                 layout = layout,
                 adj = tb_adj_matrix)
  
  if(save){
    save(tbtree, file = "2_pseudotime/tbtree.rda")
  }
  object@Tree$tbtree = tbtree
  return(object)
}

#' GetTbTreeAllpoint
#'
#' @description
#' add all point in tbtree
#'
#' @param object MGPfact object
#' @param tbtree tbtree result obtained by the function GetTbTree
#' @param pse_sdf the data frame obtained by the function GetPseSdf
#' @param labels some attribute about cells
#'
#' @export
#'
GetTbTreeAllpoint <- function(object,
                              save = TRUE,
                              labels = NULL){
  
  ### prepare
  sdf = GetMURPInfo(object)
  sdf = sdf[order(sdf$T, decreasing = FALSE), ]
  tbtree = GetTbTreeResult(object)
  P = getParams(object, "murp_number")
  L = getParams(object, "trajectory_number")
  Q = getParams(object, "murp_pc_number")
  
  mamba = object@OptimResult$pse
  murp = object@MURP
  metadata = object@MetaData
  chains = mamba@chains
  
  tb_f <- sdf[, c("T", paste0("C0_",1:L))]
  
  ### 1. 分别把每个簇的矩阵提取出来做成list
  sdata_cluster <- tapply(1:length(murp$Recommended_K_cl$cluster),
                          murp$Recommended_K_cl$cluster,
                          function(x,y){
                            exp = matrix(murp$rawdata[x,], ncol = ncol(murp$rawdata))
                            dimnames(exp) = list(rownames(murp$rawdata)[x], colnames(murp$rawdata))
                            exp })
  
  ### 2. 把每个中心点加入分别加入每个簇并命名
  sdata_cluster_centers = lapply(1:P, function(i){
    print(i)
    exp = sdata_cluster[[i]]
    centers = matrix(murp$Recommended_K_cl$centers[i,], ncol = ncol(murp$rawdata))
    tmp = rbind(exp, centers)
    # rownames(tmp) = c(rownames(exp), rownames(murp$Recommended_K_cl$centers)[i])
    ci = rownames(murp$Recommended_K_cl$centers)[i]
    rownames(tmp) = c(rownames(exp),
                      paste0("T[", ci, "]"))
    tmp
  })
  
  ### 3. 分别把每个簇的距离矩阵计算出来
  cores = 3
  dist_list <- lapply(1:P, function(i){
    exp = sdata_cluster_centers[[i]]
    dist_m = parallelDist::parDist(x = exp, method = "euclidean", diag = TRUE, upper = TRUE, threads = cores-2)
    dist_m
  })
  
  ### 4. 针对每个簇建立最小生成树
  max_v = max(unlist(dist_list))
  min_v = 1e-9
  mst_list <- lapply(1:P, function(i){
    print(i)
    dist_m = dist_list[[i]]
    if(max(dist_m)==0){
      dist_m = as.matrix(replace(dist_m, dist_m == 0, min_v))
    }else{
      dist_m = as.matrix(dist_m)
    }
    graph = graph.adjacency(dist_m, weighted=TRUE)
    mst = minimum.spanning.tree(graph)
    mst
  })
  
  ### 5. 获取所有簇的邻接矩阵
  adj_matrix_list <- lapply(1:P, function(i){
    print(i)
    mst = mst_list[[i]]
    adj = as.matrix(as_adjacency_matrix(mst))
    adj
  })
  
  ### 6. 得到所有点（包括中心点）的邻接矩阵
  adj_final <- as.matrix(do.call(Matrix::bdiag, adj_matrix_list))
  rownames(adj_final) <- unlist(lapply(adj_matrix_list, rownames))
  colnames(adj_final) <- rownames(adj_final)
  
  ### 7. 把tb-tree的结果加上
  cc = rownames(tbtree$adj)
  for(i in 1:P){
    for(j in 1:P){
      adj_final[cc[i],cc[j]] = tbtree$adj[cc[i],cc[j]]
    }
  }
  
  ### 8. 准备数据
  ### get edge
  tmp <- graph_from_adjacency_matrix(adj_final, mode = "directed")
  edge <- as.data.frame(get.edgelist(tmp), stringsAsFactors = FALSE)
  
  ### 每个簇的细胞个数/总细胞数
  cell_cluster_centers <- paste0("c", murp$Recommended_K_cl$cluster)
  rat <- table(murp$Recommended_K_cl$cluster)/length(murp$Recommended_K_cl$cluster)
  names(rat) <- paste0("T[", names(rat), "]")
  
  ### 点的gene属性
  centers <- murp$Recommended_K_cl$centers
  rownames(centers) <- paste0("T[",1:nrow(centers),"]")
  centers <- centers[rownames(tb_f),]
  sdata2 <- rbind(centers, murp$rawdata)
  
  ## 得到vertex
  vertex <- data.frame(row.names = c(rownames(tb_f), rownames(murp$rawdata)),
                       name = c(rownames(tb_f), rownames(murp$rawdata)),
                       pse = c(tb_f$T, tb_f[cell_cluster_centers,"T"]),
                       pse_na = c(tb_f$T, rep("",nrow(murp$rawdata))),
                       grp = c(rownames(tb_f), cell_cluster_centers),
                       grp_na = c(rownames(tb_f), rep("",nrow(murp$rawdata))),
                       rat = c(rat[rownames(tb_f)], rep(0,nrow(murp$rawdata))),
                       alpha = c(rep(0.8,nrow(tb_f)), rep(0.5,nrow(murp$rawdata))),
                       alpha_murp = c(rep(1,nrow(tb_f)), rep(0,nrow(murp$rawdata))),
                       size = c(tb_f$T, rep(0.01,nrow(murp$rawdata))),
                       size2 = c(tb_f$T, rep(0,nrow(murp$rawdata))),
                       sdata2)
  
  ### 计算每种细胞类型在每个MURP中所占的比例（中心点的比例设为0）
  if(!is.null(labels)){
    
    df2_list <- lapply(labels, function(lab){
      # cat(lab, "\n")
      tmp <- table(ct@MetaData[,c("murp_cluster",lab)]) %>% as.data.frame
      df <- reshape2::dcast(tmp, murp_cluster~get(lab))
      df <- apply(df, 2, as.numeric)
      rownames(df) <- df[,"murp_cluster"]
      df <- df[,-1,drop=FALSE]
      
      if(ncol(df)==1){
        df[,1] = rep(1, nrow(df))
      }else{
        df <- apply(df,1,function(x){ x/sum(x) })
        df <- t(df)
      }
      
      rownames(df) <- paste0("T[",rownames(df),"]")
      df <- df[rownames(tb_f),,drop=FALSE]
      
      tmp <- matrix(0, nr = nrow(object@assay$data_matrix), nc = ncol(df))
      rownames(tmp) <- rownames(murp$rawdata)
      df2 <- rbind(df, tmp)
    })
    df2 <- do.call(cbind, df2_list) %>% data.frame
    if(nrow(df2)==0) df2 = NULL
    
    ### 点的celltype属性
    if(is.null(labels)){
      meta_tmp = NULL
    }else if(length(labels)==1){
      meta_tmp <- c(sdf[rownames(tb_f), labels],
                    metadata[rownames(murp$rawdata),labels]) %>% data.frame
      colnames(meta_tmp) <- labels
    }else{
      meta_tmp <- rbind(sdf[rownames(tb_f), labels],
                        metadata[rownames(murp$rawdata),labels])
    }
    
    ### 合并
    vertex <- data.frame(vertex,meta_tmp,df2)
  }
  
  ## graph_tbtree_all
  graph_tbtree_all <- graph_from_data_frame(edge, directed = FALSE, vertices = vertex)
  
  ### 9. save
  s_ind <- unique(intersect(grep("T\\[",edge[,1]),grep("T\\[",edge[,2])))
  b_ind <- setdiff(1:nrow(edge), s_ind)
  edge$weight <- 10
  edge$weight[b_ind] <- 0.1
  edge$weight[s_ind] <- 1
  new_edge <- edge
  new_graph <- graph_from_data_frame(new_edge, directed = FALSE, vertices = vertex)
  start = rownames(sdf)[which.min(sdf$T)]
  fork_p = names(which(degree(tbtree$graph)>=3))
  
  ### 10. backbone
  bb <- layout_as_backbone(new_graph, keep = 0.4, backbone = TRUE)
  tbtree_all <- list(edge = new_edge,
                     vertex = vertex,
                     graph = new_graph,
                     adj = adj_final,
                     bb = bb,
                     start = start,
                     fork_p = fork_p)
  
  if(save){
    save(tbtree_all, file = "2_pseudotime/tbtree_all.rda")
  }
  object@Tree$tbtree_all = tbtree_all
  return(object)
}

#' GetTreeResult
#'
#' @param object
#'
GetBinTreeResult <- function(object){
  r = object@Tree[["bintree"]]
  return(r)
}

#' GetTbTreeResult
#'
#' @param object
#'
GetTbTreeResult <- function(object){
  r = object@Tree[["tbtree"]]
  return(r)
}

#' GetTbTreeAllResult
#'
#' @param object
#'
GetTbTreeAllResult <- function(object){
  r = object@Tree[["tbtree_all"]]
  return(r)
}

