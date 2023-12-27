#!/usr/local/bin/Rscript

task <- dyncli::main()
# suppressMessages(library(dplyr, warn.conflicts = FALSE))
# suppressMessages(library(purrr, warn.conflicts = FALSE))
# suppressMessages(library(MURP, warn.conflicts = FALSE))
# suppressMessages(library(stringr, warn.conflicts = FALSE))
# suppressMessages(library(JuliaCall, warn.conflicts = FALSE))
# suppressMessages(library(pbmcapply, warn.conflicts = FALSE))
# suppressMessages(library(doParallel, warn.conflicts = FALSE))
# suppressMessages(library(reshape, warn.conflicts = FALSE))
# suppressMessages(library(reshape2, warn.conflicts = FALSE))
# suppressMessages(library(igraph, warn.conflicts = FALSE))
# suppressMessages(library(graphlayouts, warn.conflicts = FALSE))

library(dplyr, warn.conflicts = FALSE)
library(purrr, warn.conflicts = FALSE)
library(MURP, warn.conflicts = FALSE)
library(stringr, warn.conflicts = FALSE)
library(JuliaCall, warn.conflicts = FALSE)
library(pbmcapply, warn.conflicts = FALSE)
library(doParallel, warn.conflicts = FALSE)
library(reshape, warn.conflicts = FALSE)
library(reshape2, warn.conflicts = FALSE)
library(igraph, warn.conflicts = FALSE)
library(graphlayouts, warn.conflicts = FALSE)
invisible()

# source("/share/data6/tmp/renjun/CellTrekResult/CellTrek/ti_mgpfact2/function.R")
# julia_home = "/public/home/renjun/tool/julia-1.6.6/bin"

source("/code/function.R")
julia_home = "/usr/local/bin"

counts <- as.matrix(task$counts)
parameters <- task$parameters
start_id <- task$priors$start_id
cat(start_id,"\n")
tag <- parameters$data_name
path <- parameters$save_path

## debug
# save(parameters, file = "~/parameters.rda")
# tag <- parameters$data_name
# cat("::::tag: ", tag, "\n")
# save(task, file = "~/task.rda")
# path <- parameters$save_path
# path = getwd()
# cat("::::path: ", path, "\n")

############################ 测试专用
# wrap = dataset
# expression <- wrap$expression
# counts <- wrap$counts
# parameters = list(dataset_id = "tree_10",
#                   save_path = getwd(),
#                   omega = 0.167,
#                   max_murp = 5,
#                   trajectory_number = 3,
#                   murp_pc_number = 3,
#                   murp_pca_center = as.logical(TRUE),
#                   murp_pca_scale = as.logical(TRUE),
#                   iterations = 5,
#                   chains_number = 1,
#                   trajectory_type = "consensus_tree")
# start_id  = wrap$prior_information$start_id
############################

##############################################################################
#
#                                  create
#
##############################################################################

checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

cat("\n\n normalize data \n\n")
expression = LogNormalize(t(counts)) %>% t

cat("\n\n create MGPfact object \n\n")
# hvg_num <- round(ncol(counts) * 0.1)
# seurat <- FindVariableFeatures(object = seurat, selection.method = "vst", nfeatures = hvg_num)
# expression <- seurat@assays$RNA@data[VariableFeatures(seurat),] %>% as.matrix %>% t
# expression <- seurat@assays$RNA@data %>% as.matrix %>% t
cell_info = data.frame(row.names = rownames(expression), cell_id = rownames(expression))
ct <- CreateMGPfactObject(data_matrix = as.matrix(expression), MetaData = cell_info, dir = ".")

##############################################################################
#
#                                   MURP
#
##############################################################################

cat("\n\n MURP downsampling \n\n")
if(parameters$max_murp > nrow(expression)){
  max_murp = nrow(expression)
}else{
  max_murp = parameters$max_murp
}

sink(file=NULL) 
ct =  MURPDownsampling(ct, omega = parameters$omega, max_murp = max_murp, 
                   iter = 10, seed = 723, fast = T, cores = 1,
                   pca.center = parameters$murp_pca_center, 
                   pca.scale = parameters$murp_pca_scale,
                   plot = F)
sink()

##############################################################################
#
#                            RunningmodMGPpseudoT
#
##############################################################################

cat("\n\n running pseudotime model \n\n")

# need_root = parameters$need_root
need_root = TRUE
if(need_root){
  start_id <- if (!is.null(start_id)) { sample(start_id, 1) } else { NULL }
  if(is.null(start_id)){
    murp_start_id = 999
  }else{
    ind = which(names(ct@MURP$Recommended_K_cl$cluster)==start_id)
    murp_start_id = as.vector(ct@MURP$Recommended_K_cl$cluster[ind])
  }
}else{
  murp_start_id = 999
}
# cat("roottttt:: ", murp_start_id, "\n")

murp_pca = ct@MURP$centersPCA$x
traj_number = ifelse(ncol(murp_pca)>parameters$trajectory_number,
                     parameters$trajectory_number,
                     ncol(murp_pca))

if(ncol(ct@MURP$centersPCA$x)<parameters$murp_pc_number){
  murp_pc_number = ncol(ct@MURP$centersPCA$x)
}else{
  murp_pc_number = parameters$murp_pc_number
}

SaveMURPDatToJulia(ct, murp_pc_number = murp_pc_number)
ct = SetSettings(ct, 
                 murp_pc_number = murp_pc_number, 
                 trajectory_number = traj_number, 
                 pse_optim_iterations = parameters$iterations, 
                 start_murp = murp_start_id,
                 chains_number = parameters$chains_number)
ct = RunningmodMGPpseudoT(ct, julia_home = julia_home, seed = 723, cores = 1)

##############################################################################
#
#                                GetPseSDF
#
##############################################################################

ct <- GetIterCor(ct, iteration_list = list(c(1, getParams(ct,"pse_optim_iterations"))))
ct <- GetPredT(object = ct, chains = 1:getParams(ct,"chains_number"), adjust = TRUE, filter_chain = TRUE, mean_th = 0.3)
ct <- GetPseSdf(ct, unified_direction = FALSE, rm_adjust_chain = TRUE)
sdf <- GetMURPInfo(ct)

##############################################################################
#
#                                 tbtree
#
##############################################################################

cat("\n\n construct tree \n\n")

ct <- GetBinTree(object = ct)
ct <- GetTbTree(object = ct)
ct <- GetTbTreeAllpoint(object = ct, save = F, labels = getParams(ct,"label"))

bintree_all <- GetBinTreeResult(ct)
tbtree <- GetTbTreeResult(ct)
tbtree_all <- GetTbTreeAllResult(ct)

orig_T = ct@MetaData$T %>% magrittr::set_names(rownames(ct@MetaData))
checkpoints$method_aftermethod <- as.numeric(Sys.time())

##############################################################################
#
#                         process MGPfact output
#
##############################################################################

vertex = tbtree_all$vertex
dima =  tbtree_all$bb$xy %>% data.frame
rownames(dima) = rownames(vertex)

dimred = dima[(ct@MURP$Recommended_K+1):nrow(dima),] %>%
  magrittr::set_colnames(c("comp_1", "comp_2")) %>%
  as.matrix
dimred_milestones = dima[1:ct@MURP$Recommended_K,]   %>%
  magrittr::set_colnames(c("comp_1", "comp_2")) %>%
  as.matrix
milestone_network <-
  igraph::as_data_frame(tbtree$graph) %>%
  transmute(
    from, to,
    length = sqrt(rowSums((dimred_milestones[from,,drop=F] - dimred_milestones[to,,drop=F])^2)),
    # length = sdf[to,"T"]-sdf[from,"T"],
    directed = TRUE
  )
milestone_percentages = data.frame(cell_id = rownames(ct@MetaData),
                                   milestone_id = paste0("T[", ct@MetaData$murp_cluster, "]"),
                                   percentage = 1)

dimred_segment_progressions <-
  milestone_network %>%
  dplyr::select(from, to) %>%
  mutate(percentage = map(seq_len(n()), ~ c(0, 1))) %>%
  tidyr::unnest(percentage)

dsp_names <-
  dimred_segment_progressions %>%
  {ifelse(.$percentage == 0, .$from, .$to)}
dimred_segment_points <- dimred_milestones[dsp_names, , drop = FALSE]

##############################################################################
#
#                              process output
#
##############################################################################

## linear
output_linear <-
  dynwrap::wrap_data( cell_ids = rownames(expression) ) %>%
  # dynwrap::add_dimred( dimred = as.matrix(dimred) ) %>%
  # dynwrap::add_linear_trajectory(pseudotime = set_names(orig_T, rownames(expression))) %>%
  dynwrap::add_linear_trajectory(pseudotime = orig_T) %>%
  dynwrap::add_timings( timings = checkpoints )

## consensus tree
output_consensus_tree <-
  dynwrap::wrap_data(
    cell_ids = rownames(expression)
  ) %>%
  # dynwrap::add_pseudotime(pseudotime = orig_T) %>%
  dynwrap::add_trajectory(
    milestone_ids = rownames(dimred_milestones),
    milestone_network = milestone_network,
    milestone_percentages = milestone_percentages
  ) %>%
  dynwrap::add_timings( timings = checkpoints )

## consensus tree projection
output_consensus_tree_proj <-
  dynwrap::wrap_data( cell_ids = rownames(expression) ) %>%
  dynwrap::add_dimred_projection(
    milestone_ids = rownames(dimred_milestones),
    milestone_network = milestone_network,
    dimred = dimred,
    dimred_milestones = dimred_milestones
  ) %>%
  dynwrap::add_timings( timings = checkpoints )

## binary tree 
output_binarytree_list <- lapply(seq_along(bintree_all), function(l){
  milestone_network <-
    igraph::as_data_frame(bintree_all[[l]]$graph) %>%
    transmute( from, to, 
               length = sqrt(rowSums((dimred_milestones[from, ] - dimred_milestones[to, ])^2)),
               # length = 1, 
               directed = TRUE )
  milestone_percentages = data.frame(cell_id = rownames(ct@MetaData),
                                     milestone_id = paste0("T[", ct@MetaData$murp_cluster, "]"),
                                     percentage = 1)
  output <- 
    dynwrap::wrap_data( cell_ids = rownames(expression) ) %>%
    dynwrap::add_trajectory(
      milestone_ids = rownames(dimred_milestones),
      milestone_network = milestone_network,
      milestone_percentages = milestone_percentages
    ) %>%
    dynwrap::add_timings(timings = checkpoints)
  return(output)
})

##############################################################################
#
#                                save result
#
##############################################################################

cat("\n\n save result \n\n")
# dataset_id = parameters$dataset_id
# dataset_id = gsub("\\/","_",dataset_id)
# usedr_name = paste0(dataset_id, 
#                   "_om_", parameters$omega, 
#                   "_maxmm_", parameters$max_murp, 
#                   "_trajn_", parameters$trajectory_number, 
#                   "_mpc_", parameters$murp_pc_number, 
#                   "_mpcenter_", parameters$murp_pca_center,
#                   "_mpscale_", parameters$murp_pca_scale,
#                   "_iter_", parameters$iterations, 
#                   "_ch_", parameters$chains_number)
# dir.create(usedr_base, recursive = T)
save(orig_T, ct, checkpoints, 
     output_linear, output_consensus_tree, output_consensus_tree_proj, output_binarytree_list,
     file = paste0(parameters$save_path,"/other_output.rda") )

cat("\n\n write output \n\n")
if(parameters$trajectory_type=="consensus_tree"){ output = output_consensus_tree}
if(parameters$trajectory_type=="consensus_tree_proj"){ output = output_consensus_tree_proj}
if(parameters$trajectory_type=="linear"){ output = output_linear }
if(parameters$trajectory_type=="binary_tree"){ output = output_binarytree_list[[1]]}

dyncli::write_output(output, task$output)
