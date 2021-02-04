#' Computes the robust SVD of a matrix
#' Using robRsvd
#'
#'
#' @param X Matrix. X matrix.
#' @param rank Integer. Rank of SVD decomposition
#'
#' @return List. The SVD of X.



get_svd_robustH <- function(X, rank=NULL){

  if(is.null(rank)){
    decomposition <- RobRSVD.all(X)
    decomposition
  } else{
    decomposition <- RobRSVD.all(X, nrank = rank)
    decomposition
  }

}


get_sv_threshold <- function(singular_values, rank){

  .5 * (singular_values[rank] + singular_values[rank + 1])
}


#' Truncates a robust SVD.
#'
#' Removes columns from the U, D, V matrix computed form an SVD.
#'
#'
#' @param decomposition List. List with entries 'u', 'd', and 'v'from the svd function.
#' @param rank List. List with entries 'u', 'd', and 'v'from the svd function.
#'
#' @return The trucated robust SVD of X.
truncate_svd <- function(decomposition, rank){

  if(rank==0){
    n <- dim(decomposition[['u']])[1]
    d <- dim(decomposition[['v']])[1]
    decomposition[['u']] <- matrix(0, ncol=1, nrow=n)
    decomposition[['d']] <- 0
    decomposition[['v']] <- matrix(0, ncol=1, nrow=d)
  }else{
    decomposition[['u']] <- decomposition[['u']][, 1:rank, drop=FALSE]
    decomposition[['d']] <- decomposition[['d']][1:rank]
    decomposition[['v']] <- decomposition[['v']][, 1:rank, drop=FALSE]
  }

  decomposition
}

#' Reconstruces the original matrix from its robust SVD.
#'
#' Computes UDV^T to get the approximate (or full) X matrix.
#'
#' @param decomposition List. List with entries 'u', 'd', and 'v'from the svd function.
#'
#' @return Matrix. The original matrix.
svd_reconstruction <- function(decomposition){

  # decomposition rank -- need to truncated singluar values
  r <- dim(decomposition[['u']])[2]

  decomposition[['u']]  %*%
    diag(decomposition[['d']][1:r], nrow=r, ncol=r) %*%
    t(decomposition[['v']])

}

#' Gets the wedin bounds
#'
#' @param X Matrix. The data matrix.
#' @param SVD List. The SVD decomposition of the matrix. List with entries 'u', 'd', and 'v'from the svd function.
#' @param signal_rank Integer.
#' @param num_samples Integer. Number of vectors selected for resampling procedure.


get_wedin_bound_samples <- function(X, SVD, signal_rank, num_samples=1000){

  # resample for U and V
  U_perp <- SVD[['u']][ , -(1:signal_rank)]
  U_sampled_norms <- wedin_bound_resampling(X=X,
                                            perp_basis=U_perp,
                                            right_vectors=FALSE,
                                            num_samples=num_samples)

  V_perp <- SVD[['v']][ , -(1:signal_rank)]
  V_sampled_norms <- wedin_bound_resampling(X=X,
                                            perp_basis=V_perp,
                                            right_vectors=TRUE,
                                            num_samples=num_samples)

  sigma_min <- SVD[['d']][signal_rank]
  wedin_bound_samples <- mapply(function(u, v)  min(max(u, v)/sigma_min, 1)^2, U_sampled_norms, V_sampled_norms)

  wedin_bound_samples
}

#' Resampling procedure for the wedin bound
#'
#' @param X Matrix. The data matrix.
#' @param perp_basis Matrix. Either U_perp or V_perp: the remaining left/right singluar vectors of X after estimating the signal rank.
#' @param right_vectors Boolean. Right multiplication or left multiplication.
#' @param num_samples Integer. Number of vectors selected for resampling procedure.
#' @importFrom foreach %dopar%

wedin_bound_resampling <- function(X, perp_basis, right_vectors, num_samples=1000){

  rank <- dim(perp_basis)[2]
  #resampled_norms <- rep(0, num_samples)
  numCores <- 2
  doParallel::registerDoParallel(numCores)
  resampled_norms <- foreach::foreach (s=1:num_samples) %dopar% {

    sampled_col_index <- sample.int(n=dim(perp_basis)[2],
                                    size=rank,
                                    replace=TRUE)


    perp_resampled <- perp_basis[ , sampled_col_index]

    if(right_vectors){
      resampled_projection <- X %*% perp_resampled
    } else{
      resampled_projection <- t(perp_resampled) %*% X
    }

    # operator L2 norm
    norm(resampled_projection,
         type='2')
  }

  as.numeric(resampled_norms)
}

#' Estimate the wedin bound for a data matrix.
#'
#' Samples from the random direction bound. Returns on the scale of squared singular value.
#'
#' @param n_obs The number of observations.
#' @param dims The number of features in each data matrix
#' @param num_samples Integer. Number of vectors selected for resampling procedure.
#' @importFrom stats rnorm
#'
#' @return rand_dir_samples

get_random_direction_bound_robustH <- function(n_obs, dims, num_samples=1000){

dims1 = as.list(dims)
n_blocks <- length(dims)
numCores <- 2
doParallel::registerDoParallel(numCores)

rand_dir_samples <- foreach::foreach (s=1:num_samples, .export=c("get_svd_robustH", "RobRSVD.all", "RobRSVD1")) %dopar% {

  X <- lapply(dims1, function(l) matrix(rnorm(n_obs * l, mean=0,sd=1), n_obs, l))
  rand_subspaces <- lapply(X, function(l) get_svd_robustH(l)[['u']])

  M <- do.call(cbind, rand_subspaces)
  M_svd <- get_svd_robustH(M, rank=min(dims))

  M_svd[['d']][1]^2

}

as.numeric(rand_dir_samples)
}





#' Block Scores
#'
#' Gets the block scores from the Rajive decomposition
#'
#' @param ajive_output List. The decomposition from Rajive
#' @param k Integer. The index of the data block
#' @param type Character. Joint or individual
#'
#' @return The block scores
#'
#' @examples
#' \donttest{
#'n <- 10
#'pks <- c(20, 10)
#'Y <- ajive.data.sim(K =2, rankJ = 2, rankA = c(7, 4), n = n,
#'                  pks = pks, dist.type = 1)
#'initial_signal_ranks <-  c(7, 4)
#'data.ajive <- list((Y$sim_data[[1]]), (Y$sim_data[[2]]))
#'ajive.results.robust <- Rajive(data.ajive, initial_signal_ranks)
#'get_block_scores(ajive.results.robust, 2, 'joint')
#'}
#' @export


get_block_scores <- function(ajive_output, k, type){
  if(! type  %in% c('joint', 'individual')){
    stop('type must be: joint or individual')
  }
  if (type == 'joint')   res <- ajive_output$block_decomps[[3*(k-1)+2]]$u
  if (type == 'individual')   res <- ajive_output$block_decomps[[3*(k-1)+1]]$u
  
  return(res)

  }


#' Block Loadings
#'
#' Gets the block loadings from the Rajive decomposition
#'
#' @param ajive_output List. The decomposition from Rajive
#' @param k Integer. The index of the data block
#' @param type Character. Joint or individual
#'
#' @return The block loadings
#' @examples 
#'\donttest{
#'n <- 10
#'pks <- c(20, 10)
#'Y <- ajive.data.sim(K =2, rankJ = 2, rankA = c(7, 4), n = n,
#'                  pks = pks, dist.type = 1)
#'initial_signal_ranks <-  c(7, 4)
#'data.ajive <- list((Y$sim_data[[1]]), (Y$sim_data[[2]]))
#'ajive.results.robust <- Rajive(data.ajive, initial_signal_ranks)
#'get_block_loadings(ajive.results.robust, 2, 'joint')
#'}

#'
#' @export
get_block_loadings <- function(ajive_output, k, type){
  if(! type  %in% c('joint', 'individual')){
    stop('type must be: joint or individual')
  }
 if (type == 'joint')   res <- ajive_output$block_decomps[[3*(k-1)+2]]$v
  if (type == 'individual')   res <- ajive_output$block_decomps[[3*(k-1)+1]]$v
return(res)
  }



#' Joint Rank
#'
#' Gets the joint rank from the Rajive decomposition
#'
#' @param ajive_output List. The decomposition from Rajive
#'
#' @return The joint rank
#' @examples  
#'\donttest{
#'n <- 10
#'pks <- c(20, 10)
#'Y <- ajive.data.sim(K =2, rankJ = 2, rankA = c(7, 4), n = n,
#'                  pks = pks, dist.type = 1)
#'initial_signal_ranks <-  c(7, 4)
#'data.ajive <- list((Y$sim_data[[1]]), (Y$sim_data[[2]]))
#'ajive.results.robust <- Rajive(data.ajive, initial_signal_ranks)
#'get_joint_rank(ajive.results.robust)
#'}

#' @export
get_joint_rank <- function(ajive_output){
  ajive_output$joint_rank
}

#' Individual Rank
#'
#' Gets the individual ranks from the Rajive decomposition
#'
#' @param ajive_output List. The decomposition from Rajive
#' @param k Integer. The index of the data block.
#'
#'
#' @return The individual ranks
#'
#'
#' @examples 
#' \donttest{
#'n <- 10
#'pks <- c(20, 10)
#'Y <- ajive.data.sim(K =2, rankJ = 2, rankA = c(7, 4), n = n,
#'                  pks = pks, dist.type = 1)
#'initial_signal_ranks <-  c(7, 4)
#'data.ajive <- list((Y$sim_data[[1]]), (Y$sim_data[[2]]))
#'ajive.results.robust <- Rajive(data.ajive, initial_signal_ranks)
#'get_individual_rank(ajive.results.robust, 2)
#'}

#' @export
get_individual_rank <- function(ajive_output, k){

  ajive_output$block_decomps[[3*(k-1)+1]]$rank
}


#' Decomposition Heatmaps
#'
#' Visualization of the RaJIVE decomposition, it shows heatmaps of the decomposition obtained by RaJIVE
#'
#'
#' @param blocks List. The initial data blocks.
#' @param jive_results_robust List. The RaJIVE decomposition.
#'
#' @return The heatmap of the decomposition
#'
#' @examples
#' \donttest{
#'n <- 10
#'pks <- c(20, 10)
#'Y <- ajive.data.sim(K =2, rankJ = 2, rankA = c(7, 4), n = n,
#'                  pks = pks, dist.type = 1)
#'initial_signal_ranks <-  c(7, 4)
#'data.ajive <- list((Y$sim_data[[1]]), (Y$sim_data[[2]]))
#'ajive.results.robust <- Rajive(data.ajive, initial_signal_ranks)
#'decomposition_heatmaps_robustH(data.ajive, ajive.results.robust)
#'}

#' @export




decomposition_heatmaps_robustH <- function (blocks, jive_results_robust)
{
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stop("The package 'cowplot' is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  K <- length(blocks)


heatmap_listR <- list()

heatmap_listR <- list()
for (k in 1:K) {
  heatmap_listR[[k]] <- data_heatmap(blocks[[k]], ylab = ifelse(k ==
                                                                  1, "observations", ""), show_color_bar = FALSE)
  heatmap_listR[[K + k]] <- data_heatmap(jive_results_robust$block_decomps[[3*(k-1)+2]][["full"]],
                                         ylab = ifelse(k == 1, "joint", ""), show_color_bar = FALSE)
  heatmap_listR[[2 * K + k]] <- data_heatmap(jive_results_robust$block_decomps[[3*(k-1)+1]][["full"]],
                                             ylab = ifelse(k == 1, "individual", ""),
                                             show_color_bar = FALSE)
  heatmap_listR[[3 * K + k]] <- data_heatmap(jive_results_robust$block_decomps[[3*k]],
                                             ylab = ifelse(k == 1, "noise", ""), show_color_bar = FALSE)
}
cowplot::plot_grid(plotlist = heatmap_listR, ncol = K)
}

#' Decomposition Heatmaps
#'
#' Visualization of the RaJIVE decomposition, it shows heatmaps of the decomposition obtained by RaJIVE
#'
#'
#' @param data List. The initial data blocks.
#' @param show_color_bar Boolean.
#' @param title Character.
#' @param xlab Character.
#' @param ylab Character
#'
#' @import ggplot2
#' @importFrom grDevices rainbow


data_heatmap <- function (data, show_color_bar = TRUE, title = "", xlab = "",
                          ylab = "")
{
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The package 'ggplot2' is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("reshape2", quietly = TRUE)) {
    stop("The package 'reshape2' is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  reshaped_data <- as.data.frame(reshape2::melt(data))
  colnames(reshaped_data) <- c("obs", "var", "value")
  ggplot(data = reshaped_data, aes_string(x = "var",
                                          y = "obs")) +
    geom_raster(aes_string(fill = "value"), show.legend = show_color_bar) +
    scale_fill_gradientn(colours = rainbow(10)) +
    theme(panel.background = element_blank(),  axis.line = element_blank(), legend.position = "bottom") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) + labs(title = title, x = xlab, y = ylab)
}



#' Proportions of variance explained
#'
#' Gets the variance explained by each component of the Rajive decomposition
#'
#' @param ajiveResults List. The decomposition from Rajive
#' @param blocks List. The initial data blocks
#'
#' @return The proportion of variance explained by each component
#' 
#' @examples 
#' \donttest{
#'n <- 10
#'pks <- c(20, 10)
#'Y <- ajive.data.sim(K =2, rankJ = 2, rankA = c(7, 4), n = n,
#'                  pks = pks, dist.type = 1)
#'initial_signal_ranks <-  c(7, 4)
#'data.ajive <- list((Y$sim_data[[1]]), (Y$sim_data[[2]]))
#'ajive.results.robust <- Rajive(data.ajive, initial_signal_ranks)
#'showVarExplained_robust(ajive.results.robust, data.ajive)
#'}
#'
#' @export

showVarExplained_robust <- function(ajiveResults, blocks){
  l <- length(blocks)
  # joint variance
  # joint is the second component for all 3
  VarJoint = rep(0, l)
  for (i in 1:l) VarJoint[i] = norm(as.matrix(ajiveResults$block_decomps[[3*(i-1)+2]][[1]]),
                                    type = "F")^2/norm(blocks[[i]], type = "F")^2

  # individual variances
  # individual is the first component for all 3
  VarIndiv = rep(0, l)
  for (i in 1:l) VarIndiv[i] = norm(as.matrix(ajiveResults$block_decomps[[3*(i-1)+1]][[1]]),
                                    type = "F")^2/norm(blocks[[i]], type = "F")^2

  # residual variance
  VarSubtr = 1 - VarJoint - VarIndiv

  VarProp <- list(VarJoint, VarIndiv, VarSubtr)
  names(VarProp) <- c('Joint', 'Indiv', 'Resid')
  VarProp
}
