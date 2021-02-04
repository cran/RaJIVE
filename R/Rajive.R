#' Robust Angle based Joint and Individual Variation Explained
#'
#' Computes the robust aJIVE decomposition with parallel computation.
#'
#' @param blocks List. A list of the data matrices.
#' @param initial_signal_ranks Vector. The initial signal rank estimates.
#' @param full Boolean. Whether or not to store the full J, I, E matrices or just their SVDs (set to FALSE to save memory).
#' @param n_wedin_samples Integer. Number of wedin bound samples to draw for each data matrix.
#' @param n_rand_dir_samples Integer. Number of random direction bound samples to draw.
#' @param joint_rank Integer or NA. User specified joint_rank. If NA will be estimated from data.

#' @return The aJIVE decomposition.
#'
#' @examples
#' \donttest{
#' n <- 50
#' pks <- c(100, 80, 50)
#' Y <- ajive.data.sim(K =3, rankJ = 3, rankA = c(7, 6, 4), n = n,
#'                    pks = pks, dist.type = 1)
#' initial_signal_ranks <-  c(7, 6, 4)
#' data.ajive <- list((Y$sim_data[[1]]), (Y$sim_data[[2]]), (Y$sim_data[[3]]))
#' ajive.results.robust <- Rajive(data.ajive, initial_signal_ranks)
#' }
#'
#' @export


Rajive <- function(blocks, initial_signal_ranks, full=TRUE,
                           n_wedin_samples=1000, n_rand_dir_samples=1000, joint_rank=NA)
{

  K <- length(blocks)


  # step 1: initial signal space extraction --------------------------------
  # initial estimate of signal space with SVD

  # apply svd to all element in the list
  block_svd <- parallel::mclapply(blocks, get_svd_robustH)
  # extract singular values from list
  singular_values = parallel::mclapply(block_svd, function(l) l[[1]])
  # apply get_sv_threshold
  sv_thresholds <- mapply(function(l,p)
    get_sv_threshold(l, rank = p), singular_values, initial_signal_ranks)




  # step 2: joint sapce estimation -------------------------------------------------------------

  out <- get_joint_scores_robustH(blocks, block_svd, initial_signal_ranks, sv_thresholds,
                                  n_wedin_samples=n_wedin_samples,
                                  n_rand_dir_samples=n_rand_dir_samples,
                                  joint_rank=joint_rank)
  joint_rank_sel_results <- out$rank_sel_results
  joint_scores <- out$joint_scores

  joint_rank <- dim(joint_scores)[2]


  # step 3: final decomposition -----------------------------------------------------

  block_decomps <- mapply(function(l,m)
    get_final_decomposition_robustH(l, joint_scores = joint_scores, m), blocks,
    sv_thresholds )



  jive_decomposition <- list(block_decomps=block_decomps)
  jive_decomposition[['joint_scores']] <- joint_scores
  jive_decomposition[['joint_rank']] <- joint_rank

  jive_decomposition[['joint_rank_sel']] <- joint_rank_sel_results
  jive_decomposition
}


#' The singular value threshold.
#'
#' Computes the singular value threshold for the data matrix (half way between the rank and rank + 1 singluar value).
#'
#' @param singular_values Numeric. The singular values.
#' @param rank Integer. The rank of the approximation.
#'

get_sv_threshold <- function(singular_values, rank){

.5 * (singular_values[rank] + singular_values[rank + 1])
}


#' Computes the joint scores.
#'
#' Estimate the joint rank with the wedin bound, compute the signal scores SVD, double check each joint component.
#'
#' @param blocks List. A list of the data matrices.
#' @param block_svd List. The SVD of the data blocks.
#' @param initial_signal_ranks Numeric vector. Initial signal ranks estimates.
#' @param sv_thresholds Numeric vector. The singular value thresholds from the initial signal rank estimates.
#' @param n_wedin_samples Integer. Number of wedin bound samples to draw for each data matrix.
#' @param n_rand_dir_samples Integer. Number of random direction bound samples to draw.
#' @param joint_rank Integer or NA. User specified joint_rank. If NA will be estimated from data.
#' @importFrom stats quantile
#'
#'
#'

get_joint_scores_robustH <- function(blocks, block_svd, initial_signal_ranks, sv_thresholds,
                                     n_wedin_samples=1000, n_rand_dir_samples=1000,
                                     joint_rank=NA){


  if(is.na(n_wedin_samples) & is.na(n_rand_dir_samples) & is.na(joint_rank)){
    stop('at least one of n_wedin_samples, n_rand_dir_samples, or joint_rank must not be NA',
         call.=FALSE)
  }

  K <- length(blocks)
  n_obs <- dim(blocks[[1]])[1]

  # SVD of the signal scores matrix -----------------------------------------
  signal_scores <- list()
  for(k in 1:K){
    signal_scores[[k]] <- block_svd[[k]][['u']][, 1:initial_signal_ranks[k]]
  }

  M <- do.call(cbind, signal_scores)
  M_svd <- get_svd_robustH(M, rank=min(initial_signal_ranks))


  # estimate joint rank with wedin bound and random direction bound -------------------------------------------------------------

  rank_sel_results  <- list()
  rank_sel_results[['obs_svals']] <- M_svd[['d']]

  if(is.na(joint_rank)){

    # maybe comptue wedin bound
    if(!is.na(n_wedin_samples)){

      block_wedin_samples <- t(mapply(function(l,m) get_wedin_bound_samples(l,
                                                                            m,
                                                                            signal_rank=initial_signal_ranks[k],
                                                                            num_samples=n_wedin_samples),
                                      blocks, block_svd))




      wedin_samples <-  K - colSums(block_wedin_samples)
      wedin_svsq_threshold <- quantile(wedin_samples, .05)

      rank_sel_results[['wedin']] <- list(block_wedin_samples=block_wedin_samples,
                                          wedin_samples=wedin_samples,
                                          wedin_svsq_threshold=wedin_svsq_threshold)
    } else{
      wedin_svsq_threshold <- NA
    }

    # maybe compute random direction bound
    if(!is.na(n_rand_dir_samples)){

      rand_dir_samples <- get_random_direction_bound_robustH(n_obs=n_obs, dims=initial_signal_ranks,
                                                             num_samples=n_rand_dir_samples)
      rand_dir_svsq_threshold <- quantile(rand_dir_samples, .95)

      rank_sel_results[['rand_dir']] <- list(rand_dir_samples=rand_dir_samples,
                                             rand_dir_svsq_threshold=rand_dir_svsq_threshold)

    } else {
      rand_dir_svsq_threshold <- NA
    }

    overall_sv_sq_threshold <- max(wedin_svsq_threshold, rand_dir_svsq_threshold, na.rm=TRUE)
    joint_rank_estimate <- sum(M_svd[['d']]^2 > overall_sv_sq_threshold)

    rank_sel_results[['overall_sv_sq_threshold']] <- overall_sv_sq_threshold
    rank_sel_results[['joint_rank_estimate']] <- joint_rank_estimate


  } else { # user provided joint rank
    joint_rank_estimate <- joint_rank
    rank_sel_results[['joint_rank_estimate']] <- joint_rank
  }


  # estimate joint score space ------------------------------------


  joint_scores <- M_svd[['u']][ , 1:joint_rank_estimate, drop=FALSE]

  # reconsider joint score space ------------------------------------
  # remove columns of joint_scores that have a
  # trivial projection from one of the data matrices

  to_remove <- c()
  for(k in 1:K){
    for(j in 1:joint_rank_estimate){

      score <- t(blocks[[k]]) %*% joint_scores[ , j]
      sv <- norm(score)

      if(sv < sv_thresholds[[k]]){
        print(paste('removing column', j))
        to_remove <- c(to_remove, j)
        break
      }
    }

  }
  to_keep <- setdiff(1:joint_rank_estimate, to_remove)
  joint_rank <- length(to_keep)
  joint_scores <- joint_scores[ , to_keep, drop=FALSE]

  list(joint_scores=joint_scores, rank_sel_results=rank_sel_results)
}


#' Computes the final JIVE decomposition.
#'
#' Computes X = J + I + E for a single data block and the respective SVDs.
#'
#'
#' @param X Matrix. The original data matrix.
#' @param joint_scores Matrix. The basis of the joint space (dimension n x joint_rank).
#' @param sv_threshold Numeric vector. The singular value thresholds from the initial signal rank estimates.
#' @param full Boolean. Do we compute the full J, I matrices or just svd
#'
get_final_decomposition_robustH <- function(X, joint_scores, sv_threshold, full=TRUE){

  jive_decomposition <- list()
  jive_decomposition[['individual']] <- get_individual_decomposition_robustH(X, joint_scores, sv_threshold, full)
  jive_decomposition[['joint']] <- get_joint_decomposition_robustH(X, joint_scores, full)


  if(full){
    jive_decomposition[['noise']] <- X - (jive_decomposition[['joint']][['full']] +
                                            jive_decomposition[['individual']][['full']])
  } else{
    jive_decomposition[['noise']] <- NA
  }

  jive_decomposition
}


#' Computes the individual matrix for a data block.
#'
#' @param X Matrix. The original data matrix.
#' @param joint_scores Matrix. The basis of the joint space (dimension n x joint_rank).
#' @param sv_threshold Numeric vector. The singular value thresholds from the initial signal rank estimates.
#' @param full Boolean. Do we compute the full J, I matrices or just the SVD (set to FALSE to save memory).

get_individual_decomposition_robustH <- function(X, joint_scores, sv_threshold, full=TRUE){

  X_orthog <- (diag(dim(X)[1]) - joint_scores %*% t(joint_scores)) %*% X

  indiv_decomposition <- get_svd_robustH(X_orthog)

  indiv_rank <- sum(indiv_decomposition[['d']] > sv_threshold)

  indiv_decomposition <- truncate_svd(decomposition=indiv_decomposition,
                                      rank=indiv_rank)
  if(full){
    indiv_decomposition[['full']] <- svd_reconstruction(indiv_decomposition)
  } else{
    indiv_decomposition[['full']] <- NA
  }

  indiv_decomposition[['rank']] <- indiv_rank
  indiv_decomposition
}


#' Computes the individual matrix for a data block
#'
#' @param X Matrix. The original data matrix.
#' @param joint_scores Matrix. The basis of the joint space (dimension n x joint_rank).
#' @param full Boolean. Do we compute the full J, I matrices or just the SVD (set to FALSE to save memory).

get_joint_decomposition_robustH <- function(X, joint_scores, full=TRUE){

  joint_rank <- dim(joint_scores)[2]
  J <-  joint_scores %*% t(joint_scores) %*% X

  joint_decomposition <- get_svd_robustH(J, joint_rank)

  if(full){
    joint_decomposition[['full']] <- J
  } else{
    joint_decomposition[['full']] <- NA
  }

  joint_decomposition[['rank']] <- joint_rank
  joint_decomposition

}



