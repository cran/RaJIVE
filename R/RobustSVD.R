

#' Computes the robust SVD of a matrix
#'
#'
#' @param data Matrix. X matrix.
#' @param nrank Integer. Rank of SVD decomposition
#' @param svdinit List. The standard SVD.
#' @importFrom stats median
#' @return List. The SVD of X.

RobRSVD.all <- function(data, nrank = min(dim(data)), svdinit = svd(data))
{
  data.svd1<-RobRSVD1(data, sinit = svdinit$d[1],
                      uinit = svdinit$u[,1], vinit = svdinit$v[,1])
  d = data.svd1$s
  u = data.svd1$u
  v = data.svd1$v
  Red = d*u%*%t(v)
  Rm = min(min(dim(data)), nrank)
  for(i in 1: (Rm-1)){
    data.svd1 <- RobRSVD1((data-Red), sinit = svdinit$d[1],
                        uinit = svdinit$u[,1],vinit = svdinit$v[,1])
    d = c(d,data.svd1$s)
    u = cbind(u,data.svd1$u)
    v = cbind(v,data.svd1$v)
    Red <- (u%*%diag(d)%*%t(v))
  }
  out = list(d,u,v)
  names(out) <- c("d", "u", "v")
  return(out)

}


RobRSVD1<- function (data, huberk = 1.345, niter = 1000,
                     tol = 1e-05, sinit, uinit, vinit)
{
  size_data = c(dim(data))
  m = size_data[1]
  n = size_data[2]
  sold = sinit
  vold = vinit

  uold = sold * uinit
  Appold = uold %*% t(vold)
  Rmat = data - Appold
  Rvec = c(Rmat)
  mysigma = median(abs(Rvec))/0.675
  iter = 1
  localdiff = 9999
  while (localdiff > tol & iter < niter) {
    Wmat = huberk/abs(Rmat/mysigma)
    Wmat[Wmat > 1] = 1


    uterm1 = diag(colSums(diag(c(vold^2)) %*% t(Wmat))) +
      (2 * mysigma^2) * (c(t(vold) %*% (diag(n)) %*% vold) * (diag(m)) - diag(sum(vold^2), m))
    uterm2 = (Wmat * data) %*% vold
    unew = solve(uterm1) %*% uterm2

    vterm1 = diag(colSums(diag(c(unew^2)) %*% Wmat)) +
      (2 * mysigma^2) * (c(t(unew) %*% (diag(m)) %*% unew) * (diag(n)) - diag(sum(unew^2), n))
    vterm2 = t(Wmat * data) %*% unew
    vnew = solve(vterm1) %*% vterm2

    Appnew = unew %*% t(vnew)
    Rmat = data - Appnew
    localdiff = max(abs(Appnew - Appold))
    Appold = Appnew
    uold = sqrt(sum(vnew^2)) * unew
    vold = vnew/sqrt(sum(vnew^2))
    iter = iter + 1
  }
  v = vold
  s = sqrt(sum(uold^2))
  u = uold/sqrt(sum(uold^2))
  return(list(s = s, v = v, u = u))
}
