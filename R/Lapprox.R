#' @title Compute a rank `k` approximation of a graph Laplacian
#' @name Lapprox
#' @description This function calculates the rank-`k` approximation of a graph Laplacian (or any symmetric matrix). This function performs \link[base]{eigen} decomposition on the given matrix `L` and reconstructs it using only the LAST `k` eigenvectors and eigenvalues.
#' @param LapList A list of Laplacian matrices
#' @param k A vector indicating how many eigenvectors to take from each Laplacian, i.e., the number of clusters in each view
#' @param plots Whether or not to plot the eigenvalues from the rank approximated Laplacians
#'
#' @returns An n\eqn{\times}n matrix
#' @examples
#'
#' ## Generating data with 2 and 3 distinct clusters
#' ## Note that 'clustStruct' returns a list
#' n=120; k <- c(2,3)
#' set.seed(23)
#' dd <- clustStruct(n=n, p=30, k=k, noiseDat='random')
#'
#' ## Laplacians
#' L_list <- lapply(dd, kernelLaplacian, kernel="Spectrum", plots=FALSE, verbose=FALSE)
#'
#' trueGroups(n,k)
#'
#' La <- Lapprox(L_list, k=k, plots=FALSE)
#'
#' kmeans(La$vectors[,1:4], centers=4)
#' @export
Lapprox <- function(LapList, k, plots=TRUE){

  LrL <- mapply(rankL, LapList, k, SIMPLIFY = FALSE)
  Lr <- Reduce('+', LrL)
  eLr <- eigen(Lr)

  if(plots) plot(eLr$values, col='dodgerblue3', pch=20, type="b",
                 ylab="Eigen values", main = "Eigen values of the Pooled Laplacian")

  eLr
}


rankL <- function(L, k){
  eig <- eigen(L)

  vec <- eig$vector
  ind <- (ncol(vec)-k+1):ncol(vec)
  ev <- vec[,ind, drop=FALSE]
  val <- eig$values[ind]

  sweep(ev, 2, val, "*")%*%t(ev)
}