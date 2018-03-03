#' @useDynLib wkmeans
NULL

#' Test
#'
#' @param x numeric matrix of data, or an object that can be coerced to such a matrix (such as a numeric vector or a data frame with all numeric columns).
#' @param k number of clusters.
#' @param w weights for each object (row) in x. By default, all observations are weighted the same, which yields regular kmeans.
#' @param iter_max maximum number of iterations allowed.
#' @param nstart number of random initialisation of clusters. Should be high for the clustering results to be stable.
#' @param cores numbers of cores to use in to try several initialisations in parallel (only useful if nstart > 1).
#' @return An object of class "\code{kmeans}", similar to that of function \code{\link[stats]{kmeans}}, which can use its \code{print} and \code{fitted} methods. The main difference is that the sum of squares are weighted. The list has the followinf components
#' \describe{
#' \item{cluster}{A vector of integers (from 1:k) indicating the cluster to which each point is allocated.}
#' \item{centers}{A matrix of cluster centres.}
#' \item{withindist}{Vector of (unweighted) distances to the centre of the assigned cluster (one element per input point).}
#' \item{totss}{The total weighted sum of squares.}
#' \item{withinss}{Vector of within-cluster weighted sum of squares, one component per cluster.}
#' \item{tot.withinss}{Total within-cluster weighted sum of squares, i.e. sum(withinss).}
#' \item{betweenss}{The between-cluster sum of squares, i.e. totss-tot.withinss.}
#' \item{size}{The number of points in each cluster.}
#' }
#' @examples
#' # start with fake data
#' x <- matrix(runif(1000*2), ncol=2)
#' # compute and plot kmeans clusters
#' g <- wkmeans(x, k=4, iter_max=100, nstart=100)
#' plot(x, col=factor(g$cluster), asp=1)
#' # colour points according to the squared distance from the cluster center
#' # (which is used to determine the clustering)
#' plot(x, col=heat.colors(10)[cut(g$withindist, breaks=10)], pch=16, asp=1)
#' points(x, col=factor(g$cluster))
#' points(g$centers, pch=16, col=1:4, cex=2)
#'
#' # put more weight on the left side (x < 0.5) and re-cluster
#' w <- rep(1, times=nrow(x))
#' w[x[,1]<0.5] <- 10
#' g <- wkmeans(x, k=4, iter_max=100, nstart=100, w=w)
#' plot(x, col=factor(g$cluster), asp=1)
#' # more clusters are created on the left, as expected
#' @export
#' @importFrom parallel mclapply
wkmeans <- function(x, k, w=rep(1, nrow(x)), iter_max=10, nstart=1, cores=1) {
  # check arguments
  if (length(w) != nrow(x)) {
    stop("Need as many weights as there are rows in x")
  }
  # convert input to matrix, for fortran
  x <- as.matrix(x, rownames.force=F)
  if (!is.numeric(x)) {
    stop("x is not all numeric or cannot be coerced to a numeric matrix")
  }
  # initialise sum of squares
  min_ss <- Inf
  # precompute dimensions for fortran
  n <- nrow(x)
  m <- ncol(x)

  # try several random starts
  res <- parallel::mclapply(1:nstart, function(i) {
    # pick k points as random starts
    centroids <- x[sample.int(n, k, replace=F),]
    # compute kmeans
    res <- .Fortran("wkmeans",
      # input variables
      X=x,
      k=as.integer(k),
      W=as.double(w),
      iter_max=as.integer(iter_max),
      # fortran call helpers
      n=n, m=m,            # dimensions
      centroids=centroids, # initial centroids
      clusters=rep(1L, n), # dummy cluster variable
      ss=0.,               # dummy sum of squares variable
      sqdist=rep(0., n)    # dummy squadred distance variable
    )
  }, mc.cores=min(cores, nstart))

  # keep only the best one
  sss <- sapply(res, `[[`, "ss")
  out <- res[[which.min(sss)]]

  # compute the various sum of squares
  # total sum of squares = (weighted) distance to total barycenter
  center <- apply(x, 2, mean)
  totss <- sum(sweep(x, 2, center)^2*w)
  # total (weighted) distance within each cluster
  withinss <- tapply(out$sqdist*w, out$clusters, sum)

  # format output like that of kmeans()
  res <- list(
    cluster=out$clusters,
    centers=out$centroids,
    withindist=sqrt(out$sqdist),
    totss=totss,
    withinss=withinss,
    tot.withinss=out$ss,
    betweenss=totss-out$ss,
    size=table(out$cluster)
  )
  class(res) <- c("kmeans", "list")
  return(res)
}
