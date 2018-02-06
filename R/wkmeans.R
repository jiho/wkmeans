#' @useDynLib wkmeans
NULL

#' Test
#'
#' @param x numeric matrix of data, or an object that can be coerced to such a matrix (such as a numeric vector or a data frame with all numeric columns).
#' @param k number of clusters.
#' @param w weights for each object (row) in x. By default, all observations are weighted the same, which yields regular kmeans.
#' @param iter_max maximum number of iterations allowed.
#' @param nstart number of random initialisation of clusters. Should be high for the clustering results to be stable.
#' @examples
#' # start with fake data
#' x <- matrix(runif(1000*2), ncol=2)
#' # compute and plot kmeans clusters
#' g <- wkmeans(x, k=4, iter_max=100, nstart=100)
#' plot(x, col=factor(g$cluster))
#' # put more weight on the left side (x < 0.5) and re-cluster
#' w <- rep(1, times=nrow(x))
#' w[x[,1]<0.5] <- 10
#' g <- wkmeans(x, k=4, iter_max=100, nstart=100, w=w)
#' plot(x, col=factor(g$cluster))
#' # more clusters are created on the left, as expected
#' @export
wkmeans <- function(x, k, w=rep(1, nrow(x)), iter_max=10, nstart=1) {
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
  for (i in 1:nstart) {
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
      ss=0.                # dummy sum of squares variable
    )

    # keep the result only if it is better than the previous
    if (res$ss < min_ss) {
      # message(i, ", ss=", res$ss)
      min_ss <- res$ss
      out <- res
    }
  }

  # format output like that of kmeans()
  return(list(
    cluster=out$clusters,
    centers=out$centroids,
    tot.withinss=out$ss,
    size=table(out$cluster)
  ))
}
