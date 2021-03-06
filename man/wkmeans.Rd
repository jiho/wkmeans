% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wkmeans.R
\name{wkmeans}
\alias{wkmeans}
\title{Test}
\usage{
wkmeans(x, k, w = rep(1, nrow(x)), iter_max = 10, nstart = 1, cores = 1)
}
\arguments{
\item{x}{numeric matrix of data, or an object that can be coerced to such a matrix (such as a numeric vector or a data frame with all numeric columns).}

\item{k}{number of clusters.}

\item{w}{weights for each object (row) in x. By default, all observations are weighted the same, which yields regular kmeans.}

\item{iter_max}{maximum number of iterations allowed.}

\item{nstart}{number of random initialisation of clusters. Should be high for the clustering results to be stable.}

\item{cores}{numbers of cores to use in to try several initialisations in parallel (only useful if nstart > 1).}
}
\value{
An object of class "\code{kmeans}", similar to that of function \code{\link[stats]{kmeans}}, which can use its \code{print} and \code{fitted} methods. The main difference is that the sum of squares are weighted. The list has the followinf components
\describe{
\item{cluster}{A vector of integers (from 1:k) indicating the cluster to which each point is allocated.}
\item{centers}{A matrix of cluster centres.}
\item{withindist}{Vector of (unweighted) distances to the centre of the assigned cluster (one element per input point).}
\item{totss}{The total weighted sum of squares.}
\item{withinss}{Vector of within-cluster weighted sum of squares, one component per cluster.}
\item{tot.withinss}{Total within-cluster weighted sum of squares, i.e. sum(withinss).}
\item{betweenss}{The between-cluster sum of squares, i.e. totss-tot.withinss.}
\item{size}{The number of points in each cluster.}
}
}
\description{
Test
}
\examples{
# start with fake data
x <- matrix(runif(1000*2), ncol=2)
# compute and plot kmeans clusters
g <- wkmeans(x, k=4, iter_max=100, nstart=100)
plot(x, col=factor(g$cluster), asp=1)
# colour points according to the squared distance from the cluster center
# (which is used to determine the clustering)
plot(x, col=heat.colors(10)[cut(g$withindist, breaks=10)], pch=16, asp=1)
points(x, col=factor(g$cluster))
points(g$centers, pch=16, col=1:4, cex=2)

# put more weight on the left side (x < 0.5) and re-cluster
w <- rep(1, times=nrow(x))
w[x[,1]<0.5] <- 10
g <- wkmeans(x, k=4, iter_max=100, nstart=100, w=w)
plot(x, col=factor(g$cluster), asp=1)
# more clusters are created on the left, as expected
}
