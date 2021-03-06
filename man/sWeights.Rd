% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kdeborder.R
\name{sWeights}
\alias{sWeights}
\title{Proportion of the area of a circle on a polygon's area}
\usage{
sWeights(x, h, polygon)
}
\arguments{
\item{x}{center of the circle}

\item{h}{bandwidth scalar}

\item{polygon}{polygon on which data points lie}
}
\value{
Returns the proportion of the area of a circle of center x and radius 1.759*h on the polygon's area
}
\description{
Proportion of the area of a circle on a polygon's area
}
\examples{
pol_coordinates <- matrix(c(0, 0, 1, 1, 0, 0, 1, 1, 0, 0), ncol = 2)
pol <- as(pol_coordinates, "gpc.poly")
sWeights(x = c(0, 0), h = 1/1.759, polygon = pol)
}
\seealso{
\code{\link{sCircle}} which this function wraps
}
