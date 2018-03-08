#' Traffic accidents (from the BAAC dataset)
#'
#' Traffic accidents that occurred in 2008 in two French départements - Finistère and Morbihan - and involved bodily injuries.
#' This is an extract from the BAAC dataset (Bulletins d'analyse d'accident corporel).
#' The data is filled by police forces, and most accidents have a specific location.
#' @docType data
#'
#' @usage data(acci)
#'
#' @format An list containing two elements : one for Finistère and the other for Morbihan.
#'     Each element contains a list of two elements : `points` and `polygon` which gives
#'     the coordinates (long,lat) of the accidents inside a `data.frame`, and a `data.frame` containing the
#'     bounding surface of the département (long, lat).
#'
#' @keywords datasets, accidents
#'
#' @references Charpentier, A. & Gallic, E. (2015). Kernel density estimation based on Ripley’s correction. GeoInformatica, 1-22.
#' (\href{https://link.springer.com/article/10.1007/s10707-015-0232-z}{Springer})
#'
#'
#' @examples
#' data(acci)
#' pts_finistere <- acci$finistere$points
#' pol_finistere <- acci$finistere$polygon
#' plot(pol_finistere, t="l")
#' points(pts_finistere, col = "red", pch=19)
"acci"


