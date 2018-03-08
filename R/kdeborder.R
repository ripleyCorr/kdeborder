#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

#' Creates "n" points on a circle centered in "centre" with a radius of "radius"
#'
#' @param n number of points to define the circle
#' @param centre center of the circle
#' @param radius radius of the circle
#' @return Returns a matrix with two columns corresponding to the (x,y) coordinates of the circle of centre `centre`
#'   and radius `radius`. The number of lines depends on the number of points asked using the parameter `n`.
#' @examples sCircle(n=100, centre = c(0,0), radius = 1)
sCircle <- function(n = 100, centre = c(0, 0), radius){
  theta <- seq(0, 2*pi, length = n)
  m <- cbind(cos(theta), sin(theta)) * radius
  m[, 1] <- m[, 1] + centre[1]
  m[, 2] <- m[, 2] + centre[2]
  colnames(m) <- c("x", "y")
  m
}# End of sCircle()

#' Proportion of the area of a circle on a polygon's area
#'
#' @param x center of the circle
#' @param h bandwidth scalar
#' @param polygon polygon on which data points lie
#' @seealso \code{\link{sCircle}} which this function wraps
#' @return Returns the proportion of the area of a circle of center x and radius 1.759*h on the polygon's area
#' @examples
#' pol_coordinates <- matrix(c(0, 0, 1, 1, 0, 0, 1, 1, 0, 0), ncol = 2)
#' pol <- as(pol_coordinates, "gpc.poly")
#' sWeights(x = c(0, 0), h = 1/1.759, polygon = pol)
sWeights <- function(x, h, polygon) {
  leCercle <- sCircle(centre = x, radius = 1.759*h)
  POLcercle <- as(leCercle[-nrow(leCercle),], "gpc.poly")
  area.poly(intersect(polygon, POLcercle)) / area.poly(POLcercle)
}# End of sWeights()


#' Computes an estimation of the density using Kernel Density estimator,
#' correcting for fontier effects.
#'
#' @return Returns a list whose elements are:
#'   > X:    x coordinates at wich estimate is evaluated,
#'   > Y:    y coordinates at wich estimate is evaluated,
#'   > Z:    density estimates,
#'   > ZNA:  density estimates with NA values for points outside the polygon,
#'   > H:    bandwidth matrix,
#'   > W:    vector of weights.
#' @seealso \code{\link{sWeights}} which this function wraps
#' @param U data points.
#' @param polygon polygon on which points lie.
#' @param optimal if TRUE, uses Hpi() to select the optimal bandwidth.
#' @param h only if optimal=FALSE, scalar bandwidth.
#' @param parallel if TRUE, computes the weights using clusters.
#' @param n_clusters only if n_clusters=TRUE, defines the number of clusters. (Set to NULL for automatic detection (on Unix)).
#' @references Charpentier, A. & Gallic, E. (2015). Kernel density estimation based on Ripley’s correction. GeoInformatica, 1-22.
#' (\href{https://link.springer.com/article/10.1007/s10707-015-0232-z}{Springer})
#' @examples
#' data(acci)
#' # Estimation with correction
#' smoothed_fin <- sKDE(U = acci$finistere$points, polygon = acci$finistere$polygon,
#' optimal=TRUE, parallel = FALSE)
sKDE <- function(U, polygon, optimal = TRUE, h = .1, parallel = FALSE, n_clusters = 4){
  if(!class(polygon) == "gpc.poly") polygon <- as(polygon, "gpc.poly")
  if(class(U) == "data.frame") U <- as.matrix(U)
  IND <- which(is.na(U[, 1]) == FALSE)
  U <- U[IND,]
  n <- nrow(U)
  if(optimal){
    H <- Hpi(U, binned = FALSE)
    H <- matrix(c(sqrt(H[1, 1] * H[2, 2]), 0, 0, sqrt(H[1, 1] * H[2, 2])), 2, 2)
  }
  if(!optimal){
    H <- matrix(c(h, 0, 0, h), 2, 2)
  }

  # Help function to compute weights
  poidsU <- function(i, U, h, POL){
    x <- as.numeric(U[i,])
    sWeights(x, h, POL)
  }
  # Use parallel methods to compute if the number of observation is a bit high
  # Change the number of slaves according to the number of cores your processor has
  # It is recommended to use a maximum of the number of cores minus one.
  if(parallel){
    if(is.null(n_clusters)) n_clusters <- detectCores()-1
    cl <- makeCluster(n_clusters)
    clusterEvalQ(cl, library("rgeos"))
    clusterEvalQ(cl, library("sp"))
    clusterExport(cl, c("sCircle", "sWeights"))
    OMEGA <- pblapply(1:n, poidsU, U = U, h = sqrt(H[1, 1]), POL = polygon, cl = cl)
    OMEGA <- do.call("c", OMEGA)
    stopCluster(cl)
  }else{
    OMEGA <- NULL
    for(i in 1:n){
      OMEGA <- c(OMEGA, poidsU(i, U, h = sqrt(H[1, 1]), POL = polygon))
    }
  }


  # Kernel Density Estimator
  fhat <- kde(U, H, w = 1/OMEGA,
              xmin = c(min(get.bbox(polygon)$x), min(get.bbox(polygon)$y)),
              xmax = c(max(get.bbox(polygon)$x), max(get.bbox(polygon)$y)))
  fhat$estimate <- fhat$estimate * sum(1/OMEGA) / n

  vx <- unlist(fhat$eval.points[1])
  vy <- unlist(fhat$eval.points[2])
  VX <- cbind(rep(vx, each = length(vy)))
  VY <- cbind(rep(vy, length(vx)))
  VXY <- cbind(VX, VY)
  Ind <- matrix(inside.gpc.poly(x = VX, y = VY, polyregion = polygon), length(vy), length(vx))
  f0 <- fhat
  f0$estimate[t(Ind) == 0] <- NA

  list(
    X = fhat$eval.points[[1]],
    Y = fhat$eval.points[[2]],
    Z = fhat$estimate,
    ZNA = f0$estimate,
    H = fhat$H,
    W = fhat$w)
}# End of sKDE()


#' Computes an estimation of the density using Kernel Density estimator,
#' without correcting for fontier effects.
#'
#' @param U data points.
#' @param polygon polygon on which points lie.
#' @param optimal if TRUE, uses Hpi() to select the optimal bandwidth.
#' @param h only if optimal=FALSE, scalar bandwidth.
#' @return Returns a list whose elements are:
#'   > X:    x coordinates at wich estimate is evaluated,
#'   > Y:    y coordinates at wich estimate is evaluated,
#'   > Z:    density estimates,
#'   > ZNA:  density estimates with NA values for points outside the polygon,
#'   > H:    bandwidth matrix,
#'   > W:    vector of weights.
#' @examples
#' data(acci)
#' smoothed_fin_nc <- sKDE_without_c(U = acci$finistere$points,
#'     polygon = acci$finistere$polygon, optimal=TRUE)
sKDE_without_c = function(U, polygon, optimal = TRUE, h = .1){
  polygon <- as(polygon, "gpc.poly")
  IND <- which(is.na(U[,1]) == FALSE)
  U <- U[IND,]
  n <- nrow(U)
  if(optimal){
    H <- Hpi(U,binned=FALSE)
    H <- matrix(c(sqrt(H[1, 1] * H[2, 2]), 0, 0, sqrt(H[1, 1] * H[2, 2])), 2, 2)
  }
  if(!optimal){
    H <- matrix(c(h, 0, 0, h), 2, 2)
  }

  # Kernel density estimator
  fhat <- kde(U, H,
              xmin = c(min(get.bbox(polygon)$x), min(get.bbox(polygon)$y)),
              xmax = c(max(get.bbox(polygon)$x), max(get.bbox(polygon)$y)))

  vx <- unlist(fhat$eval.points[1])
  vy <- unlist(fhat$eval.points[2])
  VX <- cbind(rep(vx, each = length(vy)))
  VY <- cbind(rep(vy, length(vx)))
  VXY <- cbind(VX,VY)
  Ind <- matrix(inside.gpc.poly(x = VX, y = VY, polyregion = polygon), length(vy), length(vx))
  f0 <- fhat
  f0$estimate[t(Ind) == 0] <- NA

  list(
    X = fhat$eval.points[[1]],
    Y = fhat$eval.points[[2]],
    Z = fhat$estimate,
    ZNA = f0$estimate,
    H = fhat$H,
    W = fhat$W)
}# End of sKDE_without_c()


#' Using the result obtained by the evaluation of the functions sKDE() or sKDE_without_c(),
#' the function plot_sKDE() creates a visualization of the kernel density estimates.
#' @param smooth       result from sKDE() or sKDE_without_c();
#' @param breaks       breaks for the legend (seq(min(smooth$Z)*.95,max(smooth$Z)*1.05,length=21) by default);
#' @param polygon      polygon on which data points lie;
#' @param coord        coordinates (long, lat) of data points;
#' @param alpha_coords transparency for data points (.8 by default);
#' @param size_coords  size for data points (.8 by default);
#' @param many_points  if TRUE, @coord must be the result of condense() (package bigvis). It is helpful when there are too many points to display (FALSE by default);
#' @param colContour   colour of the contour of the polygon ("white" by default);
#' @param colPoints    colour of the data points ("dodger blue" by default);
#' @param title        title (if provided) to give to the plot;
#' @param contour      if FALSE, contour are not plotted (TRUE by default);
#' @param round        round value for the legend (2 by default);
#' @param text_size    text size (22 by default).
#' @return a ggplot2 plot.
#' @examples
#' data(acci)
#' # Estimation with correction
#' smoothed_fin <- sKDE(U = acci$finistere$points, polygon = acci$finistere$polygon,
#'     optimal=TRUE, parallel = FALSE)
#' # Estimation without correction
#'     smoothed_fin_nc <- sKDE_without_c(U = acci$finistere$points, polygon = acci$finistere$polygon,
#'     optimal=TRUE)
#' p_acci_fin <- plot_sKDE(smooth = smoothed_fin,
#'     coord = acci$finistere$points,
#'     alpha_coords = .8,
#'     size_coords = 1,
#'     breaks = seq(min(smoothed_fin$ZNA, smoothed_fin_nc$ZNA,na.rm=TRUE)*.95,
#'     max(smoothed_fin$ZNA, smoothed_fin_nc$ZNA,na.rm=TRUE)*1.05, length=21),
#'     polygon = acci$finistere$polygon, round = 3, colContour = "black") +
#' ggtitle("With correction") +
#' coord_equal()
#' print(p_acci_fin)
plot_sKDE <- function(smooth, breaks, polygon, coord, alpha_coords = .8, size_coords = .8,
                      many_points = FALSE,
                      colContour="white",
                      colPoints="dodger blue", title, contour=TRUE,
                      round = 2, text_size = 22){

  # Get the right format for ggplot2
  obtenirMelt <- function(smoothed){
    res <- melt(smoothed$ZNA)
    res[,1] <- smoothed$X[res[,1]]
    res[,2] <- smoothed$Y[res[,2]]
    names(res) <- list("X","Y","ZNA")
    return(res)
  }

  smCont <- obtenirMelt(smooth)
  if(missing(breaks)) breaks <- seq(min(smooth$Z)*.95,max(smooth$Z)*1.05,length=21)
  smCont$colour <- cut(smCont[,"ZNA"],breaks=breaks,labels=round(breaks[-1],digits=round))
  smCont$colour2 <- as.character(cut(smCont[,"ZNA"],breaks=breaks,labels=rev(heat.colors(length(breaks)-1))))

  if(is.null(polygon$group)) polygon$group <- factor(1)

  P <- ggplot() +
    geom_polygon(data = polygon,  aes(x = long, y = lat, group = group),
                 fill = NA, col = "black") +
    geom_tile(aes(x = X, y = Y, fill = ZNA),
              alpha = .9, data = smCont[!is.na(smCont$ZNA),], na.rm=TRUE)


  lesLabels <- round(breaks,round)
  lesIndicesLabels <- floor(seq(1,length(lesLabels),length.out=5)) # Only keep 5 points for the legend values
  lesIndicesLabels[length(lesIndicesLabels)] <- length(lesLabels) # Making sure we display the last value
  lesLabels <- as.character(lesLabels[lesIndicesLabels])
  lesLabels[lesLabels=="0"] <- "0.00"

  if(contour) P <- P + geom_contour(data = smCont[!is.na(smCont$ZNA),],
                                    aes(x = X, y = Y, z = ZNA),
                                    alpha=0.6,  colour = colContour,
                                    breaks = breaks[lesIndicesLabels])
  if(many_points){
    P <- P + geom_count(data = coord, aes(x = long, y = lat, alpha = ..prop..),
                        col = "blue", size = size_coords) +
      scale_alpha_continuous(guide=FALSE)
  }else{
    P <- P + geom_point(data = coord[,c("long", "lat")], aes(x = long, y = lat),
                        alpha = alpha_coords, col = "blue", size = size_coords)
  }


  if(contour){
    # To add contour levels
    ind_level <- which(unlist(lapply(ggplot_build(P)$data, function(x) "level" %in% colnames(x))))
    tmp <- ggplot_build(P)$data[[ind_level]]
    ind <- unlist(lapply(unique(tmp$piece), function(x){
      corresp <- which(tmp$piece == x)
      corresp[round(length(corresp)/2)]
    }))
    tmp$level_r <- round(tmp$level, round)
    P <- P + geom_text(aes(label = level_r, x = x, y = y), data=tmp[ind,])
  }

  P <- P + scale_fill_gradient(name="",low='yellow', high='red',
                               breaks=breaks[lesIndicesLabels],
                               limits=range(breaks),labels=lesLabels)

  P <- P + theme(axis.text.x=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks.x=element_blank(),
                 axis.ticks.y=element_blank(),
                 axis.title=element_blank(),
                 text = element_text(size = text_size))

  P <- P + geom_polygon(data=polygon, mapping=(aes(x=long, y=lat)),
                        colour="black", fill=NA)
  # Add a title if one was provided
  if(!missing(title)) P <- P + ggtitle(title)
  P
}
