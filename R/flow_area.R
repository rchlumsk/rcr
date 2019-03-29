#' Flow Area Calculation
#'
#' Calculates the cross-sectional flow area
#'
#' Calculates the flow area for a set of specified areas within the cross-section.
#'
#' @param xs an object of class xsection
#' @param wsl an elevation at which the water level exists in the section
#' @param stns station points along the horizontal cross-section to divide the flow area. If none is supplied, the left and right bank stations only will be used.
#' @param options object of class rcr_options with options and constants for hydraulic calculations
#'
#' @return \item{flow_areas}{Table of computed cross-sectional flow area, m2}
#' @seealso \code{\link{bankfull_estimator}} to estimate the flow through the cross section using Manning's equation
#'
#' @keywords compute flow area
#' @examples
#'
#' @export flow_area

flow_area <- function(xs,wsl,stns=NA,options) {

  ## Check boundary condition
  if (class(xs)[1] != "xsection") {
    stop("xsection must be of type xsection")
  }
  if (wsl <= min(xs$zz)) {
    warning("wsl provided is less than the minimum cross-sectional elevation. All flow areas will be zero.")
  }

  if (is.na(stns)) {
    stns <- c(xs$lbs_xx,xs$rbs_xx)
  } else {
    if (!(is.numeric(stns)) | length(stns) < 1) {
      stop("stns must be NA or supplied as a numeric vector of at least one value.")
    }
  }

  min_elev <- min(xs$zz)
  depth <- wsl - min_elev

  min_dist <- min(xs$xx)
  max_dist <- max(xs$xx) # take shorter of the two series
  xx <- seq(min_dist,max_dist,options$dx)
  N <- length(xx)
  zz <- approx(x=xs$xx,y=xs$zz,xout=xx)$y
  temp <- wsl - zz
  ind <- which(temp>0)
  xx.df <- data.frame(xx,zz,temp)

  if (stns[1] != min(xs$xx)) {
    stns <- c(min(xs$xx),stns)
  }
  if (stns[length(stns)] != max(xs$xx)) {
    stns <- c(stns,max(xs$xx))
  }

  mm <- as.data.frame(matrix(NA,ncol=4,nrow=(length(stns)-1)))
  colnames(mm) <- c("ID","start","end","area")
  mm$ID <- seq(1,nrow(mm))
  mm$start <- stns[1:(length(stns)-1)]
  mm$end <- stns[2:length(stns)]

  for (i in 1:nrow(mm)) {
    mm[i,]$area <- sum(xx.df[which(xx.df$temp>0 & xx.df$xx > mm[i,]$start & xx.df$xx < mm[i,]$end),]$temp)*options$dx
  }

  print("Successfully completed flow areas. :-)")
  return("flow_areas"=mm)
}
