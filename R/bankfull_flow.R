#' Bankfull Flow Calculation
#'
#' Calculates the bankfull flow for a given cross-section. Manning's equation is used to
#' solve for the flow at a given stage, slope and cross-section properties.
#'
#' The wsl parameter sets the water surface level for the flow calculation. If not provided, the function
#' will take the lower of the two bank stations at the bankfull elevation and solve for flow at that stage.
#' Note that in this context, the default calculation will be for 'bankfull' flow, however, any water surface
#' level can be inputted to this function to determine the flow at that stage.
#'
#' A composite Manning's n value is used in estimation of the roughness value for the cross-section, which is used in
#' Manning's equation.
#'
#'
#' @param xs cross-section of object class xsection
#' @param S bed slope in cross-section for use in Manning's equation (m/m)
#' @param wsl (optional) water surface level to set as bankfull for flow calculation (m)
#' @param options object of class rcr_options with options and constants for hydraulic calculations
#'
#' @return \item{flow_stage}{data frame with flow and stage data, cms and m, respectively}
#' @seealso \code{\link{normal_depth}} to estimate the normal depth using Manning's equation
#' @seealso \code{\link{flow_area}} to calculate the flow area in cross-section for given segments
#'
#' @keywords bankfull flow Manning stage
#' @examples
#'
#' @export bankfull_flow


bankfull_flow <- function(xs,S,wsl=NA,options) {

  if (class(xs)[1] != "xsection") {
    stop("xs must be of class xsection.")
  }
  if (is.na(S)) {
    stop("S is a required parameter.")
  }
  if (class(S) != "numeric") {
    stop("S must be of class numeric.")
  }
  if (!(is.na(wsl))) {
    if (wsl < 0 ) {
      stop("wsl must be a positive value")
    }
    else if (wsl < min(xs$zz)) {
      warning("wsl is less than the minimum elevation in the channel, flow will be estimated as zero.")
    }
  }
  if (class(options)[1] != "rcr_options") {
    stop("options object must be of class rcr_options.")
  }
  if (options$dx <= 0 | class(options$dx) != "numeric") {
    stop("dx in options must be a numeric value greater than 0.")
  }

  mm <- data.frame(matrix(NA,nrow=1,ncol=37))
  colnames(mm) <- c("xsection","Flow","Flow_LOB","Flow_Main","Flow_ROB","Min_Elev","Depth","WSL",
                    "Velocity","Velocity_LOB","Velocity_Main","Velocity_ROB","K_Total",
                    "K_LOB","K_Main","K_ROB","alpha","Area","Area_LOB","Area_Main","Area_ROB",
                    "HRadius","HRadius_LOB","HRadius_Main","HRadius_ROB","WetPerimeter","WetPerimeter_LOB",
                    "WetPerimeter_Main","WetPerimeter_ROB","Energy_total","Velocity_head","Froude",
                    "Sf","Sf_Avg","Length_Effective","Head_Loss","Manning_Composite")

  ### Calculation routine
  i = 1

  mm[i,]$Min_Elev <- min(xs$zz)

  ## interpolate series (for computational purposes)
  min_elev <- min(xs$zz)
  min_dist <- min(xs$xx)
  max_dist <- max(xs$xx) # take shorter of the two series
  xx <- seq(min_dist,max_dist,options$dx)
  N <- length(xx)
  zz <- approx(x=xs$xx,y=xs$zz,xout=xx)$y

  # determine WSL based on input elevation or minimum of bank heights
  if (!(is.na(wsl))) {
    mm[i,]$WSL <- wsl
  } else {
    mm[i,]$WSL <- min(c(approx(x=xx,y=zz,xout=xs$lbs_xx)$y,approx(x=xx,y=zz,xout=xs$rbs_xx)$y))
  }

  mm[i,]$Depth <- mm[i,]$WSL - mm[i,]$Min_Elev

  # continue interpolation
  temp <- mm[i,]$WSL - zz
  ind <- which(temp>0)
  xx.df <- data.frame(xx,zz,temp)

  mm[i,]$Area <- as.numeric(sum(temp[which(temp>0)])*options$dx)
  mm[i,]$Area_Main <- sum(xx.df[which(xx.df$temp>0 & xx.df$xx > xs$lbs_xx & xx.df$xx < xs$rbs_xx),]$temp)*options$dx
  mm[i,]$Area_LOB <- sum(xx.df[which(xx.df$temp>0 & xx.df$xx < xs$lbs_xx),]$temp)*options$dx
  mm[i,]$Area_ROB <- sum(xx.df[which(xx.df$temp>0 & xx.df$xx > xs$rbs_xx),]$temp)*options$dx

  mm[i,]$WetPerimeter <- wetted_perimeter(zz,ind,options$dx)
  zz_lob <- xx.df[which(xx.df$xx < xs$lbs_xx),]$zz
  ind_lob <- which(xx.df[which(xx.df$xx < xs$lbs_xx),]$temp > 0)
  mm[i,]$WetPerimeter_LOB <- wetted_perimeter(zz_lob,ind_lob,options$dx)
  zz_rob <- xx.df[which(xx.df$xx > xs$rbs_xx),]$zz
  ind_rob <- which(xx.df[which(xx.df$xx > xs$rbs_xx),]$temp > 0)
  mm[i,]$WetPerimeter_ROB <- wetted_perimeter(zz_rob,ind_rob,options$dx)

  zz_main <- xx.df[which(xx.df$xx >= xs$lbs_xx & xx.df$xx <= xs$rbs_xx),]$zz
  ind_main <- which(xx.df[which(xx.df$xx >= xs$lbs_xx & xx.df$xx <= xs$rbs_xx),]$temp > 0)
  mm[i,]$WetPerimeter_Main <- wetted_perimeter(zz_main,ind_main,options$dx)

  mm[i,]$HRadius <- max(mm[i,]$Area/mm[i,]$WetPerimeter,0,na.rm=T)
  mm[i,]$HRadius_Main <- max(mm[i,]$Area_Main/mm[i,]$WetPerimeter_Main,0,na.rm=T)
  mm[i,]$HRadius_LOB <- max(mm[i,]$Area_LOB/mm[i,]$WetPerimeter_LOB,0,na.rm=T)
  mm[i,]$HRadius_ROB <- max(mm[i,]$Area_ROB/mm[i,]$WetPerimeter_ROB,0,na.rm=T)

  mm[i,]$Manning_Composite <- (sum(mm[i,]$WetPerimeter_LOB*xs$Manning_LOB^1.5,mm[i,]$WetPerimeter_Main*xs$Manning_Main^1.5,mm[i,]$WetPerimeter_ROB*xs$Manning_ROB^1.5) /
                                 mm[i,]$WetPerimeter)^(2/3)

  mm[i,]$Flow <- (1/mm[i,]$Manning_Composite)*mm[i,]$Area*mm[i,]$HRadius^(2/3)*S^(1/2)

  return(data.frame("Flow"=mm[i,]$Flow,"WSL"=mm[i,]$WSL))

}

