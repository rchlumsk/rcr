#' Normal Depth Calculation
#'
#' Calculates the normal depth for a given cross-section using the Manning's equation.
#'
#' @param xs cross-section of object class xsection
#' @param Q flow value through cross-section (m3/s)
#' @param S bed slope in cross-section, typically supplied as part of a boundary condition
#' @param init_WSL initial water surface level provided to initialize calculation
#' @param options object of class rcr_options with options and constants for hydraulic calculations
#'
#' @return \item{WSL}{Water surface elevation}
#' @seealso \code{\link{bankfull_estimator}} to estimate the flow through the cross section using Manning's equation
#' @seealso \code{\link{flow_area}} to calculate the flow area in cross-section for given segments
#'
#' The Manning's n used in the calculation of the normal depth is calculated as a composite Manning's n,
#' weighted using the wetted pereimeter of the left overbank, main channel, and right overbank.
#' @keywords normal depth
#' @examples
#'
#' @export normal_depth


normal_depth <- function(xs,Q,S,init_WSL=NA,options) {

  if (class(xs)[1] != "xsection") {
    stop("xs must be of class xsection.")
  }
  if (is.na(Q)) {
    stop("Q is a required parameter.")
  }
  if (class(Q) != "numeric") {
    stop("Q must be of class numeric.")
  }
  if (is.na(S)) {
    stop("S is a required parameter.")
  }
  if (class(S) != "numeric") {
    stop("S must be of class numeric.")
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

  mm$Flow <- Q
  mm[i,]$Min_Elev <- min(xs$zz)


  if (is.na(init_WSL)) {
    mm[i,]$WSL <- mm[i,]$Min_Elev + 1 # initial point, not sure what is a better starting point
  } else {
    # print("Initializing normal depth with initial WSL.")
    # print(init_WSL)
    mm[i,]$WSL <- init_WSL
  }


  for (j in 1:(options$iteration_limit_nd)) {

    mm[i,]$Depth <- mm[i,]$WSL - mm[i,]$Min_Elev

    ## interpolate series (for computational purposes)
    min_elev <- min(xs$zz)
    min_dist <- min(xs$xx)
    max_dist <- max(xs$xx) # take shorter of the two series
    xx <- seq(min_dist,max_dist,options$dx)
    N <- length(xx)
    zz <- approx(x=xs$xx,y=xs$zz,xout=xx)$y
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

    RHS <- (1/mm[i,]$Manning_Composite)*mm[i,]$Area*mm[i,]$HRadius^(2/3)*S^(1/2)

    if (abs(Q - RHS) > options$tolerance_nd) {

      # terminate if iteration limit reached
      if (j >= options$iteration_limit_nd) {
        warning("Iteration limit on normal depth calculation exceeded, terminating normal depth calculation.")
        return("WSL"=mm[i,]$WSL)
      }

      # estimate next WSL as ratio of RHS/Q to modify the depth in the channel by a portion of that ratio
      mm[i,]$WSL <- mm[i,]$Min_Elev + mm[i,]$Depth + mm[i,]$Depth*((Q-RHS)/Q)*options$next_WSL_split_nd

    } else {
      if (!(options$silent_nd)) {
        print("Normal depth estimated successfully.")
      }
      return("WSL"=mm[i,]$WSL)
    }
  }
}

