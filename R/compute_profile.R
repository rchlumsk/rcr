#' Hydraulic Profile Calculation
#'
#' Calculates the water profile using the standard step method.
#'
#' **Currently assumes subcritical regime throughout the profile, will be updated in future versions.**
#'
#' @param geometry geometry object of class geom
#' @param boundary_conditions boundary conditions for the most downstream cross-section of object bc
#' @param method method used in calculation as "subcritical", "supercritical", or "mixed" (currently only "subcritical" method is supported)
#' @param options object of class rcr_options with options and constants for hydraulic calculations
#'
#' @return \item{hydraulic_output}{Table of computed hydraulic parameters by xsection}
#'
#' @keywords compute profile HEC-RAS
#' @examples
#'
#' @export compute_profile

compute_profile <- function(geometry,Q,boundary_conditions,method="subcritical",options) {


  # geometry: object of class geom with >1 xsection objects in xsectionList

  # Q <- flow value (cms), to be replaced with flow profile object in future iterations

  # method: method for computation, one of subcritical, supercritical, mixed
  #   currently only subcritical is implemented

  #boundary_conditions: object of class bc with boundary conditions supplied for geom

  # options: object of class rcr_options with options for computation embedded (iteration limit, output amount, etc.)




  ## need lots of error checking here
  # all xsection has valid data points, non-empty xx and zz
  # all parameters filled out
  # boundary conditions make sense
  # etc etc

  # check method
  if (method != "subcritical") {
    stop("Only subcritical mode with hardcoded options currently available.")
  }

  # check geometry
  if (length(geometry$xsectionList) <= 1) {
    stop("Calculation requires at least 2 cross-sections in the geometry object.")
  }

  ## Check boundary condition
  if (class(boundary_conditions)[1] != "bc") {
    stop("boundary_conditions must be of type bc")
  }
  if (boundary_conditions$bctype != "Normal Depth") {
    stop("Boundary condition must be of type 'Normal Depth'.")
  }
  if (is.na(boundary_conditions$bcvalue) | boundary_conditions$bcvalue <= 0) {
    stop("Boundary condition value must be slope of value greater than 0.")
  }


  # check number of items in geometry list
  mm <- data.frame(matrix(NA,nrow=length(geometry$xsectionList),ncol=37))
  colnames(mm) <- c("xsection","Flow","Flow_LOB","Flow_Main","Flow_ROB","Min_Elev","Depth","WSL",
                    "Velocity","Velocity_LOB","Velocity_Main","Velocity_ROB","K_Total",
                    "K_LOB","K_Main","K_ROB","alpha","Area","Area_LOB","Area_Main","Area_ROB",
                    "HRadius","HRadius_LOB","HRadius_Main","HRadius_ROB","WetPerimeter","WetPerimeter_LOB",
                    "WetPerimeter_Main","WetPerimeter_ROB","Energy_total","Velocity_head","Froude",
                    "Sf","Sf_Avg","Length_Effective","Head_Loss","Manning_Composite")

  # check options
  if (class(options)[1] != "rcr_options") {
    stop("options object must be of class rcr_options.")
  }
  if (options$dx <= 0 | class(options$dx) != "numeric") {
    stop("dx in options must be a numeric value greater than 0.")
  }

  # update - get xsection from each item in xsectionList
  for (i in 1:length(geometry$xsectionList)) {
    mm[i,]$xsection <- geometry$xsectionList[[i]]$riverstation
  }

  ### Calculation routine
  if (!(options$silent_cp)) {
    print("Beginning backwater calculations.")
  }


  # parameters at downstream Xsection

  i = 1

  if (!(options$silent_cp)) {
    print(sprintf("Computing profile for xsection %i",i))
  }

  mm$Flow <- Q
  mm[i,]$Min_Elev <- min(geometry$xsectionList[[i]]$zz)
  mm[i,]$WSL <- rcr::normal_depth(geometry$xsectionList[[i]], mm[i,]$Flow, boundary_conditions$bcvalue,options = options) # assume boundary condition as normal depth
  mm[i,]$Depth <- mm[i,]$WSL - mm[i,]$Min_Elev

  ## interpolate series (for computational purposes)
  min_elev <- min(geometry$xsectionList[[i]]$zz)
  min_dist <- min(geometry$xsectionList[[i]]$xx)
  max_dist <- max(geometry$xsectionList[[i]]$xx) # take shorter of the two series
  xx <- seq(min_dist,max_dist,options$dx)
  N <- length(xx)
  zz <- approx(x=geometry$xsectionList[[i]]$xx,y=geometry$xsectionList[[i]]$zz,xout=xx)$y
  temp <- mm[i,]$WSL - zz
  ind <- which(temp>0)
  xx.df <- data.frame(xx,zz,temp)

  mm[i,]$Area <- as.numeric(sum(temp[which(temp>0)])*options$dx)
  mm[i,]$Area_Main <- sum(xx.df[which(xx.df$temp>0 & xx.df$xx > geometry$xsectionList[[i]]$lbs_xx & xx.df$xx < geometry$xsectionList[[i]]$rbs_xx),]$temp)*options$dx
  mm[i,]$Area_LOB <- sum(xx.df[which(xx.df$temp>0 & xx.df$xx < geometry$xsectionList[[i]]$lbs_xx),]$temp)*options$dx
  mm[i,]$Area_ROB <- sum(xx.df[which(xx.df$temp>0 & xx.df$xx > geometry$xsectionList[[i]]$rbs_xx),]$temp)*options$dx

  mm[i,]$WetPerimeter <- wetted_perimeter(zz,ind,options$dx)
  zz_lob <- xx.df[which(xx.df$xx < geometry$xsectionList[[i]]$lbs_xx),]$zz
  ind_lob <- which(xx.df[which(xx.df$xx < geometry$xsectionList[[i]]$lbs_xx),]$temp > 0)
  mm[i,]$WetPerimeter_LOB <- wetted_perimeter(zz_lob,ind_lob,options$dx)
  zz_rob <- xx.df[which(xx.df$xx > geometry$xsectionList[[i]]$rbs_xx),]$zz
  ind_rob <- which(xx.df[which(xx.df$xx > geometry$xsectionList[[i]]$rbs_xx),]$temp > 0)
  mm[i,]$WetPerimeter_ROB <- wetted_perimeter(zz_rob,ind_rob,options$dx)

  zz_main <- xx.df[which(xx.df$xx >= geometry$xsectionList[[i]]$lbs_xx & xx.df$xx <= geometry$xsectionList[[i]]$rbs_xx),]$zz
  ind_main <- which(xx.df[which(xx.df$xx >= geometry$xsectionList[[i]]$lbs_xx & xx.df$xx <= geometry$xsectionList[[i]]$rbs_xx),]$temp > 0)
  mm[i,]$WetPerimeter_Main <- wetted_perimeter(zz_main,ind_main,options$dx)

  mm[i,]$HRadius <- max(mm[i,]$Area/mm[i,]$WetPerimeter,0,na.rm=T)
  mm[i,]$HRadius_Main <- max(mm[i,]$Area_Main/mm[i,]$WetPerimeter_Main,0,na.rm=T)
  mm[i,]$HRadius_LOB <- max(mm[i,]$Area_LOB/mm[i,]$WetPerimeter_LOB,0,na.rm=T)
  mm[i,]$HRadius_ROB <- max(mm[i,]$Area_ROB/mm[i,]$WetPerimeter_ROB,0,na.rm=T)

  mm[i,]$K_LOB <- conv_calc(geometry$xsectionList[[i]]$Manning_LOB,mm[i,]$Area_LOB,mm[i,]$HRadius_LOB)
  mm[i,]$K_Main <- conv_calc(geometry$xsectionList[[i]]$Manning_Main,mm[i,]$Area_Main,mm[i,]$HRadius_Main)
  mm[i,]$K_ROB <- conv_calc(geometry$xsectionList[[i]]$Manning_ROB,mm[i,]$Area_ROB,mm[i,]$HRadius_ROB)
  mm[i,]$K_Total <- sum(mm[i,]$K_LOB,mm[i,]$K_Main,mm[i,]$K_ROB) # need to sum up all three conveyances to get total

  mm[i,]$Flow_LOB <- mm[i,]$Flow*mm[i,]$K_LOB/mm[i,]$K_Total
  mm[i,]$Flow_Main <- mm[i,]$Flow*mm[i,]$K_Main/mm[i,]$K_Total
  mm[i,]$Flow_ROB <- mm[i,]$Flow*mm[i,]$K_ROB/mm[i,]$K_Total

  mm[i,]$Velocity <- mm[i,]$Flow/mm[i,]$Area
  mm[i,]$Velocity_LOB <- max(mm[i,]$Flow_LOB/mm[i,]$Area_LOB,0,na.rm=T)
  mm[i,]$Velocity_Main <- max(mm[i,]$Flow_Main/mm[i,]$Area_Main,0,na.rm=T)
  mm[i,]$Velocity_ROB <- max(mm[i,]$Flow_ROB/mm[i,]$Area_ROB,0,na.rm=T)

  mm[i,]$alpha <- (mm[i,]$Flow_LOB*mm[i,]$Velocity_LOB^2 + mm[i,]$Flow_Main*mm[i,]$Velocity_Main^2 + mm[i,]$Flow_ROB*mm[i,]$Velocity_ROB^2) / (mm[i,]$Flow*mm[i,]$Velocity^2)
  mm[i,]$Velocity_head <- vhead_calc(mm[i,]$alpha,mm[i,]$Velocity,options$g)
  mm[i,]$Energy_total <- mm[i,]$Velocity_head + mm[i,]$WSL
  mm[i,]$Froude <- froude_calc(mm[i,]$Velocity,options$g,mm[i,]$Depth)

  mm[i,]$Manning_Composite <- (sum(mm[i,]$WetPerimeter_LOB*geometry$xsectionList[[i]]$Manning_LOB^1.5,mm[i,]$WetPerimeter_Main*geometry$xsectionList[[i]]$Manning_Main^1.5,
                                   mm[i,]$WetPerimeter_ROB*geometry$xsectionList[[i]]$Manning_ROB^1.5) / mm[i,]$WetPerimeter)^(2/3)

  mm[i,]$Sf <- (mm[i,]$Flow/mm[i,]$K_Total)^2


  # estimate upstream section

  for (i in 2:length(geometry$xsectionList)) {

    if (!(options$silent_cp)) {
      print(sprintf("Computing profile for xsection %i",i))
    }

    # initial setting of items
    mm[i,]$Min_Elev <- min(geometry$xsectionList[[i]]$zz)
    mm[i,]$WSL <- mm[i-1,]$Depth + mm[i,]$Min_Elev # best guess at upstream WSL, same depth applied to bottom bed elevation
    mm[i,]$Depth <- mm[i,]$WSL - mm[i,]$Min_Elev

    for (j in 1:options$iteration_limit_cp) {

      mm[i,]$Depth <- mm[i,]$WSL - mm[i,]$Min_Elev

      ## interpolate series (for computational purposes)
      min_elev <- min(geometry$xsectionList[[i]]$zz)
      min_dist <- min(geometry$xsectionList[[i]]$xx)
      max_dist <- max(geometry$xsectionList[[i]]$xx) # take shorter of the two series
      xx <- seq(min_dist,max_dist,options$dx)
      N <- length(xx)
      zz <- approx(x=geometry$xsectionList[[i]]$xx,y=geometry$xsectionList[[i]]$zz,xout=xx)$y
      temp <- mm[i,]$WSL - zz
      ind <- which(temp>0)
      xx.df <- data.frame(xx,zz,temp)

      mm[i,]$Area <- as.numeric(sum(temp[which(temp>0)])*options$dx)
      mm[i,]$Area_Main <- sum(xx.df[which(xx.df$temp>0 & xx.df$xx > geometry$xsectionList[[i]]$lbs_xx & xx.df$xx < geometry$xsectionList[[i]]$rbs_xx),]$temp)*options$dx
      mm[i,]$Area_LOB <- sum(xx.df[which(xx.df$temp>0 & xx.df$xx < geometry$xsectionList[[i]]$lbs_xx),]$temp)*options$dx
      mm[i,]$Area_ROB <- sum(xx.df[which(xx.df$temp>0 & xx.df$xx > geometry$xsectionList[[i]]$rbs_xx),]$temp)*options$dx

      mm[i,]$WetPerimeter <- wetted_perimeter(zz,ind,options$dx)
      zz_lob <- xx.df[which(xx.df$xx < geometry$xsectionList[[i]]$lbs_xx),]$zz
      ind_lob <- which(xx.df[which(xx.df$xx < geometry$xsectionList[[i]]$lbs_xx),]$temp > 0)
      mm[i,]$WetPerimeter_LOB <- wetted_perimeter(zz_lob,ind_lob,options$dx)
      zz_rob <- xx.df[which(xx.df$xx > geometry$xsectionList[[i]]$rbs_xx),]$zz
      ind_rob <- which(xx.df[which(xx.df$xx > geometry$xsectionList[[i]]$rbs_xx),]$temp > 0)
      mm[i,]$WetPerimeter_ROB <- wetted_perimeter(zz_rob,ind_rob,options$dx)

      zz_main <- xx.df[which(xx.df$xx >= geometry$xsectionList[[i]]$lbs_xx & xx.df$xx <= geometry$xsectionList[[i]]$rbs_xx),]$zz
      ind_main <- which(xx.df[which(xx.df$xx >= geometry$xsectionList[[i]]$lbs_xx & xx.df$xx <= geometry$xsectionList[[i]]$rbs_xx),]$temp > 0)
      mm[i,]$WetPerimeter_Main <- wetted_perimeter(zz_main,ind_main,options$dx)

      mm[i,]$HRadius <- max(mm[i,]$Area/mm[i,]$WetPerimeter,0,na.rm=T)
      mm[i,]$HRadius_Main <- max(mm[i,]$Area_Main/mm[i,]$WetPerimeter_Main,0,na.rm=T)
      mm[i,]$HRadius_LOB <- max(mm[i,]$Area_LOB/mm[i,]$WetPerimeter_LOB,0,na.rm=T)
      mm[i,]$HRadius_ROB <- max(mm[i,]$Area_ROB/mm[i,]$WetPerimeter_ROB,0,na.rm=T)

      mm[i,]$K_LOB <- conv_calc(geometry$xsectionList[[i]]$Manning_LOB,mm[i,]$Area_LOB,mm[i,]$HRadius_LOB)
      mm[i,]$K_Main <- conv_calc(geometry$xsectionList[[i]]$Manning_Main,mm[i,]$Area_Main,mm[i,]$HRadius_Main)
      mm[i,]$K_ROB <- conv_calc(geometry$xsectionList[[i]]$Manning_ROB,mm[i,]$Area_ROB,mm[i,]$HRadius_ROB)
      mm[i,]$K_Total <- sum(mm[i,]$K_LOB,mm[i,]$K_Main,mm[i,]$K_ROB) # need to sum up all three conveyances to get total

      mm[i,]$Flow_LOB <- mm[i,]$Flow*mm[i,]$K_LOB/mm[i,]$K_Total
      mm[i,]$Flow_Main <- mm[i,]$Flow*mm[i,]$K_Main/mm[i,]$K_Total
      mm[i,]$Flow_ROB <- mm[i,]$Flow*mm[i,]$K_ROB/mm[i,]$K_Total

      mm[i,]$Velocity <- mm[i,]$Flow/mm[i,]$Area
      mm[i,]$Velocity_LOB <- max(mm[i,]$Flow_LOB/mm[i,]$Area_LOB,0,na.rm=T)
      mm[i,]$Velocity_Main <- max(mm[i,]$Flow_Main/mm[i,]$Area_Main,0,na.rm=T)
      mm[i,]$Velocity_ROB <- max(mm[i,]$Flow_ROB/mm[i,]$Area_ROB,0,na.rm=T)

      mm[i,]$alpha <- (mm[i,]$Flow_LOB*mm[i,]$Velocity_LOB^2 + mm[i,]$Flow_Main*mm[i,]$Velocity_Main^2 + mm[i,]$Flow_ROB*mm[i,]$Velocity_ROB^2) / (mm[i,]$Flow*mm[i,]$Velocity^2)
      mm[i,]$Velocity_head <- vhead_calc(mm[i,]$alpha,mm[i,]$Velocity,options$g)
      mm[i,]$Energy_total <- mm[i,]$Velocity_head + mm[i,]$WSL
      mm[i,]$Froude <- froude_calc(mm[i,]$Velocity,options$g,mm[i,]$Depth)

      mm[i,]$Manning_Composite <- (sum(mm[i,]$WetPerimeter_LOB*geometry$xsectionList[[i]]$Manning_LOB^1.5,mm[i,]$WetPerimeter_Main*geometry$xsectionList[[i]]$Manning_Main^1.5,
                                       mm[i,]$WetPerimeter_ROB*geometry$xsectionList[[i]]$Manning_ROB^1.5) / mm[i,]$WetPerimeter)^(2/3)

      # calculate friction slope
      mm[i,]$Sf <- (mm[i,]$Flow/mm[i,]$K_Total)^2
      mm[i,]$Sf_Avg <- mean(mm[i,]$Sf,mm[i-1,]$Sf)

      # calculate head loss, he
      loss_coeff <- 0
      if (mm[i-1,]$Velocity > mm[i,]$Velocity) {
        # downstream velocity greater than upstream  velocity
        #   water contracting moving downstream, use contraction coeff
        loss_coeff <- geometry$xsectionList[[i]]$contraction_coeff
      } else {
        # downstream velocity less than upstream velocity,
        #   water expanding moving downstream, use expansion coeff
        loss_coeff <- geometry$xsectionList[[i]]$expansion_coeff
      }

      mm[i,]$Length_Effective <- (geometry$xsectionList[[i]]$ds_length_LOB*mm[i,]$Flow_LOB + geometry$xsectionList[[i]]$ds_length_Main*mm[i,]$Flow_Main + geometry$xsectionList[[i]]$ds_length_ROB*mm[i,]$Flow_ROB)/ (mm[i,]$Flow_LOB + mm[i,]$Flow_Main + mm[i,]$Flow_ROB)

      mm[i,]$Head_Loss <- mm[i,]$Length_Effective*mm[i,]$Sf_Avg + loss_coeff*abs(mm[i,]$alpha*mm[i,]$Velocity^2/2/options$g - mm[i-1,]$alpha*mm[i-1,]$Velocity^2/2/options$g)

      next_WSL <- mm[i-1,]$WSL + mm[i-1,]$Velocity_head + mm[i,]$Head_Loss - mm[i,]$Velocity_head

      if ( abs(next_WSL - mm[i,]$WSL) > options$tolerance_cp) {

        # print(sprintf("out of tolerance %s, %.4f assumed %.4f computed",j,mm[i,]$WSL,next_WSL))

        if (j > options$iteration_limit_cp) {
          warning(sprintf("Iteration limit %s reached, terminating iterations on xsection %s",options$iteration_limit_cp,i))
          break
        }

        if ( j == 1) {
          mm[i,]$WSL <- next_WSL
        } else if (j == 2 ) {
          mm[i,]$WSL <- mm[i,]$WSL + options$next_WSL_split_cp*(next_WSL - mm[i,]$WSL)
        } else {
          mm[i,]$WSL <- mm[i,]$WSL + 0.5*(next_WSL - mm[i,]$WSL)
        }

      } else {
        if (!(options$silent_cp)) {
          print(sprintf("Iterated on WSL within tolerance on iteration %s",j))
        }
        break
      }
    }
  }

  print("Successfully completed hydraulic calculations. :-)")
  return("hydraulic_output"=mm)
}
