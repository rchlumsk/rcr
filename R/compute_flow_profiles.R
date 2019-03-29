#' Hydraulic Profile Calculation for Multiple Flows
#'
#' Wrapper for the compute_profile function to run with multiple flow values and return an organized data frame
#'
#' **Currently assumes subcritical regime throughout the profile, will be updated in future versions.**
#'
#' @param geometry geometry object of class geom
#' @param flows vector of flow values for computation (cms)
#' @param boundary_conditions boundary conditions for the most downstream cross-section of object bc
#' @param method method used in calculation as "subcritical", "supercritical", or "mixed" (currently only "subcritical" method is supported)
#' @param options object of class rcr_options with options and constants for hydraulic calculations
#'
#' @return \item{hydraulic_output}{Table of computed hydraulic parameters by xsection}
#' @seealso \code{\link{compute_profile}} to calculate the water profile for a single flow value
#'
#' @keywords compute profiles flow
#' @examples
#'
#' @export compute_flow_profiles

compute_flow_profiles <- function(geometry,flows,boundary_conditions,method="subcritical",options) {

  if (!(is.numeric(flows)) | length(flows) < 1) {
    stop("Flows must be a numeric vector with at least one value.")
  }

  mm1 <- compute_profile(geometry,flows[1],boundary_conditions,method="subcritical",options)

  if (length(flows) > 1) {

    for (i in 2:length(flows)) {
      mm2 <- compute_profile(geometry,flows[i],boundary_conditions,method="subcritical",options)
      mm1 <- rbind(mm1,mm2)
    }
  }

  print(sprintf("Successfully completed hydraulic calculations for %i flow profiles. :-)",length(flows)))
  return("hydraulic_output"=mm1)
}
