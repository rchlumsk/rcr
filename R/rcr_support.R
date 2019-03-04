#' Support functions for hydraulic calculations in the rcr package
#'
#' Support functions used in the rcr package for calculating various hydraulic properties.
#' Includes calculations for conveyance, friction slope, energy calculation, and others.
#'
#' Note that some functions may have multiple versions for calcualting a hydraulic parameter,
#' such as conveyance or friction slope. To view the function calculation, type the function name
#' into the R console with no brackets to view the function code (e.g., 'energy_calc').
#'
#' @param Z water surface elevation (m)
#' @param y water depth (m)
#' @param v water average velocity (m/s)
#' @param n Manning's n value (weighted for channel)
#' @param Rh hydraulic radius (m)
#' @param Q flow (m3/s)
#' @param K channel conveyance
#'
#' Other parameters and various outputs as per function usage.
#'
#' @keywords rcr hydraulic properties functions calculations
#' @name rcr_support
energy_calc <- function(Z,y,v) {
  # total energy in the channel
  return(Z+y+v^2/2/9.81)
}

#' @rdname rcr_support
sf_calc <- function(n,v,Rh) {
  # friction slope from Manning's equation
  return((n*v*Rh^(3/2))^2)
}

#' @rdname rcr_support
sf_calc2 <- function(Q,K) {
  # friction slope from conveyance and flow
  return((Q/K)^2)
}

#' @rdname rcr_support
conv_calc <- function(n,A,Rh) {
  # conveyance
  return(max((1/n)*A*Rh^(2/3),0,na.rm=T))
}

#' @rdname rcr_support
vhead_calc <- function(alpha,v,g) {
  # velocity head
  return(alpha*v^2/2/g^2)
}

#' @rdname rcr_support
wetted_perimeter <- function(zz,ind,dx) {
  # wetted perimeter
  # zz, ind, dx provided within other functions - see compute_profile or normal_depth for an example
  return(max(sum(((zz[ind][-1] - zz[ind][-length(ind)])^2 + dx^2)^0.5),0,na.rm=T))
}

#' @rdname rcr_support
froude_calc <- function(v,g,d) {
  # froude number
  return(v/sqrt(g*d))
}
