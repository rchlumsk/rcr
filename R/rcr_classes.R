#' Cross-Section Object
#'
#' xsection:
#' A reference class to create a cross-section with all relevant data and static hydraulic properties.
#'
#' @field riverstation name of the cross-section as character
#' @field station station location along profile as numeric value
#' @field xx horizontal coordinates of 1D section
#' @field zz vertical elevation coordinates of 1D section, corresponding to xx
#' @field Manning vector of Manning's n values corresponding to xx (currently not used in rcr)
#' @field Manning_LOB Manning's n value for the left overbank area
#' @field Manning_Main Manning's n value for the main channel area
#' @field Manning_ROB Manning's n value for the right overbank area
#' @field lbs_xx horizontal co-orindate of the left bank station
#' @field rbs_xx horizontal co-orindate of the right bank station
#' @field ds_length_LOB downstream length of the left channel overbank
#' @field ds_length_Main downstream length of the main channel (to the next cross-section)
#' @field ds_length_ROB downstream length of the right channel overbank
#' @field contraction_coeff contraction coefficient for use in energy loss calculation
#' @field expansion_coeff expansion coefficient for use in energy loss calculation
xsection <- setRefClass("xsection",
                        field=list(riverstation="character",station="numeric",xx="numeric",zz="numeric",Manning="numeric",
                                   Manning_LOB="numeric",Manning_Main="numeric",Manning_ROB="numeric",
                                   lbs_xx="numeric",rbs_xx="numeric",ds_length_LOB="numeric",ds_length_Main="numeric",
                                   ds_length_ROB="numeric",contraction_coeff="numeric",expansion_coeff="numeric"),
                        method = list(initialize =
                                        function(..., riverstation=as.character(station),station=NA,xx=c(0),zz=c(0),Manning=c(0),
                                                 Manning_LOB=0.08,Manning_Main=0.035,Manning_ROB=0.08,
                                                 lbs_xx=0,rbs_xx=0,ds_length_LOB=0,ds_length_Main=0,
                                                 ds_length_ROB=0, contraction_coeff=0.1,expansion_coeff=0.3)
                                        {
                                          callSuper(..., riverstation = riverstation, station = station, xx = xx,zz=zz,Manning=Manning,
                                                    Manning_LOB=Manning_LOB,Manning_Main=Manning_Main,Manning_ROB=Manning_ROB,
                                                    lbs_xx=lbs_xx,rbs_xx=rbs_xx,ds_length_LOB=ds_length_LOB,ds_length_Main=ds_length_Main,
                                                    ds_length_ROB=ds_length_ROB,contraction_coeff=contraction_coeff,expansion_coeff=expansion_coeff)
                                        })
)

#' Geometry Object
#'
#' geom:
#' A reference class to create a geometry object from a collection of cross-sections of class xsection
#'
#' @field geomname name of the geometry object
#' @field xsectionList a list of xsection objects, which contains all of the cross-sections and relevant data
geom <- setRefClass("geom", field = list(geomname = "character", xsectionList = "list"))
validgeom <- function(object) {
  if (length(object$xsectionList) == 0) {
    stop("`geom` must contain at least one object of class `xsection`.")
  }
  sapply(object$xsectionList, function(x) {
    if (class(x)[1] != "xsection") {
      stop("Each element in `xsectionList` must be of class `xsection`")
    }
  })
  T
}
setValidity("geom", validgeom)



# want methods to:
# - check that all stations are different (no 2 stations the exact same)
# - add new xsections to geom
# - remove xsections from list by station/ name
# - sort from upstream - downstream or vice versa by station - # https://stackoverflow.com/questions/3309188/how-to-sort-a-listt-by-a-property-in-the-object
# - plot out basic network plot of sections along station line
# summary stats on geometry?

#
# ## Flows ----
#
# flow <- setRefClass("flow",
#                     field=list(riverstation="character",station="numeric",profile="character",reach="character",flow_cms="numeric"),
#                     method = list(initialize =
#                                     function(..., riverstation=as.character(station),station=NA,profile="",reach="",flow_cms=0)
#                                     {
#                                       callSuper(..., riverstation = riverstation, station = station, profile=profile, reach=reach, flow_cms=flow_cms)
#                                     })
# )
#
#
# flowprofile <- setRefClass("flowprofile", contains="flow",
#                            field=list(profilename="character",profiles=flow),
#                            method = list(initialize =
#                                            function(..., profilename="",profiles=NA)
#                                            {
#                                              callSuper(..., profilename=profilename,profiles=profiles)
#                                            })
# )


## Boundary Conditions ----

#' Boundary Condition Object
#'
#' bc:
#' A reference class to create a boundary conditions at a particular location in the model.
#' Note that currently only subcritical mode with normal depth boundary condition at the downstream-most
#' cross-section in the model is supported.
#'
#' @field riverstation name of the cross-section where the boundary condition is applied (character type)
#' @field station station location where the boundary condition is applied along profile (numeric type)
#' @field reach name of the reach
#' @field location description of the boundary condition location as either "Downstream", "Upstream"
#' @field bctype the type of boundary condition applied as "Normal Depth" (only supported type in current version)
#' @field bcvalue the value supplied with the boundary condition, currently the slope for normal depth computation
#' @field init_WSL the initial WSL value passed to boundary condition, used in normal depth calculation
bc <- setRefClass("bc",
                  field=list(riverstation="character",station="numeric",reach="character",location="character",
                             bctype="character",bcvalue="numeric",init_WSL="numeric"),
                  method = list(initialize =
                                  function(..., riverstation=as.character(station),station=NA,reach="",location="Downstream",
                                           bctype="Normal Depth",bcvalue=0.001,init_WSL=-99)
                                  {
                                    callSuper(..., riverstation = riverstation, station = station, reach=reach, location=location,
                                              bctype=bctype, bcvalue=bcvalue,init_WSL=init_WSL)
                                  })
)


## River Reaches ----



# ## Results object ----
#
# rcrresults <- setRefClass("rcrresults",
#                            field=list(profilename="character",flowprofiles="character",sections="character"),
#                            method = list(initialize =
#                                            function(..., profilename="",flowprofiles="",sections="")
#                                            {
#                                              callSuper(..., profilename=profilename,flowprofiles=flowprofiles,sections=sections)
#                                            })
# )

#
# ## Model Plan ----
#
# modeplan <- setRefClass("modelplan",contains=c("geom","flowprofile","bc","rcrresults"),
#                         field=list(name="character",desc="character",geom=geom,flowprofile=flowprofile,bc=bc,rcrresults=rcrresults),
#                         method = list(initialize =
#                                         function(..., name="Plan01",desc="",geom=NA,flowprofile=NA,bc=NA,rcrresults=NA)
#                                         {
#                                           callSuper(..., name=name,desc=desc,
#                                                     geom=geom,flowprofile=flowprofile,bc=bc,rcrresults=rcrresults)
#                                         })
# )



#' rcr Options Object
#'
#' rcr_options:
#' A reference class to define the computational parameters and physical constants used in the rcr package.
#' The parameters are defined with defaults when an instance of the options object is generated, however,
#' these may be modified as needed.
#'
#' @field g gravitational constant (N/kg)
#' @field dx resolution for interpolating cross-section data, defined by the horizontal space between interpolated points (m)
#' @field tolerance_cp tolerance limit for water surface calculations in compute_profile function (m)
#' @field iteration_limit_cp iteration limit for calculations between cross-sections in compute_profile
#' @field next_WSL_split_cp coefficient in assigning the next iteration of water surface elevations from the difference in previous WSL to calculated WSL
#' @field tolerance_nd tolerance limit for normal depth calculation in comparing flow to Manning's equation terms (m3/s)
#' @field iteration_limit_nd iteration limit for normal depth calculations
#' @field next_WSL_split_nd coefficient in assigning next water surface depth in normal depth calculation based on ratio of Manning's equation terms
#' @field silent_cp level of output in compute_profile function; if true, minimal output is produced
#' @field silent_nd level of output in normal_depth function; if true, minimal output is produced
rcr_options <- setRefClass("rcr_options",
                        field=list(g="numeric",dx="numeric",
                                   tolerance_cp="numeric",iteration_limit_cp="numeric",next_WSL_split_cp="numeric",
                                   tolerance_nd="numeric",iteration_limit_nd="numeric",next_WSL_split_nd="numeric",
                                   silent_cp="logical",silent_nd="logical"),
                        method = list(initialize =
                                        function(..., g=9.81,dx=0.001,
                                                 tolerance_cp=0.003,iteration_limit_cp=50,next_WSL_split_cp=0.7,
                                                 tolerance_nd=0.001,iteration_limit_nd=50,next_WSL_split_nd=0.4,
                                                 silent_cp=FALSE,silent_nd=FALSE)
                                        {
                                          callSuper(..., g=g,dx=dx,
                                                    tolerance_cp=tolerance_cp,iteration_limit_cp=iteration_limit_cp,next_WSL_split_cp=next_WSL_split_cp,
                                                    tolerance_nd=tolerance_nd,iteration_limit_nd=iteration_limit_nd,next_WSL_split_nd=next_WSL_split_nd,
                                                    silent_cp=silent_cp,silent_nd=silent_nd)
                                        })
)


