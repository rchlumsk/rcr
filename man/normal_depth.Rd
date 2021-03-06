% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normal_depth.R
\name{normal_depth}
\alias{normal_depth}
\title{Normal Depth Calculation}
\usage{
normal_depth(xs, Q, S, init_WSL = NA, options)
}
\arguments{
\item{xs}{cross-section of object class xsection}

\item{Q}{flow value through cross-section (m3/s)}

\item{S}{bed slope in cross-section, typically supplied as part of a boundary condition}

\item{init_WSL}{initial water surface level provided to initialize calculation}

\item{options}{object of class rcr_options with options and constants for hydraulic calculations}
}
\value{
\item{WSL}{Water surface elevation}
}
\description{
Calculates the normal depth for a given cross-section using the Manning's equation.
}
\examples{

}
\seealso{
\code{\link{bankfull_flow}} to estimate the flow through the cross section using Manning's equation

\code{\link{flow_area}} to calculate the flow area in cross-section for given segments

The Manning's n used in the calculation of the normal depth is calculated as a composite Manning's n,
weighted using the wetted pereimeter of the left overbank, main channel, and right overbank.
}
\keyword{depth}
\keyword{normal}
