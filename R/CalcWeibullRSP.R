#### CalcWeibullRSP function
## The function takes the following
## w - a matrix. Each row is observation and each column is questionnaire time in the interval. w equal to Inf once
# an observation is censore/had the event
## w.res - a matrix of the same dimensions as w. Equal to the x(t) at time w. For example second column is 
# second questionnaire result for all participents.
## point - scalar. The time of the risk set in the main analysis. In terms of the paper, t.
## obs.tm - the vector of observed times for all observaions. In terms of the paper, T. This is used for finding the risk set.
###
# The function returns a vector with individual predictions for P(X(t)=1|history(time t)). 
# For observations with X(a(t))=1 the above probability is 1 by definition and this is what the
# function returns for them.
# The function calculate prediction for all observations, even though predictions for observations outside 
# the risk set are not used 
#### The following functions is used: CalcAuxatPoint (R function), FindIntervalCalibCPP (cpp function)
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param w A matrix of time points when measurements on the binary covariate were obtained.
#' @param w.res A matrix of measurement results of the binary covariate. Each measurement corresponds to the time points in \code{w}
#' @param point PARAM_DESCRIPTION
#' @param weib.params A bivariate vector. Shape and scale paramters of the Weibull calibration model.  
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[stats]{pweibull}}
#' @rdname CalcWeibullRSP
#' @export 
#' @importFrom stats pweibull
CalcWeibullRSP <- function(w, w.res, point, weib.params)
{
  weib.shape <- weib.params[1]
  weib.scale <- weib.params[2]
  lr.for.lik <- CalcAuxAtPoint(w,w.res,point = point)
  a.point <- lr.for.lik$a.point
  p.point <- lr.for.lik$x.one
  prob.at.point <- stats::pweibull(point, shape = weib.shape, scale = weib.scale)
  prob.at.a.point <- stats::pweibull(a.point[p.point==0], shape = weib.shape, scale = weib.scale)
  p.point[p.point==0] <- (prob.at.point - prob.at.a.point)/(1-prob.at.a.point)
  return(p.point)
}
