## point - scalar. The time of the risk set in the main analysis. In terms of the paper, t.
## weib.params - the shape and scale parameters from the Weibull calibration fitting to the interval-cenosed time to exposure/treatment
###
# The function returns a vector with individual predictions for P(X(t)=1|history(time t)). 
# For observations with X(a(t))=1 the above probability is 1 by definition and this is what the
# function returns for them.
# The function calculates prediction for all observations, even though predictions for observations outside 
# the risk set are not used
#### The following functions is used: CalcAuxatPoint (R function)
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
#' @rdname CalcWeibullCalibP
#' @export 
#' @importFrom stats pweibull
CalcWeibullCalibP <- function(w, w.res, point, weib.params)
{
  weib.shape <- weib.params[1]
  weib.scale <- weib.params[2]
  lr.for.lik <- CalcAuxAtPoint(w,w.res,point = point)
  a.point <- lr.for.lik$a.point
  p.point <- lr.for.lik$x.one
  prob.at.point <- stats::pweibull(point, shape = weib.shape,scale = weib.scale)
  prob.at.a.point <- stats::pweibull(a.point[p.point==0], shape = weib.shape, scale = weib.scale)
  p.point[p.point==0] <- (prob.at.point - prob.at.a.point)/(1-prob.at.a.point)
  return(p.point)
}
