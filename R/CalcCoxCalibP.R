#' @title Calculating the probabilities of positive binary exposure status at a given time point using a proportional hazards calibration model 
#' @description For a given time point, calculate the probability of positive exposure value  for multiple observations (participants). 
#' The function uses the results of a proportional hazards calibration model fit, and given covariates and collected data on the history 
#' of the binary exposure for each participant. 
#' @param w A matrix of time points when measurements on the binary covariate were obtained.
#' @param w.res A matrix of measurement results of the binary covariate. Each measurement corresponds to the time points in \code{w}
#' @param point The time point at which the probabilities are estimated
#' @param fit.cox The result of \code{icenReg::ic_sp} on the interval-censored data
#' @param hz.times Times used for calculating the baseline hazard function from PH calibartion model
#' @param Q Matrix of covariates for PH calibration model
#' @return A vector of estimated probabilities of positive exposure status at time \code{point}.
# @details DETAILS
# @examples 
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
#' @rdname CalcCoxCalibP
#' @export 
CalcCoxCalibP <- function(w, w.res, point, fit.cox, hz.times,  Q)
{
  lr.for.lik <- CalcAuxAtPoint(w,w.res,point = point)
  a.point <- lr.for.lik$a.point
  p.point <- lr.for.lik$x.one
  hz <- fit.cox$hz
  Qb <- Q%*%fit.cox$b
  ## Calculate hazard for the point, first baseline hazard, then add covariates:
  interval.point <- FindIntervalCPP(point = point, w = t(as.matrix(hz.times)))
  if (interval.point==1) {base.hz.point <- hz[1]*point/hz.times[1] } else {
    if(interval.point==length(hz.times)+1) {base.hz.point <- hz[length(hz.times)]
    } else {
      base.hz.point <- hz[interval.point-1] +  (hz[interval.point]-hz[interval.point-1])*(point-hz.times[interval.point-1])/
                    (hz.times[interval.point]-hz.times[interval.point-1]) #Extrapolation
    }}
  prob.at.point <- 1-exp(-base.hz.point*exp(Qb[p.point==0,]))
  prob.at.a.point <- 1-CalcSurvFromCox(fit.cox = fit.cox,Qb = Qb[p.point==0,], points = a.point[p.point==0], hz.times = hz.times)
  p.point[p.point==0] <- (prob.at.point - prob.at.a.point)/(1-prob.at.a.point)
  return(p.point)
}
