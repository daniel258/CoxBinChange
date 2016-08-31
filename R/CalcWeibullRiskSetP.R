#### CalcWeibullRiskSetP function
### Daniel Nevo
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

CalcWeibullRiskSetP <- function(w, w.res, point, obs.tm)
{
  lr.for.fit <- as.data.frame(FindIntervalCalibCPP(w = w, wres = w.res))
  colnames(lr.for.fit) <- c("left","right")
  lr.for.fit[lr.for.fit==Inf] <- 20
  lr.for.fit[lr.for.fit==0] <- 0.0001
  lr.for.fit <- lr.for.fit[obs.tm>point,]  # Keep only observations in the risk set
  fit.weib <- tryCatch(fitdistcens(censdata = lr.for.fit, distr = "weibull"),   error=function(e) {e})
  if (inherits(fit.weib, "error")) { 
    weib.params <- FitCalibWeibull(w,w.res)
    warning(paste("In point", point, "Calibration was used instead of risk set calibration"))
  } else if (fit.weib$estimate[1]> 20 | fit.weib$estimate[2] < 1/1000)
  {
    weib.params <- FitCalibWeibull(w,w.res)
    warning(paste("In point", point, "Calibration was used instead of risk set calibration"))
  } else   {    weib.params <- fit.weib$estimate}
  weib.shape <- weib.params[1]
  weib.scale <- weib.params[2]
  lr.for.lik <- CalcAuxAtPoint(w,w.res,point = point)
  a.point <- lr.for.lik$a.point
  p.point <- lr.for.lik$x.one
  prob.at.point <- pweibull(point, shape = weib.shape, scale = weib.scale)
  prob.at.a.point <- pweibull(a.point[p.point==0], shape = weib.shape, scale = weib.scale)
  p.point[p.point==0] <- (prob.at.point - prob.at.a.point)/(1-prob.at.a.point)
  return(p.point)
}


