#### CalcWeibullCalibP function
### Daniel Nevo
## The function takes the following
## w - a matrix. Each row is observation and each column is questionnaire time in the interval. w equal to Inf once
# an observation is censore/had the event
## w.res - a matrix of the same dimensions as w. Equal to the x(t) at time w. For example second column is 
# second questionnaire result for all participents.
## point - scalar. The time of the risk set in the main analysis. In terms of the paper, t.
## weib.params - the shape and scale parameters from the Weibull calibration fitting to the interval-cenosed time to exposure/treatment
###
# The function returns a vector with individual predictions for P(X(t)=1|history(time t)). 
# For observations with X(a(t))=1 the above probability is 1 by definition and this is what the
# function returns for them.
# The function calculates prediction for all observations, even though predictions for observations outside 
# the risk set are not used
#### The following functions is used: CalcAuxatPoint (R function)

CalcWeibullCalibPderivEta1 <- function(w, w.res, point, weib.params)
{
  lr.for.lik <- CalcAuxAtPoint(w,w.res,point = point)
  weib.shape <- weib.params[1]
  weib.scale <- weib.params[2]
  a.point <- lr.for.lik$a.point
  p.point <- lr.for.lik$x.one
  surv.at.point <- pweibull(point, shape = weib.shape,scale = weib.scale, lower.tail = F)
  surv.at.a.point <- pweibull(a.point[p.point==0], shape = weib.shape, scale = weib.scale, lower.tail = F)
  deriv.eta1 <- vector(length=length(p.point))
  deriv.eta1[p.point=0] <- 1-((point^weib.shape)*surv.at.point)/((a.point[p.point==0]^weib.shape)*surv.at.a.point)
  deriv.eta1[p.point>0] <- 0
  return(deriv.eta1)
}
CalcWeibullCalibPderivEta2 <- function(w, w.res, point, weib.params)
{
  weib.shape <- weib.params[1]
  weib.scale <- weib.params[2]
  lr.for.lik <- CalcAuxAtPoint(w,w.res,point = point)
  a.point <- lr.for.lik$a.point
  p.point <- lr.for.lik$x.one
  surv.at.point <- pweibull(point, shape = weib.shape,scale = weib.scale, lower.tail = F)
  surv.at.a.point <- pweibull(a.point[p.point==0], shape = weib.shape, scale = weib.scale, lower.tail = F)
  deriv.eta2 <- vector(length=length(p.point))
  deriv.eta2[p.point=0] <- 1-(log(point/weib.scale)*(point^weib.shape)*surv.at.point)/
    (log(a.point[p.point==0]/weib.scale)*(a.point[p.point==0]^weib.shape)*surv.at.a.point)
  deriv.eta2[p.point>0] <- 0
  return(deriv.eta2)
}