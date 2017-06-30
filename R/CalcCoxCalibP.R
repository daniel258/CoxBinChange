#### CalcCoxCalibP function
### Daniel Nevo
## The function takes the following
## w - a matrix. Each row is observation and each column is questionnaire time in the interval. w equal to Inf once
# an observation is censore/had the event
## w.res - a matrix of the same dimensions as w. Equal to the x(t) at time w. For example second column is 
# second questionnaire result for all participents.
## point - scalar. The time of the risk set in the main analysis. In terms of the paper, t.
## fit.cox - the result of icenReg::ic_sp on the whole data. This used for the actual calibration
###
# The function returns a vector with individual predictions for P(X(t)=1|history(time t)). 
# For observations with X(a(t))=1 the above probability is 1 by definition and this is what the
# function returns for them.
#### The following package is needed: fitdistrplus
#### The following functions are used: CalcAuxatPoint (R function)
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
                    (hz.times[interval.point]-hz.times[interval.point-1])#Extrapolation
    }}
  #base.hz.point <- fit.cox$hz[case.times==point]
  prob.at.point <- 1-exp(-base.hz.point*exp(Qb[p.point==0,]))
  prob.at.a.point <- 1-CalcSurvFromCox(fit.cox = fit.cox,Qb = Qb[p.point==0,], points = a.point[p.point==0], hz.times = hz.times)
 #   probs = fit.cox$p_hat, Tbull = fit.cox$T_bull_Intervals,
  #                                       points = a.point[p.point==0])
  p.point[p.point==0] <- (prob.at.point - prob.at.a.point)/(1-prob.at.a.point)
  return(p.point)
}
