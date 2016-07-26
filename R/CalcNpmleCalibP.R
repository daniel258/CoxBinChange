#### CalcNpmlePX function
### Daniel Nevo
## The function takes the following
## w - a matrix. Each row is observation and each column is questionnaire time in the interval. w equal to Inf once
# an observation is censore/had the event
## w.res - a matrix of the same dimensions as w. Equal to the x(t) at time w. For example second column is 
# second questionnaire result for all participents.
## point - scalar. The time of the risk set in the main analysis. In terms of the paper, t.
## fit.npmle - the result of icenReg::ic_np on the whole data. This used for the actual calibration
###
# The function returns a vector with individual predictions for P(X(t)=1|history(time t)). 
# For observations with X(a(t))=1 the above probability is 1 by definition and this is what the
# function returns for them.
#### The following package is needed: fitdistrplus
#### The following functions are used: CalcAuxatPoint (R function), CalcSurvFromNPMLE (cpp function)
CalcNpmleCalibP <- function(w, w.res, point, fit.npmle)
{
  lr.for.lik <- CalcAuxAtPoint(w,w.res,point = point)
  a.point <- lr.for.lik$a.point
  p.point <- lr.for.lik$x.one
  prob.at.point <- 1-CalcSurvFromNPMLE(probs = fit.npmle$p_hat, Tbull = fit.npmle$T_bull_Intervals,
                                       points = point)
  prob.at.a.point <- 1-CalcSurvFromNPMLE(probs = fit.npmle$p_hat, Tbull = fit.npmle$T_bull_Intervals,
                                         points = a.point[p.point==0])
  p.point[p.point==0] <- (prob.at.point - prob.at.a.point)/(1-prob.at.a.point)
  return(p.point)
}
