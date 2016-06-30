#### CalcNpmlePX function
## June 08 2016
### Daniel Nevo
## The function takes the following
## w - a matrix. Each row is observation and each column is questioniire time in the interval. w equal to Inf once
# an observation is censore/had the main event
## w.res - a matrix of the same dimensions as w. Equal to the x(t) at time w. For example second column is 
# second questioniire result for all participents.
## point - scalar. The time of the risk set in the main analysis. In terms of the paper, t.
###
# The function returnsa vector with indeividual predictions for P(X(t)=1|history(time t)). 
# For observations with X(a(t))=1 the above probability is 1 by definition and this is what the
# function returns for them.
#### The following package is needed
# fitdistrplus
#### The following function is used
# CalcAuxatPoint 
CalcNpmleRiskSetP <- function(w, w.res, point, obs.tm)
{
  calc.lr <- CalcAuxAtPoint(w,w.res,point = point)
  df.lr <- calc.lr$df.lr
  df.lr.risk.set <- df.lr[obs.tm>=point,] # risksetca
  fit.npmple <- ic_np(cbind(left,right)~0,data = df.lr.risk.set)
  a.point <- calc.lr$a.point
  p.point <- calc.lr$x.one
  prob.at.point <- 1-CalcSurvFromNPMLE(probs = fit.npmple$p_hat, Tbull = fit.npmple$T_bull_Intervals,
                                       points = point)
  prob.at.a.point <- 1-CalcSurvFromNPMLE(probs = fit.npmple$p_hat, Tbull = fit.npmple$T_bull_Intervals,
                                         points = a.point[p.point==0])
  p.point[p.point==0] <- (prob.at.point - prob.at.a.point)/(1-prob.at.a.point)
  return(p.point)
}
