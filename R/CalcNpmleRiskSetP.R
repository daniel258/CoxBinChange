#### CalcNpmleRiskSetP function
### Daniel Nevo
## The function takes the following
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
#### The following package is needed: fitdistrplus
#### The following functions are used: CalcAuxatPoint (R function), FindIntervalCalibCPP, CalcSurvFromNPMLE (cpp functions)
CalcNpmleRiskSetP <- function(w, w.res, point, obs.tm)
{
  lr.for.fit <- as.data.frame(FindIntervalCalibCPP(w = w, wres = w.res))
  colnames(lr.for.fit) <- c("left","right")
  #lr.for.fit[lr.for.fit==Inf] <- 20
  #lr.for.fit[lr.for.fit==0] <- 0.001
  lr.for.fit <- lr.for.fit[obs.tm>point,] # Keep only observations in the risk set
  fit.npmple.rs <-  tryCatch(ic_np(cbind(left,right)~0,data = lr.for.fit),   error=function(e) {e})
  if (inherits(fit.npmple.rs, "error")) { 
    fit.npmle <- FitCalibNpmle(w = w, w.res = w.res)
    p.point <- CalcNpmleCalibP(w = w, w.res = w.res, point = point, fit.npmle = fit.npmle)
    warning(paste("In point", point, "Calibration was used instead of risk set calibration"))
    return(p.point)
    }
    lr.for.lik <- CalcAuxAtPoint(w,w.res,point = point)
    a.point <- lr.for.lik$a.point
  p.point <- lr.for.lik$x.one
  prob.at.point <- 1-CalcSurvFromNPMLE(probs = fit.npmple.rs$p_hat, Tbull = fit.npmple.rs$T_bull_Intervals,
                                       points = point)
  prob.at.a.point <- 1-CalcSurvFromNPMLE(probs = fit.npmple.rs$p_hat, Tbull = fit.npmple.rs$T_bull_Intervals,
                                         points = a.point[p.point==0])
  # If the support of the estimated distriubtion ends before someones last available questionnire, the we get a contradiction
  # because \hat{F}(a)=1 but X(a)=0. In this rare case, we just do carry forward the zero.
  prob.at.a.point[prob.at.a.point==1] <-  -prob.at.point
  # If the last questionnire result was X(a)=1 then no problem caused above since p.point==1
  p.point[p.point==0] <- (prob.at.point - prob.at.a.point)/(1-prob.at.a.point)
  return(p.point)
}
