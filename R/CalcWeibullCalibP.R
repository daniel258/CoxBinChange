#### CalcWeibullPX function
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
# The function calculate prediction for all observations, even though predictions for observations outside 
# the risk set are not used
CalcWeibullCalibP <- function(w, w.res, point, weib.params)
{
  lr.for.lik <- CalcAuxAtPoint(w,w.res,point = point)
  a.point <- lr.for.lik$a.point
  p.point <- lr.for.lik$x.one
  prob.at.point <- pweibull(point, shape = weib.params[1],scale = weib.params[2])
  prob.at.a.point <- pweibull(a.point[p.point==0], shape = weib.params[1],scale = weib.params[2])
  p.point[p.point==0] <- (prob.at.point - prob.at.a.point)/(1-prob.at.a.point)
  return(p.point)
}
