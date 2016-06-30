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
CalcWeibullRiskSetP <- function(w, w.res, point, obs.tm)
{
  calc.lr <- CalcAuxAtPoint(w,w.res,point = point)
  df.lr <- calc.lr$df.lr
  df.lr[df.lr==Inf] <- NA
  df.lr[df.lr==0] <- NA
  df.lr.risk.set <- df.lr[obs.tm>=point,] # risksetca
  df.cln <- df.lr.risk.set[apply(df.lr.risk.set,1,function(x)  sum(is.na(x)))<2,]
  #  fit.weib <- fitdistcens(censdata = df.cln, distr = "weibull")
  fit.weib <-   tryCatch(fitdistcens(censdata = df.cln, distr = "weibull"),   error=function(e) {e})
   if (inherits(fit.weib, "error")) { 
      orc.weib.params <- FitCalibWeibull(w = w, w.res = w.res)
      p.point <- CalcWeibullCalibP(w = w, w.res = w.res, point = point, weib.params = orc.weib.params)
     } else if (fit.weib$estimate[1]> 20 | fit.weib$estimate[2] < 1/1000)
       {
       orc.weib.params <- FitCalibWeibull(w = w, w.res = w.res)
       p.point <- CalcWeibullCalibP(w = w, w.res = w.res, point = point, weib.params = orc.weib.params)
       } else {
              a.point <- calc.lr$a.point
              p.point <- calc.lr$x.one
              prob.at.point <- pweibull(point, shape = fit.weib$estimate[1],scale = fit.weib$estimate[2])
              prob.at.a.point <- pweibull(a.point[p.point==0], shape = fit.weib$estimate[1],scale = fit.weib$estimate[2])
              p.point[p.point==0] <- (prob.at.point - prob.at.a.point)/(1-prob.at.a.point)
              }
  return(p.point)
}
