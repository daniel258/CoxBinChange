#### FitCalibCox function
### Daniel Nevo
## The function takes the following
## w - a matrix. Each row is observation and each column is questionnaire time in the interval. w equal to Inf once
# an observation is censore/had the event
## w.res - a matrix of the same dimensions as w. Equal to the x(t) at time w. For example second column is 
# second questionnaire result for all participents.
###
# The function returns the fit of NPMLE for interval-censored time to event data.
FitCalibCox <- function(w,w.res,df.vars)
{
lr.for.fit <- as.data.frame(FindIntervalCalibCPP(w = w, wres = w.res))
df.vars <- df.vars[!(lr.for.fit[,1]==0 & lr.for.fit[,2]==Inf),]
lr.for.fit <- lr.for.fit[!(lr.for.fit[,1]==0 & lr.for.fit[,2]==Inf),]
lr.for.fit[lr.for.fit==Inf] <- 200
lr.for.fit[lr.for.fit==0] <- 0.0001
colnames(lr.for.fit) <- c("left","right")
n.s <- nrow(lr.for.fit)
fit.cox <- fast.PH.ICsurv.EM(d1 = rep(0,n.s), d2 = rep(1,n.s), d3 = rep(0,n.s),Li = lr.for.fit[,1],
                       Ri = lr.for.fit[,2], n.int = 10, order = 2,  Xp = as.matrix(df.vars), g0 =rep(1,12), b0 = c(0,2),
                       t.seq = case.times, tol = 0.001)
if (inherits(fit.cox, "error")) { 
  return(NA) 
} else {
return(fit.cox)
}
}