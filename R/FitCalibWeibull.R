#### FitCalibWeibull function
### Daniel Nevo
## The function takes the following
## w - a matrix. Each row is observation and each column is questionnaire time in the interval. w equal to Inf once
# an observation is censore/had the event
## w.res - a matrix of the same dimensions as w. Equal to the x(t) at time w. For example second column is 
# second questionnaire result for all participents.
###
# The function returns the estimates of Weibull scale and shape paramters for interval-censored time to event data.
FitCalibWeibull <- function(w,w.res)
{
lr.for.fit <- as.data.frame(FindIntervalCalibCPP(w = w, wres = w.res))
colnames(lr.for.fit) <- c("left","right")
lr.for.fit[lr.for.fit==Inf] <- 200
lr.for.fit[lr.for.fit==0] <- 0.0001
#df.cln <- lr.for.fit[apply(lr.for.fit,1,function(x)  sum(is.na(x)))<2,]
fit.weib <- tryCatch(fitdistcens(censdata = lr.for.fit, distr = "weibull"),   error=function(e) {e})
if (inherits(fit.weib, "error")) { 
  fit.weib <- tryCatch(fitdistcens(censdata = lr.for.fit, distr = "weibull", lower = c(0, 0)),   error=function(e) {e})
  if (inherits(fit.weib, "error")) { 
  return(c(NA,NA)) }
} 
if (fit.weib$estimate[1]> 20 | fit.weib$estimate[2] < 1/1000)
{
  return(c(NA,NA)) 
}
return(fit.weib$estimate)
}