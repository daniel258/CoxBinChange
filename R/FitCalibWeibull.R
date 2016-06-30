FitCalibWeibull <- function(w,w.res)
{
lr.for.fit <- as.data.frame(FindIntervalCalibCPP(w = w, wres = w.res))
colnames(lr.for.fit) <- c("left","right")
lr.for.fit[lr.for.fit==Inf] <- 200
lr.for.fit[lr.for.fit==0] <- 0.0001
#df.cln <- lr.for.fit[apply(lr.for.fit,1,function(x)  sum(is.na(x)))<2,]
fit.weib <- tryCatch(fitdistcens(censdata = lr.for.fit, distr = "weibull"),   error=function(e) {e})
if (inherits(fit.weib, "error")) { 
  return(c(NA,NA)) 
} else if (fit.weib$estimate[1]> 20 | fit.weib$estimate[2] < 1/1000)
{
  return(c(NA,NA)) 
}
return(fit.weib$estimate)
}