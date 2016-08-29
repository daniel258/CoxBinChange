#### FitCalibNpmle function
### Daniel Nevo
## The function takes the following
## w - a matrix. Each row is observation and each column is questionnaire time in the interval. w equal to Inf once
# an observation is censore/had the event
## w.res - a matrix of the same dimensions as w. Equal to the x(t) at time w. For example second column is 
# second questionnaire result for all participents.
###
# The function returns the fit of NPMLE for interval-censored time to event data.
FitCalibNpmle <- function(w,w.res)
{
lr.for.fit <- as.data.frame(FindIntervalCalibCPP(w = w, wres = w.res))
lr.for.fit <- lr.for.fit[!(lr.for.fit[,1]==0 & lr.for.fit[,2]==Inf),]
colnames(lr.for.fit) <- c("left","right")
lr.for.fit[lr.for.fit==Inf] <- 200
lr.for.fit[lr.for.fit==0] <- 0.0001
#df.cln <- lr.for.fit[apply(lr.for.fit,1,function(x)  sum(is.na(x)))<2,]
fit.npmple <- ic_np(cbind(left,right)~0,data = lr.for.fit)
if (inherits(fit.npmple, "error")) { 
  return(NA) 
} else {
return(fit.npmple)
}
}

