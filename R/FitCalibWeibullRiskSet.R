#### FitCalibWeibull function
### Daniel Nevo
## The function takes the following
## w - a matrix. Each row is observation and each column is questionnaire time in the interval. w equal to Inf once
# an observation is censore/had the event
## w.res - a matrix of the same dimensions as w. Equal to the x(t) at time w. For example second column is 
# second questionnaire result for all participents.
###
# The function returns the estimates of Weibull scale and shape paramters for interval-censored time to event data.
FitCalibWeibullRS <- function(w, w.res, obs.tm, event)
{
r <- sum(event)
event.index <- which(event==1)
lr.for.fit.all <- as.data.frame(FindIntervalCalibCPP(w = w, wres = w.res))
weib.params <- matrix(nr = r, nc = 2)
for (j in 1:r)
{
  point <- obs.tm[event.index[j]]
  lr.for.fit <- lr.for.fit.all[obs.tm>=point,]  # Keep only observations in the risk set
  lr.for.fit <- lr.for.fit[!(lr.for.fit[,1]==0 & lr.for.fit[,2]==Inf),]
  colnames(lr.for.fit) <- c("left","right")
  lr.for.fit[lr.for.fit==Inf] <- 200
  lr.for.fit[lr.for.fit==0] <- 0.0001
  fit.weib <- tryCatch(fitdistcens(censdata = lr.for.fit, distr = "weibull"),   error=function(e) {e})
  if (inherits(fit.weib, "error")) { 
    weib.params[j,] <- FitCalibWeibull(w,w.res)
    warning(paste("In point", point, "Calibration was used instead of risk set calibration"))
    } else if (fit.weib$estimate[1]> 20 | fit.weib$estimate[2] < 1/1000)
      {
      weib.params[j,] <- FitCalibWeibull(w,w.res)
      warning(paste("In point", point, "Calibration was used instead of risk set calibration"))
      } else   {    
        weib.params[j,] <- fit.weib$estimate
        }
}
return(weib.params)
}

