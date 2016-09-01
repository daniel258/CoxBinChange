#### Functions for fitting ordinery calibration ###
# On Sep 1, 2016, this file included three calibrations: Weibull, Nonparameteric and Cox. 
# The Weibull calibrations return the shape and scale paramters
# The Cox and NP calibrations returns fitted objects.
# Sep 1, 2016: Only the Cox calibration function allows for covariates
################################################################################################################
################### Weibull ###############################################################################################
FitCalibWeibull <- function(w,w.res)
{
  lr.for.fit <- as.data.frame(FindIntervalCalibCPP(w = w, wres = w.res))
  lr.for.fit <- lr.for.fit[!(lr.for.fit[,1]==0 & lr.for.fit[,2]==Inf),]
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

################################################################################################################
################### Nonparametric (NPMLE)#####################################################################################
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
################################################################################################################
################### Cox ####################################################################################################
FitCalibCox <- function(w, w.res, Z, hz.times)
{
  lr.for.fit <- as.data.frame(FindIntervalCalibCPP(w = w, wres = w.res))
  Z <- Z[!(lr.for.fit[,1]==0 & lr.for.fit[,2]==Inf),]
  lr.for.fit <- lr.for.fit[!(lr.for.fit[,1]==0 & lr.for.fit[,2]==Inf),]
  lr.for.fit[lr.for.fit==Inf] <- 200
  lr.for.fit[lr.for.fit==0] <- 0.0001
  colnames(lr.for.fit) <- c("left","right")
  n.s <- nrow(lr.for.fit)
  fit.cox <- fast.PH.ICsurv.EM(d1 = rep(0,n.s), d2 = rep(1,n.s), d3 = rep(0,n.s),Li = lr.for.fit[,1],
                               Ri = lr.for.fit[,2], n.int = 10, order = 2,  Xp = Z, g0 =rep(1,12), b0 = c(0,2),
                               t.seq = hz.times, tol = 0.001)
  if (inherits(fit.cox, "error")) { 
    return(NA) 
  } else {
    return(fit.cox)
  }
}
######################################################