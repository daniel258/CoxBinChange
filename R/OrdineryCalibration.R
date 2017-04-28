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
  lr.for.fit[lr.for.fit==Inf] <- NA
#  lr.for.fit[lr.for.fit==0] <- 0.0001
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
  #lr.for.fit <- lr.for.fit[!(lr.for.fit[,1]==0 & lr.for.fit[,2]==Inf),]
  colnames(lr.for.fit) <- c("left","right")
  #lr.for.fit[lr.for.fit==Inf] <- 200
  #lr.for.fit[lr.for.fit==0] <- 0.0001
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
FitCalibCox <- function(w, w.res, Z, hz.times, n.int = 5, order = 2 )
{
  lr.for.fit <- as.data.frame(FindIntervalCalibCPP(w = w, wres = w.res))
  Z <- as.matrix(Z[!(lr.for.fit[,1]==0 & lr.for.fit[,2]==Inf),])
  lr.for.fit <- lr.for.fit[!(lr.for.fit[,1]==0 & lr.for.fit[,2]==Inf),]
  #lr.for.fit[lr.for.fit==Inf] <- 200
  #lr.for.fit[lr.for.fit==0] <- 0.0001
  colnames(lr.for.fit) <- c("left","right")
  n.s <- nrow(lr.for.fit)
  d1 <- lr.for.fit[,1]==0
  d3 <- lr.for.fit[,2]==Inf
  d2 <- 1 - d1 - d3
  fit.cox <- tryCatch(fast.PH.ICsurv.EM(d1 = d1, d2 = d2, d3 = d3,Li = lr.for.fit[,1],
                               Ri = lr.for.fit[,2], n.int = n.int, order = order,  Xp = Z, g0 =rep(1,n.int + order), b0 = rep(0,ncol(Z)),
                               t.seq = hz.times, tol = 0.001), error = function(e){e})
  while(inherits(fit.cox, "error") & n.int >= 2) { 
    n.int <- n.int - 1
    fit.cox <- tryCatch(fast.PH.ICsurv.EM(d1 = d1, d2 = d2, d3 = d3,Li = lr.for.fit[,1],
                                 Ri = lr.for.fit[,2], n.int = n.int, order = order,  Xp = Z, g0 =rep(1,n.int + order), b0 = rep(0,ncol(Z)),
                                 t.seq = hz.times, tol = 0.001), error = function(e){e})
  }
  if (n.int<2) {return(NA)}   else {
    ti <- c(lr.for.fit[d1 == 0,1], lr.for.fit[d3 == 0,2])
    fit.cox$knots <-   seq(min(ti) - 1e-05,  max(ti) + 1e-05, length.out = (n.int + 2))
    fit.cox$order <- order
    return(fit.cox)
  }
}
######################################################