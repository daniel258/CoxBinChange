#### FitCalibCoxRSInts function
# Unlike FitCalibCoxRS, this function only estimate the survival in certain number of points
### Daniel Nevo
## The function takes the following
## w - a matrix. Each row is observation and each column is questionnaire time in the interval. w equal to Inf once
# an observation is censore/had the event
## w.res - a matrix of the same dimensions as w. Equal to the x(t) at time w. For example second column is 
# second questionnaire result for all participents.
## pts.for.ints: the points defining the intervals (first one has to be zero)  - should be sorted from zero up
###
# The function returns a list of cox fits for interval-censored time to event data, one for each risk set.
FitCalibCoxRSInts<- function(w, w.res, Z, hz.times, n.int = 5, order = 2 , tm, event, pts.for.ints)
{
if (pts.for.ints[1] != 0) {pts.for.ints <- c(0, pts.for.ints)}
r <- length(pts.for.ints)
event.index <- which(event==1)
lr.for.fit.all <- as.data.frame(FindIntervalCalibCPP(w = w, wres = w.res))
Z.all <- Z
all.fit.cox.res <- list()
for (j in 1:r)
{
  point <- pts.for.ints[j]
  # Keep only observations in the risk set
  lr.for.fit <- lr.for.fit.all[tm>=point,]  
  Z <- Z.all[tm>=point,]
  # Take out noninformative observations
  Z <- Z[!(lr.for.fit[,1]==0 & lr.for.fit[,2]==Inf),]
  lr.for.fit <- lr.for.fit[!(lr.for.fit[,1]==0 & lr.for.fit[,2]==Inf),]
  #
  colnames(lr.for.fit) <- c("left","right")
  d1 <- lr.for.fit[,1]==0
  d3 <- lr.for.fit[,2]==Inf
  d2 <- 1 - d1 - d3
  fit.cox.point <- tryCatch(fast.PH.ICsurv.EM(d1 = d1, d2 = d2, d3 = d3,Li = lr.for.fit[,1],
                                        Ri = lr.for.fit[,2], n.int = n.int, order = order,  Xp = Z, g0 =rep(1,n.int + order), b0 = rep(0,ncol(Z)),
                                        t.seq = hz.times, tol = 0.001), error = function(e){e})
  # while(inherits(fit.cox.point, "error") & n.int >= 2) { 
  #   n.int <- n.int - 1
  #   fit.cox.point <- tryCatch(fast.PH.ICsurv.EM(d1 = d1, d2 = d2, d3 = d3,Li = lr.for.fit[,1],
  #                                         Ri = lr.for.fit[,2], n.int = n.int, order = order,  Xp = Z, g0 =rep(1,n.int + order), b0 = rep(0,ncol(Z)),
  #                                         t.seq = hz.times, tol = 0.001), error = function(e){e})
  # }
  if (inherits(fit.cox.point, "error")) {
    fit.cox.point <- FitCalibCox(w = w, w.res = w.res, Z = Z, hz.times = hz.times, n.int = n.int, order = order)
    warning(paste("In point", point, "Calibration was used instead of risk set calibration"))
  }   else {
    ti <- c(lr.for.fit[d1 == 0,1], lr.for.fit[d3 == 0,2])
    fit.cox.point$knots <-   seq(min(ti) - 1e-05,  max(ti) + 1e-05, length.out = (n.int + 2))
    fit.cox.point$order <- order
  }
  all.fit.cox.res[[j]] <- fit.cox.point
}
return(all.fit.cox.res)
}
