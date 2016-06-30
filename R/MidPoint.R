#### CalcAuxAtPoint function
## June 08 2016
### Daniel Nevo
## The function takes the following
## w - a matrix. Each row is observation and each column is questioniire time in the interval. w equal to Inf once
# an observation is censore/had the main event
## w.res - a matrix of the same dimensions as w. Equal to the x(t) at time w. For example second column is 
# second questioniire result for all participents.
## point - scalar. The time of the risk set in the main analysis. In terms of the paper, t.
###
# The function returns a list with theree objects:
# df.lr: a data frame of the intervals which the evnet (aspirin takin is contianed)
# x.one : one/zero. Equals to iff X(point)=1
# a.point: time of ast questionnire, i.e., a(t) in terms of the paper
MidPoint <- function(w, w.res, point) 
  {
    lr.for.fit <- as.data.frame(FindIntervalCalibCPP(w = w, wres = w.res))
    colnames(lr.for.fit) <- c("left","right")
    tm.v <- vector(length=nrow(w)) #time of chhange
    tm.v[lr.for.fit$right==Inf] <- Inf
    tm.v[lr.for.fit$right<Inf] <- (lr.for.fit$right[lr.for.fit$right<Inf]+lr.for.fit$left[lr.for.fit$right<Inf])/2
    x <- ifelse(tm.v<point,1,0)
    return(x)
  }