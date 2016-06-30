#### CarryFow function
## June 08 2016
### Daniel Nevo
## The function takes the following
## w - a matrix. Each row is observation and each column is questioniire time in the interval. w equal to Inf once
# an observation is censore/had the main event
## w.res - a matrix of the same dimensions as w. Equal to the x(t) at time w. For example second column is 
# second questioniire result for all participents.
## point - scalar. The time of the risk set in the main analysis. In terms of the paper, t.
###
# The function returns a vector
# x.carry - vector of x(a(t))
CarryForward <- function(w, w.res, point) {
  n.sample <- nrow(w)
  interval.w <-   FindIntervalCPP(point = point, w =w)
#  interval.w[is.na(interval.w)] <- ncol(w) + 1 
   x.carry <- vector(length = n.sample)
   for (j in 1:n.sample)
   {
     if (interval.w[j]==1)
     {
       x.carry[j] <- 0
     } else  if (any(w.res[j,1:(interval.w[j]-1)]==1))
     {
       x.carry[j] <- 1
    } else
     {
       x.carry[j] <- 0
     }} 
   #   else
   #   {
   #     right.for.surv[j] <- w[j,ncol(w)]
   #     left.for.surv[j] <- w[j,ncol(w)-1]
   #   }
   # return(cbind(left.for.surv,right.for.surv))
   return(x.carry)
  return(interval.w)
}
