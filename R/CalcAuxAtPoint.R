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
# x.one : one/zero. Equals to 1 iff X(point)=1
# a.point: time of last questionnire, i.e., a(t) in terms of the paper
CalcAuxAtPoint <- function(w, w.res, point) {
n.sample <- nrow(w)
interval.w <-   FindIntervalCPP(point = point, w =w)

right.for.surv <- left.for.surv <- x.one <- a.point <- vector(length = n.sample)

for (j in 1:n.sample)
{
  if (interval.w[j]==1)
  {
    right.for.surv[j] <- Inf
    left.for.surv[j] <- 0
    a.point[j] <- 0
  } else  if (any(w.res[j,1:(interval.w[j]-1)]==1))
  {
    right.for.surv[j] <- w[j, Position(f = function(x) x==1,x=w.res[j,])]
    left.for.surv[j] <- ifelse(right.for.surv[j]==w[j,1], 0, 
                             w[j, Position(f = function(x) x==1,x=w.res[j,])-1])
    x.one[j] <- 1
    a.point[j] <- Inf 
  } else
    # if (any(w[j,1:interval.w[j]]==Inf))
  {
    right.for.surv[j] <- Inf
    left.for.surv[j] <- ifelse(w[j,1]==Inf, 0, 
                             w[j, interval.w[j]-1])
    a.point[j] <- left.for.surv[j]
                              
  }} 
#   else
#   {
#     right.for.surv[j] <- w[j,ncol(w)]
#     left.for.surv[j] <- w[j,ncol(w)-1]
#   }
# return(cbind(left.for.surv,right.for.surv))
ret.list <- list(df.lr = data.frame(left = left.for.surv, right = right.for.surv), x.one = x.one, a.point= a.point)
 return(ret.list)
}
