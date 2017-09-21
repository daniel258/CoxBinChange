#### FitCalibCoxRS function
### Daniel Nevo
## The function takes the following
## w - a matrix. Each row is observation and each column is questionnaire time in the interval. w equal to Inf once
# an observation is censore/had the event
## w.res - a matrix of the same dimensions as w. Equal to the x(t) at time w. For example second column is 
# second questionnaire result for all participents.
###
# The function returns a list of cox fits for interval-censored time to event data, one for each risk set.
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param w PARAM_DESCRIPTION
#' @param w.res PARAM_DESCRIPTION
#' @param Q PARAM_DESCRIPTION
#' @param hz.times PARAM_DESCRIPTION
#' @param n.int PARAM_DESCRIPTION, Default: 5
#' @param order PARAM_DESCRIPTION, Default: 2
#' @param event PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[ICsurv]{fast.PH.ICsurv.EM}}
#' @rdname FitCalibCoxRS
#' @export 
#' @importFrom ICsurv fast.PH.ICsurv.EM
FitCalibCoxRS <- function(w, w.res, Q, hz.times, n.int = 5, order = 2 , event)
{
r <- sum(event)
event.index <- which(event==1)
lr.for.fit.all <- as.data.frame(FindIntervalCalibCPP(w = w, wres = w.res))
Q.all <- Q
all.fit.cox.res <- list()
for (j in 1:r)
{
  point <- obs.tm[event.index[j]]
  # Keep only observations in the risk set
  lr.for.fit <- lr.for.fit.all[obs.tm>=point,]  
  Q <- Q.all[obs.tm>=point,]
  # Take out noninformative observations
  Q <- Q[!(lr.for.fit[,1]==0 & lr.for.fit[,2]==Inf),]
  lr.for.fit <- lr.for.fit[!(lr.for.fit[,1]==0 & lr.for.fit[,2]==Inf),]
  #
  colnames(lr.for.fit) <- c("left","right")
  d1 <- lr.for.fit[,1]==0
  d3 <- lr.for.fit[,2]==Inf
  d2 <- 1 - d1 - d3
  fit.cox.point <- tryCatch(ICsurv::fast.PH.ICsurv.EM(d1 = d1, d2 = d2, d3 = d3,Li = lr.for.fit[,1],
                                        Ri = lr.for.fit[,2], n.int = n.int, order = order,  Xp = Q, g0 =rep(1,n.int + order), b0 = rep(0,ncol(Q)),
                                        t.seq = hz.times, tol = 0.001), error = function(e){e})
  while(inherits(fit.cox.point, "error") & n.int >= 2) { 
    n.int <- n.int - 1
    fit.cox.point <- tryCatch(ICsurv::fast.PH.ICsurv.EM(d1 = d1, d2 = d2, d3 = d3,Li = lr.for.fit[,1],
                                          Ri = lr.for.fit[,2], n.int = n.int, order = order,  Xp = Q, g0 =rep(1,n.int + order), b0 = rep(0,ncol(Q)),
                                          t.seq = hz.times, tol = 0.001), error = function(e){e})
  }
  if (n.int<2) {
    fit.cox.point <- FitCalibCox(w = w, w.res = w.res, Q = Q, hz.times = hz.times, n.int = n.int, order = order)
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
