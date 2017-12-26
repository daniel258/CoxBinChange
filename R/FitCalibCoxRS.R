#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param w A matrix of time points when measurements on the binary covariate were obtained.
#' @param w.res A matrix of measurement results of the binary covariate. Each measurement corresponds to the time points in \code{w}
#' @param Q Matrix of covariates for PH calibration model
#' @param hz.times Times used for calculating the baseline hazard function from PH calibartion model
#' @param n.int The number of interior knots to be used, see \code{ICsurv::fast.PH.ICsurv.EM}, Default: 5
#' @param order the order of the basis functions. See \code{ICsurv::fast.PH.ICsurv.EM}, Default: 2
#' @param event Vector of censoring indicators. \code{1} for event \code{0} for censored
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
