#### FitCalibWeibullRS function
### Daniel Nevo
## The function takes the following
## w - a matrix. Each row is observation and each column is questionnaire time in the interval. w equal to Inf once
# an observation is censore/had the event
## w.res - a matrix of the same dimensions as w. Equal to the x(t) at time w. For example second column is 
# second questionnaire result for all participents.
###
# The function returns the estimates of Weibull scale and shape paramters for interval-censored time to event data.
#' @title FitCalibWeibullRS
#' @description Fit a parametric Weibull calibration model for time-to-exposure
#' @param w A matrix of measurement times: each row is a participant, each entry is meausrment time, see example below
#' @param w.res A matrix of exposure measurments corresponding to the matrix w
#' @param tm PARAM_DESCRIPTION
#' @param event PARAM_DESCRIPTION
#' @param lower PARAM_DESCRIPTION, Default: 1e-04
#' @param uppper PARAM_DESCRIPTION, Default: 200
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[fitdistrplus]{fitdistcens}}

#'  \code{\link[fitdistrplus]{fitdistcens}}
#' @rdname FitCalibWeibullRS
#' @export 
#' @importFrom fitdistrplus fitdistcens
#' @importFrom fitdistrplus fitdistcens
FitCalibWeibullRS <- function(w, w.res, tm, event, lower = 0.0001, upper = 200)
{
r <- sum(event)
event.index <- which(event==1)
lr.for.fit.all <- as.data.frame(FindIntervalCalibCPP(w = w, wres = w.res))
weib.params <- matrix(nr = r, nc = 2)
for (j in 1:r)
{
  point <- tm[event.index[j]]
  lr.for.fit <- lr.for.fit.all[tm>=point, ]  # Keep only observations in the risk set
  lr.for.fit <- lr.for.fit[!(lr.for.fit[, 1]==0 & lr.for.fit[, 2]==Inf),]
  colnames(lr.for.fit) <- c("left","right")
  lr.for.fit[lr.for.fit==Inf] <- upper
  lr.for.fit[lr.for.fit==0] <- lower
  fit.weib <- tryCatch(fitdistrplus::fitdistcens(censdata = lr.for.fit, distr = "weibull"),   error=function(e) {e})
  if (inherits(fit.weib, "error")) { 
    weib.params[j,] <- FitCalibWeibull(w, w.res)
    warning(paste("In point", point, "Calibration was used instead of risk set calibration"))
    } else if (fit.weib$estimate[1] > 20 | fit.weib$estimate[2] < 1/1000)
      {
      weib.params[j,] <- FitCalibWeibull(w, w.res)
      warning(paste("In point", point, "Calibration was used instead of risk set calibration"))
      } else   {    
        weib.params[j,] <- fit.weib$estimate
        }
}
return(weib.params)
}

