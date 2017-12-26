
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param fit.cox PARAM_DESCRIPTION
#' @param Qb PARAM_DESCRIPTION
#' @param points PARAM_DESCRIPTION
#' @param hz.times PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
# @rdname CalcSurvFromCox
CalcSurvFromCox <- function(fit.cox, Qb, points, hz.times)
{
  n.points <- length(points)
  hz <- fit.cox$hz
  surv.probs <- base<- vector(length = n.points)
  intervals <- sapply(points,FindIntervalCPP, w = t(as.matrix(hz.times)))
  for (i in 1:n.points)
  {
    if (intervals[i]==1)
    {
      base[i] <- hz[1]*points[i]/hz.times[1]
    } else
      if(intervals[i]==length(hz.times)+1) {base[i] <- hz[length(hz.times)]
      } else {
    base[i] <- (hz[intervals[i]-1] +  (hz[intervals[i]]-hz[intervals[i]-1])*(points[i]-hz.times[intervals[i]-1])/
                                                                        (hz.times[intervals[i]]-hz.times[intervals[i]-1]))#Extrapolation
      }
  }
  surv.probs <- exp(-base*exp(Qb))
}
