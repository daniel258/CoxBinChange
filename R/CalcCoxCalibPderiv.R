#### CalcCoxCalibPderiv function
### NEED TO WRITE EXPLANATION
### Daniel Nevo
## The function takes the following
## w - a matrix. Each row is observation and each column is questionnaire time in the interval. w equal to Inf once
# an observation is censore/had the event
## w.res - a matrix of the same dimensions as w. Equal to the x(t) at time w. For example second column is 
# second questionnaire result for all participents.
## point - scalar. The time of the risk set in the main analysis. In terms of the paper, t.
## weib.params - the shape and scale parameters from the Weibull calibration fitting to the interval-cenosed time to exposure/treatment
###
# The function returns a vector with individual predictions for P(X(t)=1|history(time t)). 
# For observations with X(a(t))=1 the above probability is 1 by definition and this is what the
# function returns for them.
# The function calculates prediction for all observations, even though predictions for observations outside 
# the risk set are not used
#### The following function is used: CalcAuxatPoint (R function)
CalcCoxCalibPderiv <- function(w, w.res, point, fit.cox, hz.times, Z)
{
  eta.b <- fit.cox$b
  eta.g <- fit.cox$g
  knots <- fit.cox$knots
  order <- fit.cox$order
  lr.for.lik <- CalcAuxAtPoint(w,w.res,point = point)
  a.point <- lr.for.lik$a.point
  p.point <- lr.for.lik$x.one
  hz <- fit.cox$hz
  Zb <- Z%*%eta.b
  exp.Zb <- exp(Zb)
  ## Calculate hazard for the point, first baseline hazard, then add covariates:
  interval.point <- FindIntervalCPP(point = point, w = t(as.matrix(hz.times)))
  if (interval.point==1) {base.hz.point <- hz[1]*point/hz.times[1] } else {
    if(interval.point==length(hz.times)+1) {base.hz.point <- hz[length(hz.times)]
    } else {
      base.hz.point <- hz[interval.point-1] +  (hz[interval.point]-hz[interval.point-1])*(point-hz.times[interval.point-1])/
        (hz.times[interval.point]-hz.times[interval.point-1])#Extrapolation
    }}
  
  H.point <- base.hz.point*exp.Zb[p.point==0]
  S.at.point <- exp(-H.point)
  S.at.a.point <- CalcSurvFromCox(fit.cox = fit.cox,Zb = Zb[p.point==0,], points = a.point[p.point==0], hz.times = hz.times)
  H.a.point <- -log(S.at.a.point)
  
  b.point <- t(Ispline(x = point, order = order, knots = knots))
  b.a.point <- t(Ispline(x = a.point[p.point==0], order = order, knots = knots))
  
  deriv.eta <- matrix(nr = length(p.point), nc = length(eta.b) + length(eta.g),0) 
  deriv.eta[p.point==0,1:ncol(Z)] <- (S.at.point/S.at.a.point)*(H.point-H.a.point)*Z[p.point==0,]
  deriv.eta[p.point==0,(ncol(Z)+1):ncol(deriv.eta)] <- (S.at.point/S.at.a.point)*(as.vector(b.point)-b.a.point)*exp.Zb[p.point==0]
  return(deriv.eta)
}


CalcCoxCalibPderivRSInsts <- function(w, w.res, point, fit.cox.rs.ints, hz.times, Z,  pts.for.ints, tm)
{
  interval <- findInterval(point, pts.for.ints)
  fit.cox.int <- fit.cox.rs.ints[[interval]]
  eta.b <- fit.cox.int$b
  eta.g <- fit.cox.int$g
  knots <- fit.cox.int$knots
  order <- fit.cox.int$order
  in.risk.set <- tm >= point
  lr.for.lik <- CalcAuxAtPoint(w,w.res,point = point)
  a.point <- lr.for.lik$a.point
  p.point <- lr.for.lik$x.one
  hz <- fit.cox.int$hz
  Zb <- Z%*%eta.b
  exp.Zb <- exp(Zb)
  ## Calculate hazard for the point, first baseline hazard, then add covariates:
  interval.point <- FindIntervalCPP(point = point, w = t(as.matrix(hz.times)))
  if (interval.point==1) {base.hz.point <- hz[1]*point/hz.times[1] } else {
    if(interval.point==length(hz.times)+1) {base.hz.point <- hz[length(hz.times)]
    } else {
      base.hz.point <- hz[interval.point-1] +  (hz[interval.point]-hz[interval.point-1])*(point-hz.times[interval.point-1])/
        (hz.times[interval.point]-hz.times[interval.point-1])#Extrapolation
    }}
  
  H.point <- base.hz.point*exp.Zb[p.point==0 & in.risk.set]
  S.at.point <- exp(-H.point)
  S.at.a.point <- CalcSurvFromCox(fit.cox = fit.cox.int,Zb = Zb[p.point==0 & in.risk.set,], points = a.point[p.point==0 & in.risk.set], 
                                  hz.times = hz.times)
  H.a.point <- -log(S.at.a.point)
  
  b.point <- t(Ispline(x = point, order = order, knots = knots))
  b.a.point <- t(Ispline(x = a.point[p.point==0 & in.risk.set], order = order, knots = knots))
  
  deriv.eta.ints <- matrix(nr = length(p.point), nc = length(eta.b) + length(eta.g),0) 
  deriv.eta.ints[p.point==0 & in.risk.set, 1:ncol(Z)] <- (S.at.point/S.at.a.point)*(H.point-H.a.point)*Z[p.point==0 & in.risk.set,]
  deriv.eta.ints[p.point==0 & in.risk.set, (ncol(Z)+1):ncol(deriv.eta.ints)] <- (S.at.point/S.at.a.point)*(as.vector(b.point)-b.a.point)*
                                                                                exp.Zb[p.point==0  & in.risk.set]
  
  deriv.eta <- matrix(nr = length(p.point), nc = (length(eta.b) + length(eta.g))*length(fit.cox.rs.ints),0) 
  deriv.eta[,((interval-1)*(length(eta.b) + length(eta.g))+1):(interval*(length(eta.b) + length(eta.g)))] <- deriv.eta.ints
  
  return(deriv.eta)
}



# CalcCoxPderiv <- function(w, w.res, ps, fit.cox, Z)
# {
#   eta.b <- fit.cox$b
#   eta.g <- fit.cox$g
#   exp.Zb <- exp(Z%*%eta.b)
#   lr.for.lik <- CalcAuxAtPoint(w,w.res,point = point)
#   a.point <- lr.for.lik$a.point
#   p.point <- lr.for.lik$x.one
#   surv.at.point <- exp(-bpoint*exp.Zb)
#     pweibull(point, shape = weib.shape,scale = weib.scale, lower.tail = F)
#   surv.at.a.point <- pweibull(a.point[p.point==0], shape = weib.shape, scale = weib.scale, lower.tail = F)
#   deriv <- vector(length=length(p.point))
#   deriv.scale[p.point==0] <- -(weib.shape/(weib.scale^(weib.shape + 1))) * (point^weib.shape - a.point[p.point==0]^weib.shape) * (surv.at.point/surv.at.a.point)
#   deriv.scale[p.point>0] <- 0
#   #deriv.scale[a.point==0] <- 0
#   return(deriv.scale)
# }
# CalcWeibullCalibPderivShape <- function(w, w.res, point, weib.params)
# {
#   weib.shape <- weib.params[1]
#   weib.scale <- weib.params[2]
#   lr.for.lik <- CalcAuxAtPoint(w,w.res,point = point)
#   a.point <- lr.for.lik$a.point
#   p.point <- lr.for.lik$x.one
#   surv.at.point <- pweibull(point, shape = weib.shape,scale = weib.scale, lower.tail = F)
#   surv.at.a.point <- pweibull(a.point[p.point==0], shape = weib.shape, scale = weib.scale, lower.tail = F)
#   deriv.shape <- vector(length=length(p.point))
#   deriv.shape[p.point==0] <- (1/(weib.scale^weib.shape)) * (log(a.point[p.point==0]/weib.scale)*a.point[p.point==0]^weib.shape - 
#                                                              log(point/weib.scale)*point^weib.shape ) * (surv.at.point/surv.at.a.point)
#   deriv.shape[p.point>0] <- 0
#   deriv.shape[a.point==0] <- -(1/(weib.scale^weib.shape)) * log(point/weib.scale)*point^weib.shape * surv.at.point
#   return(deriv.shape)
# }
# 
# 
# # #########################################################################################################################
# # CalcWeibullCalibPderivShapeRS <- function(w, w.res, obs.tm, event, weib.rs.params)
# # {
# #   r <- sum(event)
# #   n <- length(event)
# #   event.index <- which(event==1)
# #   deriv.shape <- matrix(nr = n, nc = r)
# #   for (j in 1:r)
# #   {
# #   point <- obs.tm[event.index[j]]
# #   weib.rs.shape <- weib.rs.params[j, 1]
# #   weib.rs.scale <- weib.rs.params[j, 2]
# #   in.risk.set <- obs.tm>=point
# #   lr.for.lik <- CalcAuxAtPoint(w, w.res, point = point)
# #   a.point <- lr.for.lik$a.point
# #   p.point <- lr.for.lik$x.one
# #   surv.at.point <- pweibull(point, shape = weib.rs.shape, scale = weib.rs.scale, lower.tail = F)
# #   surv.at.a.point <- pweibull(a.point[p.point==0 & in.risk.set], shape = weib.rs.shape, scale = weib.rs.scale, lower.tail = F)
# #   deriv.shape[p.point==0 & in.risk.set, j] <- (1/(weib.rs.scale^weib.rs.shape)) * (log(a.point[p.point==0 & in.risk.set]/weib.rs.scale) * 
# #                                                 a.point[p.point==0 & in.risk.set]^weib.rs.shape - 
# #                                                   log(point/weib.rs.scale)*point^weib.rs.shape ) *
# #                                                                (surv.at.point/surv.at.a.point)
# #   deriv.shape[p.point>0, j] <- 0
# #   deriv.shape[a.point==0 & in.risk.set, j] <- -(1/(weib.rs.scale^weib.rs.shape)) * log(point/weib.rs.scale)*point^weib.rs.shape * surv.at.point
# #   deriv.shape[!in.risk.set, j] <- 0
# #   }
# #   return(t(deriv.shape))
# # }
# # 
# # #########################################################################################################################
# # 
# # 
# # CalcWeibullCalibPderivScaleRS <- function(w, w.res, obs.tm, event, weib.rs.params)
# # {
# #   r <- sum(event)
# #   n <- length(event)
# #   event.index <- which(event==1)
# #   deriv.scale <- matrix(nr = n, nc = r)
# #   for (j in 1:r)
# #   {
#     point <- obs.tm[event.index[j]]
#     weib.rs.shape <- weib.rs.params[j, 1]
#     weib.rs.scale <- weib.rs.params[j, 2]
#     in.risk.set <- obs.tm>=point
#     lr.for.lik <- CalcAuxAtPoint(w, w.res, point = point)
#     a.point <- lr.for.lik$a.point
#     p.point <- lr.for.lik$x.one
#     surv.at.point <- pweibull(point, shape = weib.rs.shape,scale = weib.rs.scale, lower.tail = F)
#     surv.at.a.point <- pweibull(a.point[p.point==0 & in.risk.set], shape = weib.rs.shape, scale = weib.rs.scale, lower.tail = F)
#     deriv.scale[p.point==0 & in.risk.set, j] <- -(weib.rs.shape/(weib.rs.scale^(weib.rs.shape + 1))) * 
#     (point^weib.rs.shape - a.point[p.point==0 & in.risk.set]^weib.rs.shape) * (surv.at.point/surv.at.a.point)
#     deriv.scale[p.point>0 , j] <- 0
#    # deriv.scale[a.point==0 , j] <- 0
#     deriv.scale[!in.risk.set, j] <- 0
#   #deriv.eta1[a.point==0] <- 0
#   }
#   return(t(deriv.scale))
# }