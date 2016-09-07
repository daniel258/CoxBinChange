###


CalcVarTheta <- function(beta, etas, tm, event, ps, ps.deriv.shape, ps.deriv.scale, w, w.res)
{
  n <- length(tm)
  b.vec <- Calcb(beta = beta, tm = tm, event = event, ps = ps)
  nabla.eta.shape.Ubeta <- CalcUbetabeeta(beta = beta, tm = tm, event = event, ps = ps, psDeriv = ps.deriv.shape)
  nabla.eta.scale.Ubeta <- CalcUbetabeeta(beta = beta, tm = tm, event = event, ps = ps, psDeriv = ps.deriv.scale)
  nabla.eta.Ubeta <- c(nabla.eta.shape.Ubeta, nabla.eta.scale.Ubeta)/n
  hess.etas.l.v <- (hessian(func = ICweibLik, x = etas, w = w, w.res = w.res))
  grad.eta.pers <- ICweibGrad(etas = etas, w = w, w.res = w.res)
  r.vec <- b.vec - nabla.eta.Ubeta%*%solve(hess.etas.l.v)%*%t(grad.eta.pers)
  meat <- mean(r.vec^2)  # since beta is one-dimensional here 
  bread <- myFmyHess(beta, tm, event, ps)/n
  var.beta <- (meat/(bread^2))/n
  return(var.beta)
}

CalcVarThetaRS <- function(beta, etas.matrix, tm, event, ps.rs, ps.deriv.shape.rs, ps.deriv.scale.rs, w, w.res)
{
  n <- length(tm)
  b.vec <- Calcb(beta = beta, tm = tm, event = event, ps = ps.rs)
  nabla.etas.shape.Ubeta <- CalcUbetabeetaRS(beta = beta, tm = tm, event = event, ps = ps.rs, psDeriv = ps.deriv.shape.rs)
  nabla.etas.scale.Ubeta <- CalcUbetabeetaRS(beta = beta, tm = tm, event = event, ps = ps.rs, psDeriv = ps.deriv.scale.rs)
  nabla.etas.Ubeta <- c(rbind(nabla.etas.shape.Ubeta, nabla.etas.scale.Ubeta))/n

  hess.eta.inv <- ICweibHessSolvedRS(etas.matrix = etas.matrix, w = w, w.res = w.res, tm = tm, event = event)
  grad.eta.pers <- ICweibGradRS(etas = etas.matrix, w = w, w.res = w.res,  tm = tm, event = event)
   r.vec <- b.vec - nabla.etas.Ubeta%*%hess.eta.inv%*%t(grad.eta.pers)
  meat <- mean(r.vec^2)  # since beta is one-dimensional here 
  bread <- myFmyHess(beta, tm, event, ps.rs)/n
  var.beta <- (meat/(bread^2))/n
  return(var.beta)
}
CalcVarEta <- function(etas,  w, w.res)
{
  n <- nrow(w)
  hess.etas.l.v <- (hessian(func = ICweibLik, x = etas, w = w, w.res = w.res))
  grad.eta.pers <- ICweibGrad(etas = etas, w = w, w.res = w.res)
  grad.eta <- 0
  for (j in 1:nrow(grad.eta.pers))
  {
    grad.eta <- grad.eta + grad.eta.pers[j,]%*%t(grad.eta.pers[j,])
  }    
  var.eta <-  solve(hess.etas.l.v)%*%(grad.eta)%*%solve(hess.etas.l.v)
    return(var.eta)
}
CalcVarNpmle <- function(tm, event, w, w.res, BS = 100)
{
  n <- length(tm)
  beta.np.calib.bs <- vector(length = BS)
  for (j in 1:BS)
  {
  #  cat("j = ", j)
    indices <- sample.int(n = n, size = n, replace = T)
    tm.bs <- tm[indices]
    event.bs <- event[indices]
    case.times.bs <- tm.bs[event.bs==1]
    w.bs <- w[indices,]
    w.res.bs <- w.res[indices,]
    fit.npmle.bs <- FitCalibNpmle(w = w.bs, w.res = w.res.bs)
    px.np.bs <- t(sapply(case.times.bs, CalcNpmleCalibP, w = w.bs, w.res =  w.res.bs, fit.npmle = fit.npmle.bs))
    beta.np.calib.bs[j] <- optimize(f = myF,  tm = tm.bs, event = event.bs, ps = px.np.bs, 
                             interval = c(-50,50), maximum = T)$maximum
  }
  v.hat.npmle <- var(beta.np.calib.bs)
  return(v.hat.npmle)
}
CalcVarNpmleRS <- function(tm, event, w, w.res, BS = 100)
{
  n <- length(tm)
  beta.np.calib.rs.bs <- vector(length = BS)
  for (j in 1:BS)
  {
    #  cat("j = ", j)
    indices <- sample.int(n = n, size = n, replace = T)
    tm.bs <- tm[indices]
    event.bs <- event[indices]
    case.times.bs <- tm.bs[event.bs==1]
    w.bs <- w[indices,]
    w.res.bs <- w.res[indices,]
    px.np.rs.bs <- t(sapply(case.times.bs, CalcNpmleRiskSetP, w = w.bs, w.res =  w.res.bs, obs.tm = tm.bs))
    beta.np.calib.rs.bs[j] <- optimize(f = myF,  tm = tm.bs, event = event.bs, ps = px.np.rs.bs, 
                                    interval = c(-50,50), maximum = T)$maximum
  }
  v.hat.npmle <- var(beta.np.calib.bs)
  return(v.hat.npmle)
}
