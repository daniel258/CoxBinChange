## Daniel Nevo 
#CalcGradEtaPersRSInts <- function(d1, d2, d3, Li, Ri,  knots, order, eta.g, eta.b, Z, pts.for.ints)
CalcGradEtaPersRSInts <- function(d1, d2, d3, Li, Ri, Z, fit.cox.rs.ints, pts.for.ints, tm, n.etas.per.fit)
{
  n <- length(Ri)
  n.fits <- length(pts.for.ints)
  #fit.cox.int.one <- fit.cox.rs.ints[[1]]
  #eta.b.one <- fit.cox.int.one$b
  #eta.g.one <- fit.cox.int.one$g
  #n.g <- length(eta.g.one)
  #n.b <- length(eta.b.one)
  #n.pars.ints <- n.b + n.g
  deriv.ell.etas <- matrix(nr = n, nc = sum(n.etas.per.fit), 0)
  for (j in 1:n.fits)
  {
  point <- pts.for.ints[j]
  fit.cox.int <- fit.cox.rs.ints[[j]]
  eta.b <- fit.cox.int$b
  eta.g <- fit.cox.int$g
  n.g <- length(eta.g)
  n.b <- length(eta.b)
  knots <- fit.cox.int$knots
  order <- fit.cox.int$order
  in.risk.set <- tm >= point
  n.set <- sum(in.risk.set)
  Li.int <- Li[in.risk.set]
  Ri.int <- Ri[in.risk.set]
  d1.int <- d1[in.risk.set]
  d2.int <- d2[in.risk.set]
  d3.int <- d3[in.risk.set]
  Z.int <- Z[in.risk.set,]
  expZb <- as.vector(exp(Z.int%*%eta.b))
  # Portion of the code are taken from the ICsurv package
  ti <- c(Li.int[d1.int == 0], Ri.int[d3.int == 0])
  ti.max <- max(ti) + 1e-05
  ti.min <- min(ti) - 1e-05
  bRi <- t(Ispline(x = Ri.int, order = order, knots = knots))
  bLi <- t(Ispline(x = Li.int, order = order, knots = knots))
  GRi <- as.vector(bRi %*% eta.g)
  GLi <- as.vector(bLi %*% eta.g)
  HRi <-  as.vector(GRi*expZb)
  HLi <-  as.vector(GLi*expZb) 
  SRi <- exp(-HRi)
  SLi <- exp(-HLi)
  FRi <- 1-SRi
  FLi <- 1-SLi
  
  term.deriv.etab.d1 <- Z.int*(SRi*HRi/FRi)
  term.deriv.etab.d2 <- Z.int*(SRi*HRi -SLi*HLi)/(SLi-SRi)
  term.deriv.etab.d3 <- -Z.int*HLi
  
  term.deriv.etag.d1 <- bRi*(SRi*expZb/FRi)
  term.deriv.etag.d2 <- expZb*(SRi*bRi -SLi*bLi)/(SLi-SRi)
  term.deriv.etag.d3 <- -bLi*expZb
  
  deriv.ell.etag.int <- matrix(nr = n.set, nc =  n.g)
  deriv.ell.etab.int <- matrix(nr = n.set, nc =  n.b)
  
  deriv.ell.etab.int[d1.int==1,] <- term.deriv.etab.d1[d1.int==1]
  deriv.ell.etab.int[d2.int==1,] <- term.deriv.etab.d2[d2.int==1]
  deriv.ell.etab.int[d3.int==1,] <- term.deriv.etab.d3[d3.int==1]
  
  deriv.ell.etag.int[d1.int==1,] <- term.deriv.etag.d1[d1.int==1]
  deriv.ell.etag.int[d2.int==1,] <- term.deriv.etag.d2[d2.int==1]
  deriv.ell.etag.int[d3.int==1,] <- term.deriv.etag.d3[d3.int==1]
  
  deriv.ell.etas.int <- cbind(deriv.ell.etab.int,deriv.ell.etag.int)
  if (j > 1) {
  deriv.ell.etas[in.risk.set, (sum(n.etas.per.fit[1:(j-1)]) + 1):sum(n.etas.per.fit[1:j])] <- deriv.ell.etas.int 
  } else {
  deriv.ell.etas[in.risk.set, 1:n.etas.per.fit[1]] <- deriv.ell.etas.int 
  }
  }
return(deriv.ell.etas)
}