## Daniel Nevo 
CalcGradEtaPers <- function(d1, d2, d3, Li, Ri,  knots, order, eta.g, eta.b, Z)
{
  n <- length(Ri)
  n.g <- length(eta.g)
  n.b <- length(eta.b)
  expZb <- as.vector(exp(Z%*%eta.b))
  
  # Portion of the code are taken from the ICsurv package
  ti <- c(Li[d1 == 0], Ri[d3 == 0])
  ti.max <- max(ti) + 1e-05
  ti.min <- min(ti) - 1e-05
  bRi <- t(Ispline(x = Ri, order = order, knots = knots))
  bLi <- t(Ispline(x = Li, order = order, knots = knots))
  GRi <- as.vector(bRi %*% eta.g)
  GLi <- as.vector(bLi %*% eta.g)
  HRi <-  as.vector(GRi*expZb )
  HLi <-  as.vector(GLi*expZb) 
  SRi <- exp(-HRi)
  SLi <- exp(-HLi)
  FRi <- 1-SRi
  FLi <- 1-SLi
  
  term.deriv.etab.d1 <- Z*(SRi*HRi/FRi)
  term.deriv.etab.d2 <- Z*(SRi*HRi -SLi*HLi)/(SLi-SRi)
  term.deriv.etab.d3 <- -Z*HLi
  
  term.deriv.etag.d1 <- bRi*(SRi*expZb/FRi)
  term.deriv.etag.d2 <- expZb*(SRi*bRi -SLi*bLi)/(SLi-SRi)
  term.deriv.etag.d3 <- -bLi*expZb
  
  deriv.ell.etag <- matrix(nr = n, nc =  n.g)
  deriv.ell.etab <- matrix(nr = n, nc =  n.b)
  
  deriv.ell.etab[d1==1,] <- term.deriv.etab.d1[d1==1]
  deriv.ell.etab[d2==1,] <- term.deriv.etab.d2[d2==1]
  deriv.ell.etab[d3==1,] <- term.deriv.etab.d3[d3==1]
  
  deriv.ell.etag[d1==1,] <- term.deriv.etag.d1[d1==1]
  deriv.ell.etag[d2==1,] <- term.deriv.etag.d2[d2==1]
  deriv.ell.etag[d3==1,] <- term.deriv.etag.d3[d3==1]
  
  deriv.ell.etas <- cbind(deriv.ell.etab,deriv.ell.etag)
  return(deriv.ell.etas)
}