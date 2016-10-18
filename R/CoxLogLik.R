#### CoxLogLik function
CoxLogLik <- function(beta.gamma, tm, event, ps, Q)
{
  beta <- beta.gamma[1]
  gamma <- beta.gamma[-1]
  Qgamma <- Q%*%gamma
  CoxLogLikCpp(beta = beta, tm = tm, event = event, ps = ps, Qgamma = Qgamma)
}