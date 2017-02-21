#include <RcppArmadillo.h>
using namespace Rcpp;


//
// tm - event time
// event - censoring indicator (1 event 0 no event)
// ps - probabilities
// theta- a vector of coefficients, the first is beta

//data is not assumed to be sorted by time

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::mat CoxLogLikHess(arma::vec theta, arma::vec tm, arma::vec event, arma::mat ps, arma::mat Q) {
  int n = tm.size();
  int sumD = sum(event);
  int nPars = theta.size();
  int iCaseNum=-1;
  int jCaseNum=-1;
  
  double beta = theta[0];
  double NablaBetaUbeta=0;
  double DenomNablaGammaUgamma = 0;
  
  arma::vec gamma = theta.subvec(1,nPars-1);
  arma::vec Szero =  arma::zeros(sumD);
  arma::vec a =  arma::zeros(sumD);
  arma::vec GamQ = Q * gamma;
  arma::vec ExpGamQ = exp(Q * gamma);
  
  arma::mat Sone = arma::zeros(sumD,nPars);
  arma::mat StwoGammaBeta =  arma::zeros(sumD,nPars-1);
  arma::mat contrib=1 + ps*(exp(beta)-1);
  arma::mat nu = 1 + ps*(exp(beta)-1);
  arma::mat nuDerivBeta= ps*exp(beta);
  arma::mat nua = 1 - ps;
  arma::mat FirstBetaTermNumer = (nuDerivBeta % (1-ps));
  arma::mat FirstBetaTermDenom = nu%nu;
  arma::mat FirstBetaTerm = FirstBetaTermNumer/FirstBetaTermDenom;
  arma::mat NablaGammaUgamma = arma::zeros(nPars-1,nPars-1);
  arma::mat NumerNablaGammaUgamma = arma::zeros(nPars-1,nPars-1);
  arma::mat NumerNablaGammaUbeta  = arma::zeros(nPars-1,1);
  arma::mat NablaGammaUbeta  = arma::zeros(1,nPars-1);
  arma::mat Hess = arma::zeros(nPars,nPars);
  
  arma::cube StwoGamma = arma::zeros(nPars-1,nPars-1,sumD);
  // First a loop to calculate all sums (e.g., S1 and S1)
  for (int i = 0; i < n; ++i) {
    if (event[i]) {
     iCaseNum += 1;
     Szero[iCaseNum] += nu(iCaseNum,i)*ExpGamQ[i];
     Sone(iCaseNum,0) += nuDerivBeta(iCaseNum,i)*ExpGamQ[i];
     a[iCaseNum] += nua(iCaseNum,i)*ExpGamQ[i];
     arma::mat Qi = Q(i,arma::span::all);
     StwoGammaBeta(iCaseNum,arma::span::all) +=  Qi*ExpGamQ[i]*nuDerivBeta(iCaseNum,i);
     StwoGamma.slice(iCaseNum) += Qi.t() * Qi * Szero[iCaseNum];
     Sone(iCaseNum, arma::span(1,nPars-1)) +=   nu(iCaseNum,i)*ExpGamQ[i]*Qi;
     for(int j = 0; j < n; ++j) {
       if (tm[j]>tm[i]) {
        Szero[iCaseNum] += nu(iCaseNum,j)*ExpGamQ[j];
        Sone(iCaseNum,0) += nuDerivBeta(iCaseNum,j)*ExpGamQ[j];
        a[iCaseNum] += nua(iCaseNum,j)*ExpGamQ[j];
        arma::mat Qj = Q(j,arma::span::all);
        StwoGammaBeta(iCaseNum,arma::span::all) +=  Qj*ExpGamQ[j]*nuDerivBeta(iCaseNum,j);
        StwoGamma.slice(iCaseNum) += Qj.t() * Qj * nu(iCaseNum,j)*ExpGamQ[j];
        Sone(iCaseNum, arma::span(1,nPars-1)) +=   nu(iCaseNum,j)*ExpGamQ[j]*Qj;
       }
      }
    }
  }
  for (int k = 0; k < n; ++k) {
     if (event[k]) {
       jCaseNum += 1;
       NablaBetaUbeta += FirstBetaTerm(jCaseNum,k) - (Sone(jCaseNum,0)*a[jCaseNum])/(Szero[jCaseNum]*Szero[jCaseNum]);
       arma::mat SoneTerm = Sone(jCaseNum,arma::span(1,nPars-1));
       NumerNablaGammaUgamma = Szero[jCaseNum]*StwoGamma.slice(jCaseNum) -   SoneTerm.t()*SoneTerm;
       DenomNablaGammaUgamma = Szero[jCaseNum]*Szero[jCaseNum];
       NablaGammaUgamma -= NumerNablaGammaUgamma/DenomNablaGammaUgamma;
       NumerNablaGammaUbeta = Szero[jCaseNum]*StwoGammaBeta(jCaseNum,arma::span::all) - SoneTerm*Sone(jCaseNum,0);
       NablaGammaUbeta -= NumerNablaGammaUbeta/DenomNablaGammaUgamma;
       }         
  }
  Hess(0,0) = NablaBetaUbeta;
  Hess(arma::span(1,nPars-1),arma::span(1,nPars-1)) = NablaGammaUgamma; 
  Hess(0,arma::span(1,nPars-1)) = NablaGammaUbeta;
  Hess(arma::span(1,nPars-1),0) = NablaGammaUbeta.t();
  return Hess;
}


