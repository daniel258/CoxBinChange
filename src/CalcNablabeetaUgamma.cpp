#include <RcppArmadillo.h>
using namespace Rcpp;

// tm - event time
// event - censoring indicator (1 event 0 no event)
// ps - probabilities
// psderiv - A cube (R array maybe) the last dimention is the eta paramter so the first slice is the matrix of derivatives with 
// respect to eta_1, for all people and all event times.
// beta- a value
//data is not assumed to be sorted by time
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat CalcNablabeetaUgamma(arma::vec theta, arma::vec tm, arma::vec event, arma::mat ps, arma::mat Q, arma::mat psDeriv) {
  int n = tm.size();
  int sumD = sum(event);
  int nPars = theta.size();
//  int nEta = psDeriv.n_slices;
  int iCaseNum=-1;
  int jCaseNum=-1;
  // double FirstTerm=0;
  // double SecondTerm=0;
  // double FirstSumType=0;
  // double SecondSumType=0;
  // double ThirdSumType=0;
  double beta = theta[0];
  
  //arma::vec all=0;
  //double all=0;
  arma::vec gamma = theta.subvec(1,nPars-1);
  arma::vec Szero =  arma::zeros(sumD);
  arma::vec a =  arma::zeros(sumD);
  arma::vec ExpGamQ = exp(Q * gamma);
  arma::vec ExpGamQbeta = exp(Q * gamma + beta);
  arma::mat deriv = arma::zeros(1,nPars-1);
  arma::mat SecondTermNumer = arma::zeros(1,nPars-1);
  arma::mat SzeroEta = arma::zeros(sumD,nPars-1);
  arma::mat SoneEta = arma::zeros(sumD,nPars-1);
  arma::mat StwoGammaBeta =  arma::zeros(sumD,nPars-1);
  arma::mat nu = 1 + ps*(exp(beta)-1);


  // First a loop to calculate all sums (e.g., S0 and S1) that do not involve psDeriv
  for (int i = 0; i < n; ++i) {
    if (event[i]) {
      iCaseNum += 1;
      Szero[iCaseNum] += nu(iCaseNum,i)*ExpGamQ[i];
      arma::mat Qi = Q(i,arma::span::all);
      SzeroEta(iCaseNum, arma::span::all) +=   nu(iCaseNum,i)*ExpGamQbeta[i];
     SoneEta(iCaseNum, arma::span::all) +=   psDeriv(iCaseNum,i)*ExpGamQ[i]*Qi;
  StwoGammaBeta(iCaseNum,arma::span::all) +=  Qi*ExpGamQbeta[i]*psDeriv(iCaseNum,i);
  for(int j = 0; j < n; ++j) {
    if (tm[j]>tm[i]) {
      Szero[iCaseNum] += nu(iCaseNum,j)*ExpGamQ[j];
      arma::mat Qj = Q(j,arma::span::all);
    SzeroEta(iCaseNum, arma::span::all) +=   nu(iCaseNum,j)*ExpGamQbeta[j];
    SoneEta(iCaseNum, arma::span::all) +=   psDeriv(iCaseNum,j)*ExpGamQ[j]*Qj;
      StwoGammaBeta(iCaseNum,arma::span::all) +=  Qj*ExpGamQbeta[j]*psDeriv(iCaseNum,j);
    }
   }
    }
  }

 for (int k = 0; k < n; ++k)
   {
     if (event[k]) {
       jCaseNum += 1;
       SecondTermNumer(0, arma::span::all) = SoneEta(jCaseNum, arma::span::all)%SzeroEta(jCaseNum, arma::span::all);
       deriv(0, arma::span::all) -= (StwoGammaBeta(jCaseNum,arma::span::all)*Szero[jCaseNum]- SecondTermNumer)/(Szero[jCaseNum]*Szero[jCaseNum]);
       }
     }
     return deriv;
   }

