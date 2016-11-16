#include <RcppArmadillo.h>
using namespace Rcpp;


// tm - event time
// event - censoring indicator (1 event 0 no event)
// ps - Matrix of probabilities: ps_{ij}=P(X_i(t_j)=1|history at time t_j) i=1,..,n j only go over the cases
// Q - n\times SOMETHING matrix with all covariates but X.
// beta- a value corresponding to the effect of X.
// gamma - a vector with the coefficents of Q.

//data is not assumed to be sorted by time
// The function returns b_i in Daniel's notation, a matrix, rows are for observations, columns for paramters.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat CalcbQ(arma::vec theta, arma::vec tm, arma::vec event, arma::mat ps, arma::mat Q){
  int n = tm.size();
  int sumD = sum(event);
  int nPars = theta.size();
  int iCaseNum=-1;
  int jCaseNum=-1;
  int kCaseNum=-1;
  
  double beta = theta[0];
  
  arma::vec gamma = theta.subvec(1,nPars-1);
  arma::vec Szero =  arma::zeros(sumD);
//  arma::vec a =  arma::zeros(sumD);
  arma::vec ExpGamQ = exp(Q * gamma);
  
  arma::mat Sone = arma::zeros(sumD,nPars);
  arma::mat StwoGammaBeta =  arma::zeros(sumD,nPars-1);
  arma::mat contrib=1 + ps*(exp(beta)-1);
  arma::mat nu = 1 + ps*(exp(beta)-1);
  arma::mat nuDerivBeta= ps*exp(beta);
  arma::mat nua = 1 - ps;
  arma::mat b = arma::zeros(n,nPars);
  
  // First a loop to calculate all sums (e.g., S1 and S1)
  for (int i = 0; i < n; ++i) {
    if (event[i]) {
      iCaseNum += 1;
      Szero[iCaseNum] += nu(iCaseNum,i)*ExpGamQ[i];
      Sone(iCaseNum,0) += nuDerivBeta(iCaseNum,i)*ExpGamQ[i];
  //    a[iCaseNum] += nua(iCaseNum,i)*ExpGamQ[i];
      arma::mat Qi = Q(i,arma::span::all);
      Sone(iCaseNum, arma::span(1,nPars-1)) +=   nu(iCaseNum,i)*ExpGamQ[i]*Qi;
      for(int j = 0; j < n; ++j) {
        if (tm[j]>tm[i]) {
          Szero[iCaseNum] += nu(iCaseNum,j)*ExpGamQ[j];
          Sone(iCaseNum,0) += nuDerivBeta(iCaseNum,j)*ExpGamQ[j];
  //        a[iCaseNum] += nua(iCaseNum,j)*ExpGamQ[j];
          arma::mat Qj = Q(j,arma::span::all);
          Sone(iCaseNum, arma::span(1,nPars-1)) +=   nu(iCaseNum,j)*ExpGamQ[j]*Qj;
        }
      }
    }
  }
  for (int k = 0; k < n; ++k) {
    if (event[k]) {
      jCaseNum += 1;
      b(k,0) += nuDerivBeta(jCaseNum,k)/nu(jCaseNum,k) - Sone(jCaseNum,0)/Szero[jCaseNum];
      b(k,arma::span(1,nPars-1)) += Q(k,arma::span::all) - Sone(jCaseNum,arma::span(1,nPars-1))/Szero[jCaseNum];
      }
    kCaseNum = -1;
    for(int l = 0; l < n; ++l) {
      if (event[l]) {
        kCaseNum +=1;
        if (tm[k]>tm[l]) {
          b(k,0) -= (n*ExpGamQ[k]*nu(kCaseNum,k)/Szero[kCaseNum]) * (nuDerivBeta(kCaseNum,k)/nu(kCaseNum,k) -
          Sone(kCaseNum,0)/Szero[kCaseNum])/n;
          b(k,arma::span(1,nPars-1)) -= (n*ExpGamQ[k]*nu(kCaseNum,k)/Szero[kCaseNum]) * (Q(k,arma::span::all) -
              Sone(kCaseNum,arma::span(1,nPars-1))/Szero[kCaseNum])/n;
          }
        }
      }
    }
  return b;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

// /*** R
// timesTwo(42)
// */
