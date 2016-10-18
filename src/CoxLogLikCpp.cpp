#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
// tm - event time
// event - censoring indicator (1 event 0 no event)
// ps - probabilities
// beta- a value

//data is not assumed to be sorted by time
// [[Rcpp::export]]
double CoxLogLikCpp(double beta, NumericVector tm, LogicalVector event, NumericMatrix ps, 
                 NumericVector Qgamma) {
  int n = tm.size();
  double denom=0;
  double logDenom=0;
  double logNumer=0;
  double logLik=0;
  int iCaseNum=-1;
  NumericMatrix contrib=1 + ps*(exp(beta)-1);
  NumericVector conribQ= exp(Qgamma);
  for (int i = 0; i < n; ++i)
  {
    if (event[i]) {
      iCaseNum += 1;
      logNumer += log(contrib(iCaseNum,i)) + Qgamma[i];
      denom = contrib(iCaseNum,i)*conribQ[i];
      for(int j = 0; j < n; ++j) {
       if (tm[j]>tm[i]) {
        denom += contrib(iCaseNum,j)*conribQ[j];
        }
        }
     logDenom += log(denom);
    } }
  logLik = logNumer - logDenom; 
  
//NumericVector contrib=ps*exp(beta)+1-ps;
//List aa; aa["logLik"] = logLik;
//  return aa;
    return logLik;
  }


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

// /*** R
// timesTwo(42)
// */
