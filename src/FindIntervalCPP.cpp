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

//data is assumed to be sorted by time
// [[Rcpp::export]]
IntegerVector FindIntervalCPP(double point, NumericMatrix w) {
  int nrow = w.nrow();
  int ncol = w.ncol();
  
  IntegerVector intervalW(nrow);
  bool CondLocation;
    for (int i = 0; i < nrow; ++i)
  {
      int j=0;
      CondLocation = true;
      while(CondLocation == true)
      {
      if (w(i,j) > point)
      {
        if(j==0)
        {
          intervalW[i] = 1;
          CondLocation = false;
        } else
          {     
          intervalW[i] = j + 1;
          CondLocation = false;
          }
        } else if (j==ncol) {
        intervalW[i] = ncol+1;
        CondLocation = false;}
      j += 1;
      }
  }
        
    return intervalW;
  }


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

// /*** R
// timesTwo(42)
// */
