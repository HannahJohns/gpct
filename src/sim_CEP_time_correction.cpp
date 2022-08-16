#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double correctTimes_row(const NumericVector& times,
                        const IntegerMatrix& censorPath,
                        int var, int step, int stepLimit)
{
  if(step > stepLimit)
  {
    // Detects infinite recursion, is caused by cycles in the censorPath matrix
    return -999;
  }

  int nVars = censorPath.ncol();

  double censoredTime = std::numeric_limits<double>::infinity();

  for(int i=0; i<nVars; i++)
  {
    if(censorPath(i,var) == 1)
    {
      double dependsTime = correctTimes_row(times, censorPath,i, step+1, stepLimit);
      if(dependsTime < censoredTime) censoredTime = dependsTime;
    }
  }

  return std::min(censoredTime,times[var]);
}



// [[Rcpp::export]]
void correctTimes(NumericMatrix& out,
                  const NumericMatrix& times,
                  const IntegerMatrix& censorPath){


  int nObs = times.nrow();
  int nVars = times.ncol();

  for(int i=0;i<nObs; i++)
  {
    NumericVector thisRow = times(i,_);

    for(int j=0; j<nVars;j++)
    {
      out(i,j) = correctTimes_row(thisRow,censorPath,j,0,nVars);
    }

  }

  //return(out);
}



