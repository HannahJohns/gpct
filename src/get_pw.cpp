#include "header.hpp"

struct PW_Numeric_Getter : public Worker
{

  RMatrix<int> out;
  const RMatrix<double> outcome;
  const RVector<double> discriminant;

  PW_Numeric_Getter(IntegerMatrix& out,
                 const NumericMatrix& outcome,
                 const NumericVector& discriminant) :
    out(out),
    outcome(outcome),
    discriminant(discriminant)
  {}

  void operator()(std::size_t begin, std::size_t end)
  {

    int nObs = outcome.nrow();

    for(int i=begin; i<end; i++)
    {
      out[i+i*nObs] = 0;

      for(int j=i+1; j<nObs; j++)
      {
        out[i+j*nObs] = indicatorFunction_numeric(i, j, outcome, discriminant);

        out[j+i*nObs] = -out[i+j*nObs];
      }
    }
  }
};


struct PW_Surv_Getter : public Worker
{

  RMatrix<int> out;
  const RMatrix<int> outcome_exists;
  const RMatrix<double> outcome;
  const RVector<double> direction;
  const bool censored_as_tie;
  const double minDiff;

  PW_Surv_Getter(IntegerMatrix& out,
                 const IntegerMatrix& outcome_exists,
                 const NumericMatrix& outcome,
                 const NumericVector& direction,
                 const bool censored_as_tie,
                 const double minDiff) :
    out(out),
    outcome_exists(outcome_exists),
    outcome(outcome),
    direction(direction),
    censored_as_tie(censored_as_tie),
    minDiff(minDiff)
  {}

  void operator()(std::size_t begin, std::size_t end)
  {

    int nObs = outcome_exists.nrow();
    for(int i=begin; i<end; i++)
    {
      out[i+i*nObs] = 0;

      for(int j=i+1; j<nObs; j++)
      {
        out[i+j*nObs] = indicatorFunction_surv(i, j, outcome_exists, outcome,
                                               direction, censored_as_tie, minDiff);

        out[j+i*nObs] = -out[i+j*nObs];
      }
    }
  }
};


// [[Rcpp::export]]
IntegerMatrix get_pw_numeric(const NumericMatrix& outcome,
                             const NumericVector& discriminant)
{

  int nObs = outcome.nrow();
  IntegerMatrix out(nObs,nObs);

  PW_Numeric_Getter pw(out,outcome,discriminant);
  parallelFor(0,nObs,pw);

  return out;
}



// [[Rcpp::export]]
IntegerMatrix get_pw_surv(const IntegerMatrix& outcome_exists,
                          const NumericMatrix& outcome,
                          const NumericVector& direction,
                          const bool censored_as_tie,
                          const double minDiff)
{

  int nObs = outcome_exists.nrow();
  IntegerMatrix out(nObs,nObs);

  PW_Surv_Getter pw_surv(out,outcome_exists,outcome,
                         direction, censored_as_tie,
                         minDiff);

  parallelFor(0,nObs,pw_surv);

  return out;
}


