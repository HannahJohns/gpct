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


struct PW_Scalar_Getter : public Worker
{

  RMatrix<int> out;
  const RMatrix<double> outcome;
  const RVector<int> varType;
  const std::vector<std::vector<std::vector<double>>> comparisonDetails;


  PW_Scalar_Getter(IntegerMatrix& out,
                   const NumericMatrix& outcome,
                   const IntegerVector& varType,
                   const std::vector<std::vector<std::vector<double>>>& comparisonDetails
  ) :
    out(out),
    outcome(outcome),
    varType(varType),
    comparisonDetails(comparisonDetails)
  {}

  void operator()(std::size_t begin, std::size_t end)
  {

    int nObs = outcome.nrow();

    for(int i=begin; i<end; i++)
    {
      out[i+i*nObs] = 0;

      for(int j=i+1; j<nObs; j++)
      {
        out[i+j*nObs] = indicatorFunction_orderedScalar(i, j,
                                                        outcome,
                                                        varType,
                                                        comparisonDetails
        );

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
IntegerMatrix get_pw_scalar(const NumericMatrix& outcome,
                            const IntegerVector& varType,
                            const List& comparisonDetails
                            )
{


  int nObs = outcome.nrow();
  int nVars = outcome.ncol();

  IntegerMatrix out(nObs,nObs);

  std::vector<std::vector<std::vector<double>>> comparisonDetails_in;

  // Take list and feed it to std data structures
  for(int var_i=0; var_i<nVars; var_i++)
  {

    NumericMatrix tmp_mat = comparisonDetails[var_i];

    std::vector<std::vector<double>> thisMatrix;
    if(varType[var_i]==0)
    {
      std::vector<double> thisColumn;
      thisColumn.push_back(tmp_mat(0,0));
      thisMatrix.push_back(thisColumn);
    }
    else if(varType[var_i]==1)
    {
      int nClasses = tmp_mat.nrow();
      for(int i=0; i<nClasses; i++)
      {
        std::vector<double> thisColumn;
        for(int j=0; j<nClasses; j++)
        {
          thisColumn.push_back(tmp_mat(i,j));
        }
        thisMatrix.push_back(thisColumn);
      }
    }
    comparisonDetails_in.push_back(thisMatrix);
  }


  PW_Scalar_Getter pw(out,outcome,varType,comparisonDetails_in);
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


