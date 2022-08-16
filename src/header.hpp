// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

// Sign function
//https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
template <typename T> int sgn(T val);

namespace constants
{
  const int INCOMPARABLE_CODE = -9999;
}

int indicatorFunction_numeric(int i, int j,
                              const RMatrix<double>& outcome,
                              const RVector<double>& discriminant
);

int indicatorFunction_orderedScalar(int i, int j,
                                    const RMatrix<double>& outcome,
                                    const RVector<int>& varType,
                                    const std::vector<std::vector<std::vector<double>>>& comparisonDetails
                                   );

int indicatorFunction_surv(int i, int j,
                           const RMatrix<int>& outcome_exists,
                           const RMatrix<double>& outcome,
                           const RVector<double>& direction,
                           const bool censored_as_tie,
                           const double minDiff
);
