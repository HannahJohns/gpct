// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>

using namespace RcppParallel;
using namespace Rcpp;

struct WinRatioGetter : public Worker
{

  double sum_conc;
  double sum_tie;
  double sum_disc;
  double sum_weight;

  const RMatrix<int> pairwise;
  const RVector<double> group;
  const RVector<double> weight;
  const double tol;

  WinRatioGetter(const IntegerMatrix& pairwise,
             const NumericVector& group,
             const NumericVector& weight,
             const double tol) :
    sum_conc(0),
    sum_tie(0),
    sum_disc(0),
    sum_weight(0),
    pairwise(pairwise),
    group(group),
    weight(weight),
    tol(tol)
  {}

  WinRatioGetter(const WinRatioGetter& wr, Split) :
    sum_conc(0),
    sum_tie(0),
    sum_disc(0),
    sum_weight(0),
    pairwise(wr.pairwise),
    group(wr.group),
    weight(wr.weight),
    tol(wr.tol)
  {}

  void operator()(std::size_t begin, std::size_t end)
  {
    int nObs = pairwise.nrow();

    for(int i=begin; i<end; i++)
    {
      for(int j=0; j<nObs; j++)
      {
        int indicator = 99; // Arbitrary code, larger than 1

        if(group[i] > group[j] + tol)
        {
          indicator=pairwise[i+j*nObs];
        }
        else if (group[j] > group[i] + tol)
        {
          indicator=pairwise[j+i*nObs];
        }

        switch (indicator)
        {
        case -1:
          // Discordant Pair
          sum_disc += weight[i]*weight[j];
          sum_weight += weight[i]*weight[j];
          break;
        case 0:
          // Tied Pair
          sum_tie += weight[i]*weight[j];
          sum_weight += weight[i]*weight[j];
          break;
        case 1:
          // Concordant Pair
          sum_conc += weight[i]*weight[j];
          sum_weight += weight[i]*weight[j];
          break;
        default:
          // Not a valid set of pairs, add weights but don't add Rd
          sum_weight += weight[i]*weight[j];
        break;
        }
      }
    }
  }

  void join(const WinRatioGetter& rhs)
  {
    sum_conc += rhs.sum_conc;
    sum_tie += rhs.sum_tie;
    sum_disc += rhs.sum_disc;
    sum_weight += rhs.sum_weight;
  }
};




// [[Rcpp::export]]
List get_permutationTest(const IntegerMatrix& pairwise,
                      const NumericVector& group,
                      const NumericVector& weight,
                      const double tol,
                      bool splitTies,
                      int nSims,
                      const IntegerMatrix& shuffles)
{

  // THIS ONLY WORKS IF WEIGHTS ARE ALL 1!

  // Get probability of concordant,discordant or tied for each record

  int nObs = pairwise.nrow();

  WinRatioGetter wr(pairwise,group,weight, tol);
  parallelReduce(0,nObs,wr);

  double result=(splitTies?
                   (wr.sum_conc + 0.5*wr.sum_tie)/(wr.sum_disc + 0.5*wr.sum_tie) :
                   (wr.sum_conc)/(wr.sum_disc)
                   );

  NumericVector result_null(nSims);
  for(int i=0; i<nSims;i++)
  {

    // Create a copy of group
    NumericVector group_null(nObs);
    for(int j=0;j<nObs-1;j++)
    {
      group_null(j) = group(j);
    }

    // Fisher-Yates shuffle of group vector
    // Random values are provided from R so we can control the seed
    // shuffles is generated through something like
    // do.call("cbind",lapply(0:(n-2),function(i){floor(runif(nSims,min=i,max=n))})

    for(int j=0;j<nObs-1;j++)
    {
      // Exchange group_null(j) with group_null(shuffles(i,j))
      double tmp = group_null(j);
      group_null(j) = group_null(shuffles(i,j));
      group_null(shuffles(i,j))=tmp;
    }

    WinRatioGetter wr_null(pairwise,group_null,weight, tol);
    parallelReduce(0,nObs,wr_null);

    result_null(i)=(splitTies?
                     (wr_null.sum_conc + 0.5*wr_null.sum_tie)/(wr_null.sum_disc + 0.5*wr_null.sum_tie) :
                     (wr_null.sum_conc)/(wr_null.sum_disc)
                   );

  }

  List out = List::create(
    Named("result")=result,
    Named("result_null")=result_null
    );

  return out;
}





