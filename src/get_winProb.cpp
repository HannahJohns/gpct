// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>

using namespace RcppParallel;
using namespace Rcpp;

struct WinProb_getter : public Worker
{

  const bool population;
  RMatrix<double> out;
  const RMatrix<int> pairwise;
  const RVector<double> weight;

  WinProb_getter(NumericMatrix& out,
                 const IntegerMatrix& pairwise,
                 const NumericVector& weight,
                 const bool population
                ) :
    out(out),
    pairwise(pairwise),
    weight(weight),
    population(population)
  {}

  void operator()(std::size_t begin, std::size_t end)
  {
    int nObs = pairwise.nrow();

    for(int i=begin; i<end; i++)
    {
      out[i] = out[i+nObs] = out[i+2*nObs] = 0;

      double sum_weights = 0;
      for(int j=0; j<nObs; j++)
      {

        int indicator=pairwise[i+j*nObs];

        // If we're running as a sample rather than a population, don't compare
        // to ourself.

        // NOTE: This will break terribly if weights are less than 1
        if( (!population) & (i==j) )
        {
          switch (indicator)
          {
          case -1:
            // Loss
            out[i+nObs] = out[i+nObs] + (weight[j]-1);
            sum_weights = sum_weights + (weight[j]-1);
            break;
          case 0:
            // Tie
            out[i+2*nObs] = out[i+2*nObs] + (weight[j]-1);
            sum_weights = sum_weights + (weight[j]-1);
            break;
          case 1:
            // Win
            out[i] = out[i] + (weight[j]-1);
            sum_weights = sum_weights + (weight[j]-1);
            break;
          default:
            // Not a valid set of pairs, add weights but don't add Rd
            sum_weights = sum_weights + (weight[j]-1);
          break;
          }

        }
        else
        {
          switch (indicator)
          {
          case -1:
            // Loss
            out[i+nObs] = out[i+nObs] + weight[j];
            sum_weights = sum_weights + weight[j];
            break;
          case 0:
            // Tie
            out[i+2*nObs] = out[i+2*nObs] + weight[j];
            sum_weights = sum_weights + weight[j];
            break;
          case 1:
            // Win
            out[i] = out[i] + weight[j];
            sum_weights = sum_weights + weight[j];
            break;
          default:
            // Not a valid set of pairs, add weights but don't add Rd
            sum_weights = sum_weights + weight[j];
          break;
          }
        }


      } // End for j

      // Standardize Rs, Rd and Rt
      for(int k=0; k<3; k++) out[i+k*nObs] = out[i+k*nObs]/sum_weights;

    } // End for i
  }

};


// [[Rcpp::export]]
NumericMatrix winProb(const IntegerMatrix& pairwise,
                      const NumericVector& weight,
                      const bool population)
{
  // Individual win probability

  int nObs = pairwise.nrow();

  NumericMatrix out(nObs,3);

  WinProb_getter winProb(out,pairwise,weight,population);

  parallelFor(0,nObs,winProb);

  return out;
}
