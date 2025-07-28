// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>

using namespace RcppParallel;
using namespace Rcpp;

struct RSDTgetter : public Worker
{

  RMatrix<double> out;
  const RMatrix<int> pairwise;
  const RVector<double> group;
  const RVector<double> weight;
  const double tol;

  RSDTgetter(NumericMatrix& out,
             const IntegerMatrix& pairwise,
             const NumericVector& group,
             const NumericVector& weight,
             const double tol) :
    out(out),
    pairwise(pairwise),
    group(group),
    weight(weight),
    tol(tol)
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
        int indicator = 99; // Arbitrary code, not -1, 0 or 1

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
          out[i+nObs] = out[i+nObs] + weight[j];
          sum_weights = sum_weights + weight[j];
          break;
        case 0:
          // Tied Pair
          out[i+2*nObs] = out[i+2*nObs] + weight[j];
          sum_weights = sum_weights + weight[j];
          break;
        case 1:
          // Concordant Pair
          out[i] = out[i] + weight[j];
          sum_weights = sum_weights + weight[j];
          break;
        default:
          // Not a valid set of pairs, add weights but don't add Rd
          sum_weights = sum_weights + weight[j];
        break;
        }
      }

      // Standardize Rs, Rd and Rt
      for(int k=0; k<3; k++) out[i+k*nObs] = out[i+k*nObs]/sum_weights;
    }

  }
};



// [[Rcpp::export]]
NumericMatrix getRsdt(const IntegerMatrix& pairwise,
                      const NumericVector& group,
                      const NumericVector& weight,
                      const double tol)
{
  // Get probability of concordant,discordant or tied for each record

  int nObs = pairwise.nrow();

  NumericMatrix out(nObs,3);

  RSDTgetter rsdt(out,pairwise,group,weight,tol);

  parallelFor(0,nObs,rsdt);

  return out;
}



/*******************************************************************************
 *
 *     Rs Rd and Rt under the null hypothesis of no effect
 *     Generalisation of O'brien's WMW odds method
 *
 *******************************************************************************/

struct NullRSDTgetter : public Worker
{

  RMatrix<double> Rs;
  RMatrix<double> Rd;
  RMatrix<double> Rt;
  RMatrix<double> P; // Probability of being in this cell

  const RMatrix<double> winProb_X;
  const RMatrix<double> winProb_Y;
  const RVector<double> weight_X;
  const RVector<double> weight_Y;

  NullRSDTgetter(NumericMatrix& Rs,
                 NumericMatrix& Rd,
                 NumericMatrix& Rt,
                 NumericMatrix& P,
                 const NumericMatrix& winProb_X,
                 const NumericMatrix& winProb_Y,
                 const NumericVector& weight_X,
                 const NumericVector& weight_Y
  ) :
    Rs(Rs),
    Rd(Rd),
    Rt(Rt),
    P(P),
    winProb_X(winProb_X),
    winProb_Y(winProb_Y),
    weight_X(weight_X),
    weight_Y(weight_Y)
  {}

  void operator()(std::size_t begin, std::size_t end)
  {
    int nX = winProb_X.nrow();
    int nY = winProb_Y.nrow();

    // Get marginal values for
    for(int i=begin; i<end; i++) // for X
    {
      for(int j=0; j<nY; j++) // for Y
      {
        // Rs <- probWin_Y[,1] * probWin_X[k,1] + probWin_Y[,2] * probWin_X[k,2]
        // Rd <- probWin_Y[,1] * probWin_X[k,2] + probWin_Y[,2] * probWin_X[k,1]
        // Rt <- probWin_Y[,3] * (1-probWin_X[k,3])

        Rs[i+nX*j] = winProb_X[i+nX*0] * winProb_Y[j+nY*0] + winProb_X[i+nX*1] * winProb_Y[j+nY*1];
        Rd[i+nX*j] = winProb_X[i+nX*0] * winProb_Y[j+nY*1] + winProb_X[i+nX*1] * winProb_Y[j+nY*0];
        Rt[i+nX*j] = (1-winProb_X[i+nX*2]) * winProb_Y[j+nY*2];
        P[i+nX*j] = weight_X[i] * weight_Y[j];
      }
    }
  }
};

// [[Rcpp::export]]
List get_nullSE(const NumericMatrix& winProb_X,
                const NumericMatrix& winProb_Y,
                const NumericVector& weight_X,
                const NumericVector& weight_Y
)
{

  int nX = winProb_X.nrow();
  int nY = winProb_Y.nrow();

  NumericMatrix Rs(nX,nY);
  NumericMatrix Rd(nX,nY);
  NumericMatrix Rt(nX,nY);
  NumericMatrix P(nX,nY);

  NullRSDTgetter nullRSDT(Rs,Rd,Rt,P,winProb_X,winProb_Y,weight_X,weight_Y);
  parallelFor(0,nX,nullRSDT);


  double sum_weight = 0;
  for(int i=0;i<nX;i++) sum_weight += weight_X(i);
  P = P/sum_weight;

  List out = List::create(Named("Rs") = Rs,
                          Named("Rd") = Rd,
                          Named("Rt") = Rt,
                          Named("P") = P
  );

  return out;
}


/*******************************************************************************
*
*     Probability of two correlated outcomes, used in genodds derived from
*     Kendal's tau derived
*
*******************************************************************************/


struct RCCgetter : public Worker
{

  RMatrix<double> sum_Ncc;

  const RMatrix<int> pairwise;
  const RVector<double> group;
  const RVector<double> weight;
  const double tol;

  RCCgetter(NumericMatrix& sum_Ncc,
            const IntegerMatrix& pairwise,
            const NumericVector& group,
            const NumericVector& weight,
            const double tol) :
    sum_Ncc(sum_Ncc),
    pairwise(pairwise),
    group(group),
    weight(weight),
    tol(tol)
  {}

  void operator()(std::size_t begin, std::size_t end)
  {
    int nObs = pairwise.nrow();

    for(int i=begin; i<end; i++)
    {
      double sum_weights = 0;

      for(int j=0; j<nObs;j++)
      {
        if(i==j) continue;

        for(int k=j+1;k<nObs;k++)
        {
          if(i==k) continue;

          int indicator_1 = 99; // Arbitrary code, larger than 1
          if(group[i] > group[j] + tol)
          {
            indicator_1=pairwise[i+j*nObs];
          }
          else if (group[j] > group[i] + tol)
          {
            indicator_1=pairwise[j+i*nObs];
          }

          int indicator_2 = 99; // Arbitrary code, larger than 1
          if(group[i] > group[k] + tol)
          {
            indicator_2=pairwise[i+k*nObs];
          }
          else if (group[k] > group[i] + tol)
          {
            indicator_2=pairwise[k+i*nObs];
          }

          sum_weights = sum_weights + weight[j]*weight[k];

          if((indicator_1 == 1) & (indicator_2 == 1) )
          {
            sum_Ncc[i] = sum_Ncc[i] +  weight[j]*weight[k];
          }
          else if( ( (indicator_1 == 1) & (indicator_2 == 0) ) |
                   ( (indicator_1 == 0) & (indicator_2 == 1) )
                 )
          {
            sum_Ncc[i+nObs] = sum_Ncc[i+nObs] + weight[j]*weight[k];
          }
          else if( (indicator_1 == 0) & (indicator_2 == 0))
          {
            sum_Ncc[i+2*nObs] = sum_Ncc[i+2*nObs] +  weight[j]*weight[k];
          }
        }
      }

      sum_Ncc[i] = sum_Ncc[i]/sum_weights;
      sum_Ncc[i+nObs] = sum_Ncc[i+nObs]/sum_weights;
      sum_Ncc[i+2*nObs] = sum_Ncc[i+2*nObs]/sum_weights;
    }
  }
};



// [[Rcpp::export]]
NumericMatrix get_Rcc(const IntegerMatrix& pairwise,
                      const NumericVector& group,
                      const NumericVector& weight,
                      const double tol)
{

  // Get probability of:
  // Probability of two concordant, one concordant one tied, and two tied

  int nObs = pairwise.nrow();

  NumericMatrix out(nObs,3);

  RCCgetter rcc_calc(out, pairwise, group, weight, tol);

  parallelFor(0, nObs, rcc_calc);

  return out;

}
