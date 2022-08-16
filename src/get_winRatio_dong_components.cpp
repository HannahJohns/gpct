// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>

using namespace RcppParallel;
using namespace Rcpp;


struct WinRatioStatsGetter_ntnc : public Worker
{

  // Calculates Nt and Nc for Dong's win ratio formula

  double nt; // Number of treatment wins (equation 4)
  double nc; // Number of control wins (analogous to 4)

  const RMatrix<int> pairwise;
  const RVector<int> treatList;
  const RVector<int> controlList;
  const bool split;


  WinRatioStatsGetter_ntnc(const IntegerMatrix& pairwise,
                           const IntegerVector& treatList,
                           const IntegerVector& controlList,
                           const bool split) :
    nt(0),
    nc(0),
    pairwise(pairwise),
    treatList(treatList),
    controlList(controlList),
    split(split)
  {}

  WinRatioStatsGetter_ntnc(WinRatioStatsGetter_ntnc& getter, Split) :
    nt(0),
    nc(0),
    pairwise(getter.pairwise),
    treatList(getter.treatList),
    controlList(getter.controlList),
    split(getter.split)
  {}


  void operator()(std::size_t begin, std::size_t end)
  {
    int Nt = treatList.size(); // Count of treatments
    int Nc = controlList.size(); // Count of controls
    int nObs = pairwise.nrow();

    // First, get nt and nc
    for(int i=begin; i<end; i++)
    {
      for(int j=0; j<Nc; j++)
      {
        int result = pairwise[treatList[i]+controlList[j]*nObs];

        switch(result)
        {
        case -1:
          nc++;
          break;
        case 0:
          if(split) {nc+=0.5; nt+=0.5;} // If we split ties, add half to each of nt and nc
          break;
        case 1:
          nt++;
          break;
        break;
        }
      }
    }
  }

  void join(const WinRatioStatsGetter_ntnc& rhs)
  {
    nt += rhs.nt;
    nc += rhs.nc;
  }

};


struct WinRatioStatsGetter_stderr : public Worker
{

  // Values we sum over

  double s_t1; // Equation 5d
  double s_t2; // Equation 5e
  double s_c1; // Analogous to 5d
  double s_c2; // Analogous to 5d

  double s_t1null; // Equation 8b
  double s_t2null; // Equation 8c
  double s_c1null; // Analogous to 8b
  double s_c2null; // Analogous to 8c

  double cov1; // Not in paper explicitly
  double cov2; // Not in paper explicitly
  double cov1null; // Equation 11b
  double cov2null; // Equation 11c


  // Results calculated from the first step

  double theta_k;
  double theta_l;
  double theta_knull;
  double theta_lnull;

  const RMatrix<int> pairwise;
  const RVector<int> treatList;
  const RVector<int> controlList;
  const bool split;


  WinRatioStatsGetter_stderr(const IntegerMatrix& pairwise,
                             const IntegerVector& treatList,
                             const IntegerVector& controlList,
                             const bool split,
                             double theta_k,
                             double theta_l,
                             double theta_knull,
                             double theta_lnull
  ):
    pairwise(pairwise),
    treatList(treatList),
    controlList(controlList),
    split(split),
    theta_k(theta_k),
    theta_l(theta_l),
    theta_knull(theta_knull),
    theta_lnull(theta_lnull),
    s_t1(0),
    s_t2(0),
    s_c1(0),
    s_c2(0),
    s_t1null(0),
    s_t2null(0),
    s_c1null(0),
    s_c2null(0),
    cov1(0),
    cov2(0),
    cov1null(0),
    cov2null(0)
  {}

  WinRatioStatsGetter_stderr(WinRatioStatsGetter_stderr& getter, Split):
    pairwise(getter.pairwise),
    treatList(getter.treatList),
    controlList(getter.controlList),
    split(getter.split),
    theta_k(getter.theta_k),
    theta_l(getter.theta_l),
    theta_knull(getter.theta_knull),
    theta_lnull(getter.theta_lnull),
    s_t1(0),
    s_t2(0),
    s_c1(0),
    s_c2(0),
    s_t1null(0),
    s_t2null(0),
    s_c1null(0),
    s_c2null(0),
    cov1(0),
    cov2(0),
    cov1null(0),
    cov2null(0)
  {}

  void operator()(std::size_t begin, std::size_t end)
  {

    int Nt = treatList.size(); // Count of treatments
    int Nc = controlList.size(); // Count of controls
    int nObs = pairwise.nrow();

    // All part 1 values are of the form sum_i(sum_j(sum_j'=/=j(....)))
    for(int i=begin;i<end;i++)
    {
      for(int j=0;j<Nc;j++)
      {
        for(int jj=0;jj<Nc;jj++) // jj is j'
        {
          if(j!=jj)
          {

            // Get results of kernel functions

            // K(i,j) and L(i,j)
            int result =  pairwise[treatList[i]+controlList[j]*nObs];

            double K=0;
            double L=0;

            switch(result)
            {
            case -1:
              K=0;
              L=1;
              break;
            case 0:
              if(split) {K=0.5; L=0.5;} // If we split ties, It's half a win/half a loss
              break;
            case 1:
              K=1;
              L=0;
              break;
            default:
              break;
            }

            // K(i,j') and L(i,j')
            int resultPrime =  pairwise[treatList[i]+controlList[jj]*nObs];

            double KPrime=0;
            double LPrime=0;

            switch(resultPrime)
            {
            case -1:
              KPrime=0;
              LPrime=1;
              break;
            case 0:
              if(split) {KPrime=0.5; LPrime=0.5;} // If we split ties, It's half a win/half a loss
              break;
            case 1:
              KPrime=1;
              LPrime=0;
              break;
            default:
              break;
            }

            // Add contribution to triple sums

            s_t1 += (K - theta_k)*(KPrime - theta_k);
            s_c1 += (L - theta_l)*(LPrime - theta_l);
            cov1 += (K - theta_k)*(LPrime - theta_l);

            s_t1null += (K - theta_knull)*(KPrime - theta_knull);
            s_c1null += (L - theta_lnull)*(LPrime - theta_lnull);
            cov1null += (K - theta_knull)*(LPrime - theta_lnull);

          }
        }
      }
    }

    // All part 2 values are of the form sum_j(sum_i(sum_i'=/=i(....)))

    // These first two outer sums are independent so we can just
    // swap them around. This lets us calculate both parts.

    for(int i=begin;i<end;i++)
    {
      for(int j=0;j<Nc;j++)
      {
        for(int ii=0;ii<Nt;ii++)
        {
          if(i!=ii)
          {
            // Get results of kernel functions

            // K(i,j) and L(i,j)
            int result = pairwise[treatList[i]+controlList[j]*nObs];

            double K=0;
            double L=0;

            switch(result)
            {
            case -1:
              K=0;
              L=1;
              break;
            case 0:
              if(split) {K=0.5; L=0.5;} // If we split ties, It's half a win/half a loss
              break;
            case 1:
              K=1;
              L=0;
              break;
            default:
              break;
            }

            // K(i',j) and L(i',j)
            int resultPrime = pairwise[treatList[ii]+controlList[j]*nObs];

            double KPrime=0;
            double LPrime=0;

            switch(resultPrime)
            {
            case -1:
              KPrime=0;
              LPrime=1;
              break;
            case 0:
              if(split) {KPrime=0.5; LPrime=0.5;} // If we split ties, It's half a win/half a loss
              break;
            case 1:
              KPrime=1;
              LPrime=0;
              break;
            default:
              break;
            }

            // Add contribution to triple sums

            s_t2 += (K - theta_k)*(KPrime - theta_k);
            s_c2 += (L - theta_l)*(LPrime - theta_l);
            cov2 += (K - theta_k)*(LPrime - theta_l);

            s_t2null += (K - theta_knull)*(KPrime - theta_knull);
            s_c2null += (L - theta_lnull)*(LPrime - theta_lnull);
            cov2null += (K - theta_knull)*(LPrime - theta_lnull);

          }
        }
      }
    }

  }

  void join(const WinRatioStatsGetter_stderr& rhs)
  {

    s_t1 += rhs.s_t1;
    s_t2 += rhs.s_t2;
    s_c1 += rhs.s_c1;
    s_c2 += rhs.s_c2;

    s_t1null += rhs.s_t1null;
    s_t2null += rhs.s_t2null;
    s_c1null += rhs.s_c1null;
    s_c2null += rhs.s_c2null;

    cov1 += rhs.cov1;
    cov2 += rhs.cov2;
    cov1null += rhs.cov1null;
    cov2null += rhs.cov2null;

  }


};



// [[Rcpp::export]]
void get_winRatioStats(NumericVector& out,
                       const IntegerMatrix& pairwise,
                       const IntegerVector& treatList,
                       const IntegerVector& controlList,
                       const bool split
                       )
{

  // Calculate various statistics used in win ratio calculations (Dong's 2016 paper)


  int Nt = treatList.size(); // Count of treatments
  int Nc = controlList.size(); // Count of controls
  int nObs = pairwise.nrow();


  WinRatioStatsGetter_ntnc ntnc(pairwise,
                                treatList,
                                controlList,
                                split);
  parallelReduce(0,Nt,ntnc);

  double nt=ntnc.nt;
  double nc=ntnc.nc;

  // From nt and nc, get theta_k and theta_l (U statistics)
  // as well as theta_k0 and theta_l0 (U statistics under null)

  double theta_k = nt/((double)Nt*Nc); //Equation 5b
  double theta_l = nc/((double)Nt*Nc); //Analogous to 5b

  double theta_knull = (nt + nc)/((double) 2*Nt*Nc); // Equation 6
  double theta_lnull = theta_knull; // Analogous to Equation 6

  // With point estimates done, we need to get standard deviations

  // Default values in case we can't calculate the standard error
  double st = 0;
  double sc = 0;
  double stnull = 0;
  double scnull = 0;
  double cov = 0;
  double covnull = 0;

  if(Nt > 1 && Nc >1)
  {
    WinRatioStatsGetter_stderr wr_stderr(pairwise,
                                       treatList,
                                       controlList,
                                       split,
                                       theta_k,
                                       theta_l,
                                       theta_knull,
                                       theta_lnull);

    parallelReduce(0,Nt,wr_stderr);

    // Sums are calculated, extract values we care about
    // and add front multiplication terms

    double s_t1 =  Nt*Nc/(Nc-1) * wr_stderr.s_t1;
    double s_c1 =  Nt*Nc/(Nc-1) * wr_stderr.s_c1;
    double cov1 =  Nt*Nc/(Nc-1) * wr_stderr.cov1;

    double s_t1null =  Nt*Nc/(Nc-1) * wr_stderr.s_t1null;
    double s_c1null =  Nt*Nc/(Nc-1) * wr_stderr.s_c1null;
    double cov1null =  Nt*Nc/(Nc-1) * wr_stderr.cov1null;

    double s_t2 =  Nt*Nc/(Nt-1) * wr_stderr.s_t2;
    double s_c2 =  Nt*Nc/(Nt-1) * wr_stderr.s_c2;
    double cov2 =  Nt*Nc/(Nt-1) * wr_stderr.cov2;

    double s_t2null =  Nt*Nc/(Nt-1) * wr_stderr.s_t2null;
    double s_c2null =  Nt*Nc/(Nt-1) * wr_stderr.s_c2null;
    double cov2null =  Nt*Nc/(Nt-1) * wr_stderr.cov2null;

    // Combine parts together to get variance and covariance

    st = s_t1/Nt+s_t2/Nc;
    sc = s_c1/Nt+s_c2/Nc;
    stnull = s_t1null/Nt+s_t2null/Nc;
    scnull = s_c1null/Nt+s_c2null/Nc;
    cov = cov1/Nt+cov2/Nc;
    covnull = cov1null/Nt+cov2null/Nc;
  }


  // Pass these values back to R

  out(0) = nt;
  out(1) = nc;
  out(2) = st;
  out(3) = sc;
  out(4) = cov;
  out(5) = stnull;
  out(6) = scnull;
  out(7) = covnull;
}
