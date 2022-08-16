#include "header.hpp"


// Sign function
//https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}



// Pairwise comparison scores for composite survival data,
// returns codes as follows:

// -1: j beats i
// 1: i beats j
// 0: i and j are tied

// Notes on usage for Dong's (2016) win ratio standard error method:
// To get kernel K, use std::max(indicatorFunction_surv(),0)
// To get kernel L, use std::max(-indicatorFunction_surv(),0)

int indicatorFunction_numeric(int i, int j,
                           const RMatrix<double>& outcome,
                           const RVector<double>& discriminant
)
{

  int result = 0;

  int nVars = outcome.ncol();
  int nObs = outcome.nrow();

  // Run from column 0 to the right in order of importance
  for(int column=0; column < nVars; column++)
  {
    if(outcome[i+column*nObs] > outcome[j+column*nObs] + discriminant[column])
    {
  //    Rprintf("(%d,%d): %f > %f + %f, assign 1\n", i,j, outcome[i+column*nObs],outcome[j+column*nObs], discriminant[column]);
      result = 1;
      break;
    }
    else if (outcome[j+column*nObs] > outcome[i+column*nObs] + discriminant[column])
    {
//      Rprintf("(%d,%d): %f > %f + %f, assign -1\n", i,j, outcome[j+column*nObs], outcome[i+column*nObs], discriminant[column]);
      result=-1;
      break;
    }
  }
  return result;
}

int indicatorFunction_orderedScalar(int i, int j,
                                    const RMatrix<double>& outcome,
                                    const RVector<int>& varType,
                                    const std::vector<std::vector<std::vector<double>>>& comparisonDetails
)
{
  int result = 0;
  int nObs = outcome.nrow();
  int nVars = outcome.ncol();

  for(int column=0; column < nVars; column++)
  {
    if(varType[column]==0) // Continuous
    {
      // comparisonDetails[column] contains a single number giving the discriminant threshold
      // for superiority

      double i_val = outcome[i+column*nObs];
      double j_val = outcome[j+column*nObs];

      if(i_val > j_val + comparisonDetails[column][0][0])
      {
        result = 1;
        break;
      }
      else if (j_val > i_val + comparisonDetails[column][0][0])
      {
        result=-1;
        break;
      }
    }
    else if(varType[column]==1) // Categorical
    {

      int nCat = comparisonDetails[column].size();

      // comparisonDetails[column] gives a nCat x nCat table with rules about
      // wins, losses and stopping criteria

      int i_val = std::floor(outcome[i+column*nObs]+0.2);
      int j_val = std::floor(outcome[j+column*nObs]+0.2);

      int code = std::floor(comparisonDetails[column][i_val][j_val]+0.1);

      if(code==1){
        // i wins
        result=1;
        break;
      } else if (code==-1){
        // j wins
        result=-1;
        break;
      } else if (code==9){
        // no winner and stop attempting to break ties
        break;
      }
    }
  }

  return result;
}

int indicatorFunction_surv(int i, int j,
                           const RMatrix<int>& outcome_exists,
                           const RMatrix<double>& outcome,
                           const RVector<double>& direction,
                           const bool censored_as_tie,
                           const double minDiff
)
{

  int result = 0;

  int nVars = outcome_exists.ncol();
  int nObs = outcome_exists.nrow();

  // Run from column 0 to the right in order of importance
  for(int column=0; column < nVars; column++)
  {

    //Rcout << "Checking Column " << column << std::endl;

    // Does at least one outcome exist for this pair?
    if(outcome_exists(i,column) || outcome_exists(j,column))
    {

      //Rcout << "A win may exist!" << std::endl;

      //If we have an event in both, we can just compare directly
      if(outcome_exists(i,column) && outcome_exists(j,column))
      {
        //Rcout << "We can compare directly on outcome" << std::endl;

        if(abs(outcome(i,column)-outcome(j,column)) > minDiff)
        {
          result = direction[column] * sgn(outcome(i,column)-outcome(j,column));
        }
        else
        {
          // Results are close enoguh together to be considered a tie
          result = 0;
        }
        break;
      }
      else
      {

        // Win or tie depends on what is censored where

        // If exists for i (code 1), then j is censored (code 0)
        int whichCensored = ( outcome_exists(i,column) ? 0 : 1);

        if((whichCensored==1) & (outcome(i,column) >= outcome(j,column)))
        {
          // Rcout << "i has censor time after j's event, i wins" << std::endl;

          //# If i was censored and the censor time was after j's event,
          //# then i had an event after j. i wins if later events are good, return 1

          if(abs(outcome(i,column)-outcome(j,column)) > minDiff)
          {
            result = direction[column];
          }
          else
          {
            // Results are close enoguh together to be considered a tie
            result = 0;
          }
          break;

        }
        else if((whichCensored==0) & (outcome(i,column) <= outcome(j,column)))
        {

          //Rcout << "j has censor time after i's event, j wins" << std::endl;

          //# If j was censored and the censor time was after i's event, then j
          //  had an event after i's event. j wins if later are events are good,
          //  return -1

          if(abs(outcome(i,column)-outcome(j,column)) > minDiff)
          {
            result = -1 * direction[column];
          }
          else
          {
            // Results are close enoguh together to be considered a tie
            result = 0;
          }
          break;

        }
      }

    }


  }
  return result;
}
