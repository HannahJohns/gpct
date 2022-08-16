// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>

using namespace RcppParallel;
using namespace Rcpp;


struct InconsistencyCounter : public Worker
{

  const RMatrix<int> pairwise;
  const RVector<double> weight;

  const int m;

  // For some reason this wouldn't compile if it was an Rvector
  // so we'll use std::vector and convert back later
  std::vector<int> count_involved;
  double count;
  double count_comparable;
  double sum_weight;
  bool includeTies;


  // #########################################################################

  // Recursive DFS search for finding the maximum path between two points.
  // Used for cycle detection by setting the start and end points to be the
  // same.

  // Note that this will not detect sub-cycles that don't fall on this path.

  // If this algorithm starts and ends at 0, it will see
  // cycle A, but not cycle B.

  /*

   0 --- 1 -----> 3-----4
   |     |        ^     |
   |  A  |        |  B  |
   \     |        |     |
    ---- 2 -----> 6-----5

   */


  // If we've identified two paths of the same number of hops,
  // we store the one with the smallest number of tied/bidirectional edges.

  // This means that if algorithm starts and ends at 0, it will see cycle B not cycle A.

  /*

    1 -- 0--> 3
     \ A | B /
      \  |  /
       \ | /
        -2-
   */



  struct PathCount {
    int number_hops;
    int number_ties;
  };

  PathCount dfs_maxPathLength(const RMatrix<int> pairwise,
                              int prev_startPoint,
                              int startPoint, int target,
                              std::vector<int>& unused_nodes,
                              bool ignore_direction,
                              int depth
                              )
  {

    // Rprintf("\n");
    // for(int i=0;i<depth;i++) Rprintf("\t");
    // Rprintf("Searching from %d to %d: ", startPoint, target);

    int nObs = pairwise.nrow();

    PathCount out;
    out.number_hops = 0;
    out.number_ties = 0;


    // Check if target can be reached from current start point without backtracking. OLD CONDITION IS unused_nodes.size()==0
    if(!ignore_direction)
    {
      if(prev_startPoint != target)
      {
        out.number_hops = ( (startPoint != target) &&
          (pairwise[startPoint+nObs*target] == 0 ||
          pairwise[startPoint+nObs*target] == 1)
        ) ? 1 : 0;

        out.number_ties = ( (startPoint != target) &&
          (pairwise[startPoint+nObs*target] == 0)
        ) ? 1 : 0;
      }
    }
    else
    {
      if(prev_startPoint != target)
      {
        out.number_hops = ( (startPoint != target) &&
          (pairwise[startPoint+nObs*target] == -1 ||
          pairwise[startPoint+nObs*target] == 0 ||
          pairwise[startPoint+nObs*target] == 1)
        ) ? 1 : 0;

        out.number_ties = ( (startPoint != target) &&
          (pairwise[startPoint+nObs*target] == 0)
        ) ? 1 : 0;
      }
    }




    if(unused_nodes.size()>0)
    {
      // Iterate through unused_nodes such that startPoint -> unused_nodes[i]

      // unused_nodes[i] becomes the new startPoint and we remove it from unused_nodes

      for(int i=0; i<unused_nodes.size();i++)
      {

        int new_startPoint = unused_nodes[i];

        // Skip this start point if it's not connected to the current startPoint
        if(!ignore_direction)
        {
          if(!(pairwise[startPoint+nObs*new_startPoint] == 0 ||
             pairwise[startPoint+nObs*new_startPoint] == 1
          )) continue;
        }
        else
        {
          if(!(pairwise[startPoint+nObs*new_startPoint] == -1 ||
             pairwise[startPoint+nObs*new_startPoint] == 0 ||
             pairwise[startPoint+nObs*new_startPoint] == 1
          )) continue;
        }

        std::vector<int> new_unused_nodes;
        for(int j=0; j<unused_nodes.size();j++)
        {
          if(j==i) continue;
          new_unused_nodes.push_back(unused_nodes[j]);
        }

        PathCount new_out = dfs_maxPathLength(pairwise,
                                              startPoint,
                                              new_startPoint,
                                              target,
                                              new_unused_nodes,
                                              ignore_direction,
                                              depth+1);

        // Rprintf("\n");
        // for(int i=0;i<depth;i++) Rprintf("\t ");
        // Rprintf("Next level down was %d hops with %d ties", new_out.number_hops, new_out.number_ties);

        // If a path was found, increment the number of hops required
        // to get there by 1
        if(new_out.number_hops>0){
          new_out.number_hops++;

          if(
             (pairwise[startPoint+nObs*new_startPoint] == 0)
          )
          {
            new_out.number_ties ++;
          }
        }

        // Only track the largest path we've seen

        if(new_out.number_hops > out.number_hops)
        {
          out.number_hops = new_out.number_hops;
          out.number_ties = new_out.number_ties;
        }
        else if(new_out.number_hops == out.number_hops &&
                new_out.number_ties < out.number_ties
                )
        {
          // If we're matched on the number of hops, store the
          // cycle with the least number of ties

          out.number_hops = new_out.number_hops;
          out.number_ties = new_out.number_ties;
        }
      }
    }

    // Rprintf("\n");
    // for(int i=0;i<depth;i++) Rprintf("\t");
    // Rprintf("Done in %d hops with %d ties", out.number_hops, out.number_ties);

    return out;
  }


  InconsistencyCounter(IntegerMatrix& pairwise, NumericVector& weight, bool includeTies, int m) :
    pairwise(pairwise),
    weight(weight),
    count(0),
    count_comparable(0),
    sum_weight(0),
    count_involved(pairwise.nrow(),0),
    includeTies(includeTies),
    m(m)
  {}

  InconsistencyCounter(InconsistencyCounter& counter, Split) :
    pairwise(counter.pairwise),
    weight(counter.weight),
    count(0),
    sum_weight(0),
    count_comparable(0),
    count_involved(pairwise.nrow(),0),
    includeTies(counter.includeTies),
    m(counter.m)
  {}

  void operator()(std::size_t begin, std::size_t end){

    int N = pairwise.nrow();

    for(int i=begin; i<end; i++)
    {

      // Rprintf("\n======================================\n\n");

      // Rprintf("N= %d, i=%d", N,i);

      // We have an m-tuple of (i, positions[0], ... positions[m-2])
      // where to prevent duplication we require each position to be larger
      // than the last.

      // Start with (i,i+1,i+2,...,i+m-1) and increment from the last position

      // roll-over by starting from the final position, scanning until we
      // stop finding illegal values, increment the allowed value
      // and reconstruct the rest of the m-tuplet as described above.

      std::vector<int> positions(m-1);
      positions[0] = i+1;
      for(int j=1; j<(m-1); j++)
      {
        positions[j] = positions[j-1]+1;
      }

      // Largest index values allowed in each position
      std::vector<int> OOB(m-1);
      OOB[m-2] = N-1;
      for(int j=m-3; j>=0; j--)
      {
        OOB[j] = OOB[j+1]-1;
      }

     // Rprintf("\n\n OOB: ");
     // for(int k=0;k<(m-1);k++) Rprintf("%d ",OOB[k]);
     // Rprintf("\n");


      int iter = 0;
      while(positions[0]<=OOB[0])
      {

        // Evaluate this m-tuplet

        // Rprintf("\n %d ", i);
        // for(int k=0;k<(m-1);k++) Rprintf("%d ",positions[k]);

        // DFS from i and attempt to construct a cycle


        // Count number of cycles
        PathCount maxCycleLength = this->dfs_maxPathLength(pairwise,
                                                           i, i, i,
                                                           positions,
                                                           false,
                                                           0);

        // Count number of cycles that can be formed assuming comparability
        PathCount maxCycleConnected = this->dfs_maxPathLength(pairwise,
                                                              i, i, i,
                                                              positions,
                                                              true,
                                                              0);

        // Rprintf("\t %d %d %d", maxCycleLength.number_hops, maxCycleLength.number_ties, maxCycleConnected.number_hops);


        double thisWeight = weight[i];
        for(int k=0;k<(m-1);k++) thisWeight *= weight[positions[k]];

        this->sum_weight += thisWeight;

        // At present only look at cycles of length m
        // i.e. don't consider sub-cycles

        if(maxCycleConnected.number_hops == m)
        {
          this->count_comparable += thisWeight;
        }

        if(maxCycleLength.number_hops==m & maxCycleLength.number_ties<m )
        {
          this->count += thisWeight;

          this->count_involved[i]++;
          for(int k=0;k<(m-1);k++) this->count_involved[positions[k]]++;
        }



        /* * * * * * * * * * * * * * * * *

         Increment to the next m-tuplet
         Start from back and identify which
         position to increment

         * * * * * * * * * * * * * * * * */


        int j=m-2; // Final position in vector
        while(positions[j]+1>OOB[j] && j>=0) j--;

        // Rprintf(" | Incrementing %d", j);

        // All positions are at their maximum value,
        // combinations have been exhausted.
        if(j<0)
        {
          // Rprintf(" j=%d ",j);
          break;
        }

        // The j'th position is incrementable and needs updated
        positions[j]++;

        // Reconstruct the rest of the m-tuplet post-increment.
        j++;
        while(j<(m-1))
        {
          positions[j] = positions[j-1]+1;
          j++;
        }

        // Failsafe for debugging
        // iter++;
        // if(iter > 200)
        // {
        //   Rprintf(" !!! ");
        //   break;
        // }


      }
    }
  }


  void join(const InconsistencyCounter& rhs)
  {
    count += rhs.count;
    sum_weight += rhs.sum_weight;
    count_comparable += rhs.count_comparable;

    int N = pairwise.nrow();
    for(int i=0; i<N; i++)
    {
      count_involved[i] += rhs.count_involved[i];
    }
  }

};



// [[Rcpp::export]]
List countInconsistency(IntegerMatrix x, NumericVector weights, bool includeTies, int m)
{

  InconsistencyCounter counter(x,weights, includeTies, m);

  int N = x.nrow();
  parallelReduce(0, N-m+1, counter);

  IntegerVector count_individual(N);
  for(int i=0; i<N; i++) count_individual(i) = counter.count_involved[i];

  List out = List::create(
                     Named("count") = wrap(counter.count),
                     Named("count_connected") = wrap(counter.count_comparable),
                     Named("totalWeight") = wrap(counter.sum_weight),
                     Named("individual") = wrap(count_individual)
  );

  return out;

}
