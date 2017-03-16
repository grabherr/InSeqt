//#define FORCE_DEBUG
#include "OptiMapAlignUnit.h"

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandArg<int> kCmmd("-k","seed size", 6);
  commandArg<int>  wCmmd("-w","Wiggle i.e. Error Tolerance in finding overlaps", 5);
  commandArg<int>  coreCmmd("-n","Number of Cores to run with", 2);
  commandLineParser P(argc,argv);
  P.SetDescription("Find overlaps in restriction maps.");
  P.registerArg(fileCmmd);
  P.registerArg(kCmmd);
  P.registerArg(wCmmd);
  P.registerArg(coreCmmd);
 
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  int seedSize    = P.GetIntValueFor(kCmmd);
  int wiggle      = P.GetIntValueFor(wCmmd);
  int numOfCores  = P.GetIntValueFor(coreCmmd);

  OptiMapAlignUnit omaUnit; 

  omp_set_num_threads(numOfCores); //The sort functions use OpenMP

  // 1. Load Optical Reads
  omaUnit.LoadReads(fileName, seedSize);

  // 2. Build Optimers and find those that share a seed as cadidates for overlap detection 
  OverlapCandids lapCandids;
  omaUnit.FindLapCandids(seedSize, lapCandids);
 
  // 3. Take the overlap candidates and refine to remove false positives
  OverlapCandids finalOverlaps;
  omaUnit.FinalOverlaps(lapCandids, wiggle, finalOverlaps);

  return 0;
}
