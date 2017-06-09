#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include <ctime>
#include "OptiMapAlignUnit.h"


int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input fasta file");
  commandArg<int> readCountCmmd("-c","number of reads in input fasta file");
  commandArg<int> kCmmd("-k","seed size", 6);
  commandArg<int>  wCmmd("-w","Wiggle i.e. Error Tolerance in finding overlaps", 5);
  commandArg<int>  coreCmmd("-n","Number of Cores to run with", 2);
  commandLineParser P(argc,argv);
  P.SetDescription("Find overlaps in restriction maps.");
  P.registerArg(fileCmmd);
  P.registerArg(readCountCmmd);
  P.registerArg(kCmmd);
  P.registerArg(wCmmd);
  P.registerArg(coreCmmd);
 
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  int readCnt     = P.GetIntValueFor(readCountCmmd);
  int seedSize    = P.GetIntValueFor(kCmmd);
  int wiggle      = P.GetIntValueFor(wCmmd);
  int numOfCores  = P.GetIntValueFor(coreCmmd);

  OptiMapAlignUnit omaUnit; 
  clock_t clock1_optiLoad, clock2_overlapCand, clock3_finalOverlaps, clock4_done;

  omp_set_num_threads(numOfCores); //The sort functions use OpenMP

  // 1a. Populate the motifs 
  omaUnit.GenerateMotifs(4, 100);
  // 1b. Construct Restriction-site Reads
  clock1_optiLoad = clock();
  omaUnit.MakeRSites(fileName, readCnt);

  // 2. Build Optimers and find those that share a seed as cadidates for overlap detection 
  clock2_overlapCand = clock();
  OverlapCandids lapCandids;
//  omaUnit.FindLapCandids(seedSize, lapCandids);
 
  // 3. Take the overlap candidates and refine to remove false positives
  clock3_finalOverlaps = clock();
  OverlapCandids finalOverlaps;
//  omaUnit.FinalOverlaps(lapCandids, wiggle, finalOverlaps);
  clock4_done = clock();


  cout << "Report runtime duration: " << endl;
  cout << " Loading Optical Reads: "  
       << ((double) (clock2_overlapCand-clock1_optiLoad) / CLOCKS_PER_SEC)
       << endl;
  cout << " Finding Overlap Candidates: " 
       << ((double) (clock3_finalOverlaps-clock2_overlapCand) / CLOCKS_PER_SEC)
       << endl;
  cout << " Refining and finalizing overlaps: " 
       << ((double) (clock4_done - clock3_finalOverlaps) / CLOCKS_PER_SEC)
       << endl;

  return 0;
}
