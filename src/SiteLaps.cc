#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include <ctime>
#include "RestSiteAlignUnit.h"


int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input fasta file");
  commandArg<int> readCountCmmd("-c","number of reads in input fasta file");
  commandArg<int> kCmmd("-k","seed size", 6);
  commandArg<int> motifCmmd("-m","Motif Length", 4);
  commandArg<int>  wCmmd("-w","Wiggle i.e. Error Tolerance in finding overlaps", 5);
  commandArg<int>  coreCmmd("-n","Number of Cores to run with", 2);
  commandArg<string> appLogCmmd("-L","Application logging file","application.log");
  commandLineParser P(argc,argv);
  P.SetDescription("Find overlaps in restriction maps.");
  P.registerArg(fileCmmd);
  P.registerArg(readCountCmmd);
  P.registerArg(kCmmd);
  P.registerArg(motifCmmd);
  P.registerArg(wCmmd);
  P.registerArg(coreCmmd);
 
  P.parse();
  
  string fileName  = P.GetStringValueFor(fileCmmd);
  int readCnt      = P.GetIntValueFor(readCountCmmd);
  int seedSize     = P.GetIntValueFor(kCmmd);
  int motifLen     = P.GetIntValueFor(motifCmmd);
  int wiggle       = P.GetIntValueFor(wCmmd);
  int numOfCores   = P.GetIntValueFor(coreCmmd);
    string logFile = P.GetStringValueFor(appLogCmmd);

  FILE* pFile               = fopen(logFile.c_str(), "w");
  Output2FILE::Stream()     = pFile;
  FILELog::ReportingLevel() = logINFO; 
  FILELog::ReportingLevel() = logDEBUG4; 
#if defined(FORCE_DEBUG)
//    FILELog::ReportingLevel() = logDEBUG4; 
#endif

  RestSiteMapper rsMapper;
  clock_t clock1_optiLoad, clock2_overlapCand, clock3_finalOverlaps, clock4_done;

  omp_set_num_threads(numOfCores); //The sort functions use OpenMP

  // 1a. Populate the motifs 
  rsMapper.GenerateMotifs(motifLen, 100);
  // 1b. Construct Restriction-site Reads
  clock1_optiLoad = clock();

  // 2. Build Optimers and find those that share a seed as cadidates for overlap detection 
  MatchCandids finalOverlaps;
  rsMapper.FindMatches(fileName, readCnt, 0, finalOverlaps); 
  clock2_overlapCand = clock();
 
  // 3. Take the overlap candidates and refine to remove false positives
  clock3_finalOverlaps = clock();
//  rsMapper.FinalOverlaps(lapCandids, wiggle, finalOverlaps);
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
