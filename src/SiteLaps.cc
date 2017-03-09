//#define FORCE_DEBUG
#include "OptiMapAlignUnit.h"

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandArg<int> kCmmd("-k","seed size", 6);
  commandArg<int> wCmmd("-w","wiggle", 5);
  commandLineParser P(argc,argv);
  P.SetDescription("Find overlaps in restriction maps.");
  P.registerArg(fileCmmd);
  P.registerArg(kCmmd);
  P.registerArg(wCmmd);
 
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  int k = P.GetIntValueFor(kCmmd);
  int w = P.GetIntValueFor(wCmmd);
  OptiMapAlignUnit omaUnit; 
  // 1. Load Optical Reads
  omaUnit.LoadReads(fileName, k);
  // 2. Build Optimers and find those that share a seed as cadidates for overlap detection 
  omaUnit.FindCandidLaps(k);
  
  return 0;
}
