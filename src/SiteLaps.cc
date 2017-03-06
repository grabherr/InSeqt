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
  // 2. Build optimers from opti reads with the given seed size
  omaUnit.BuildOptimers(k);
  // 3. Create Links for linking optimers
  //ORLinks links(optiReads.NumReads());
  // 4. Pool linked optimers together
  omaUnit.PoolReads();
  
  return 0;
}
