#define FORCE_DEBUG

#include <string>
#include "ryggrad/src/base/CommandLineParser.h"
#include "ryggrad/src/base/FileParser.h"
#include "ryggrad/src/general/DNAVector.h"
#include "ryggrad/src/aligns/KmerAlignCore.h"

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandLineParser P(argc,argv);
  P.SetDescription("Testing the file parser.");
  P.registerArg(fileCmmd);

  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);

  vecDNAVector dna;
  
  dna.Read(fileName);

  int k = 6;

  TranslateBasesToNumberExact trans;
  trans.SetSize(k);

  int n = trans.GetBoundValue();

  svec<int> counts;
  double div = 0.;
  counts.resize(n, 0);

  int i, j;

  for (i=0; i<dna.isize(); i++) {
    const DNAVector & d = dna[i];
    for (j=0; j<=d.isize()-k; j++) {
      int idx = trans.BasesToNumber(d, j);
      if (idx >= 0) {
	counts[idx]++;
	div++;
      }
    }
  }

  for (i=0; i<counts.isize(); i++)
    cout << ((double)counts[i])/div << endl;

  return 0;
}
