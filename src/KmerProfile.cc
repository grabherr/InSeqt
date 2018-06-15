#define FORCE_DEBUG
#include <string>
#include "ryggrad/src/base/CommandLineParser.h"
#include "ryggrad/src/base/FileParser.h"
#include "ryggrad/src/base/FileParser.h"
#include "ryggrad/src/general/DNAVector.h"
#include "ryggrad/src/aligns/KmerAlignCore.h"
#include <math.h>

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-f","input fasta file");
  commandArg<int> kCmmd("-k","word size");
  commandLineParser P(argc,argv);
  P.SetDescription("Computes a k-mer profile.");
  P.registerArg(fileCmmd);
  P.registerArg(kCmmd);
 
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  int k = P.GetIntValueFor(kCmmd);
  int i, j;

  TranslateBasesToNumberExact trans;
  trans.SetSize(k); // Use 15-mers
  svec<double> freq;
  int n = trans.GetBoundValue();

  freq.resize(n, 0.);

  
  vecDNAVector dna;
  dna.Read(fileName);

  double total = 0;

  for (i=0; i<dna.isize(); i++) {
    const DNAVector & d = dna[i];
    for (j=0; j<=d.isize()-k; j++) {
      int index = trans.BasesToNumber(d, j);
      if (index >= 0) {
	freq[index] += 1.;
	total += 1.;
      }
    }
  }
  cout << k << endl;
  double len = 0.;
  for (i=0; i<freq.isize(); i++)
    len += freq[i]*freq[i];
  len = sqrt(len);
  
  for (i=0; i<freq.isize(); i++)
    cout << freq[i]/len << endl;
  
  return 0;
}
