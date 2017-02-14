#define FORCE_DEBUG
#include <string>
#include "ryggrad/src/base/CommandLineParser.h"
#include "ryggrad/src/base/FileParser.h"
#include "ryggrad/src/base/FileParser.h"
#include "ryggrad/src/general/DNAVector.h"
#include "ryggrad/src/aligns/KmerAlignCore.h"
#include <math.h>

double Distance(const svec<double> & a, const svec<double> & b)
{
  double d = 0.;
  int i;
  for (i=0; i<a.isize(); i++) {
    //cout << a[i] << " " << b[i] << endl;
    //d += (a[i]-b[i])*(a[i]-b[i]);
    d += a[i]*b[i];
  }

  return sqrt(d);
}

int main( int argc, char** argv )
{

  commandArg<string> fCmmd("-f","input file");
  commandArg<string> fileCmmd("-i","k-mer profile (run KmerProfile to generate)");
  commandLineParser P(argc,argv);
  P.SetDescription("Scores against a k-mer profile.");
  P.registerArg(fileCmmd);
  P.registerArg(fCmmd);
 
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string fName = P.GetStringValueFor(fCmmd);

  vecDNAVector dna;
  dna.Read(fName);
  
  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);

 
  int i, j;

  parser.ParseLine();
  int k = parser.AsInt(0);

  TranslateBasesToNumberExact trans;
  trans.SetSize(k); // Use 15-mers

  int n = trans.GetBoundValue();
  svec<double> freq;
  freq.resize(n, 0.);

  i = 0;
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    freq[i] = parser.AsFloat(0);
    i++;
  }

  
  for (i=0; i<dna.isize(); i++) {
    const DNAVector & d = dna[i];
    //cout << "Check " << dna.Name(i) << endl;
    svec<double> tmp;
    tmp.resize(freq.isize(), 0.);
    double total = 0.;
    for (j=0; j<=d.isize()-k; j++) {
      int index = trans.BasesToNumber(d, j);
      if (index >= 0) {
	tmp[index] += 1.;
	total += 1.;
      }
    }
    double len = 0.;
    for (j=0; j<freq.isize(); j++)
      len += tmp[j]*tmp[j];
    len = sqrt(len);
    for (j=0; j<tmp.isize(); j++)
      tmp[j] /= len;
    
    cout << dna.Name(i) << " " << Distance(tmp, freq) << endl;
  }
 
  return 0;
}
