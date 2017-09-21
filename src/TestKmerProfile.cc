#include <string>
#include "ryggrad/src/base/CommandLineParser.h"
#include "ryggrad/src/base/FileParser.h"
#include "ryggrad/src/general/DNAVector.h"
#include "ryggrad/src/aligns/KmerAlignCore.h"


void Read(svec<double> & a, const string & fileName)
{
  FlatFileParser parser;
  
  parser.Open(fileName);

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    a.push_back(parser.AsFloat(0));
  }
}


int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandArg<string> r1Cmmd("-r1","input file");
  commandArg<string> r2Cmmd("-r2","input file");
  commandArg<string> r3Cmmd("-r3","input file");
  commandLineParser P(argc,argv);
  P.SetDescription("Testing the file parser.");
  P.registerArg(fileCmmd);
  P.registerArg(r1Cmmd);
  P.registerArg(r2Cmmd);
  P.registerArg(r3Cmmd);

  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);

  int k = 6;

  TranslateBasesToNumberExact trans;
  trans.SetSize(k);

  int n = trans.GetBoundValue();

  svec<double> r1, r2, r3;

  Read(r1, P.GetStringValueFor(r1Cmmd));
  Read(r2, P.GetStringValueFor(r2Cmmd));
  Read(r3, P.GetStringValueFor(r3Cmmd));

  vecDNAVector dna;
  
  dna.Read(fileName);



  int i, j;

  for (i=0; i<dna.isize(); i++) {
    const DNAVector & d = dna[i];
    svec<double> counts;
    counts.resize(n, 0.);
    double div = 0.;
    for (j=0; j<=d.isize()-k; j++) {
      int idx = trans.BasesToNumber(d, j);
      if (idx >= 0) {
	counts[idx] += 1;
	div++;
      }
    }
    for (j=0; j<n; j++)
      counts[j] /= div;


    double d1 = 0.;
    double d2 = 0.;
    double d3 = 0.;

    for (j=0; j<n; j++) {
      d1 += (counts[j] - r1[j])*(counts[j] - r1[j]);
      d2 += (counts[j] - r2[j])*(counts[j] - r2[j]);
      d3 += (counts[j] - r3[j])*(counts[j] - r3[j]);
    }

    //cout << dna.Name(i) << "\t" << d1 <<  "\t" << d2 <<  "\t" << d3 << endl;
    cout << dna.Name(i) << "\t" << d1 <<  "\t" << d2 << "\t" << d.isize() << "\t";
    cout << 1000*(d2 - d1) << "\t";
    if (d1 < d2) 
      cout << "nuclear";
    else
      cout << "chloroplast";
      
    cout << endl;

  }


  return 0;
}
