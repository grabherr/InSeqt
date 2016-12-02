#include <string>
#include "ryggrad/src/base/CommandLineParser.h"
#include "ryggrad/src/general/DNAVector.h"
#include "ryggrad/src/base/FileParser.h"
#include "src/DataRead.h"
 
int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input list of read names");
  commandArg<string> fastCmmd("-f","input fastq or fasta file");
  commandArg<string> outCmmd("-o","output file");
  commandArg<bool> rCmmd("-r","remove these sequences, don't keep", false);
  
   
  commandLineParser P(argc,argv);
  P.SetDescription("Extracts a subset of reads give a list of read names.");
  P.registerArg(fileCmmd);
  P.registerArg(fastCmmd);
  P.registerArg(outCmmd);
  P.registerArg(rCmmd);

  
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string fastName = P.GetStringValueFor(fastCmmd);
  string outName = P.GetStringValueFor(outCmmd);
  bool bRem = P.GetBoolValueFor(rCmmd);

 

  vecDNAVector dna;
  vecDNAVector out;
  ReadDNA(dna, fastName);
  
  int i, j;
  FlatFileParser parser;
  
  parser.Open(fileName);
  svec<string> s;
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    s.push_back(parser.Line());
  }
  Sort(s);

  for (i=0; i<dna.isize(); i++) {
    const string & n = dna.Name(i);
    int index = BinSearch(s, n);
    if (bRem) {
      if (index >= 0)  {
	continue;
      } else {
	out.push_back(dna[i], dna.Name(i));
      }
    } else {
      if (index >= 0)  {
	out.push_back(dna[i], dna.Name(i));
      } else {
 	continue;
      }
    }
  }
    
    out.Write(outName);
  return 0;
}
