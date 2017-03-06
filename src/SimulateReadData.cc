#define FORCE_DEBUG
#include <string>
#include "ryggrad/src/base/CommandLineParser.h"
#include "ryggrad/src/base/FileParser.h"
#include "ryggrad/src/general/DNAVector.h"
#include "ryggrad/src/base/RandomStuff.h"

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandArg<int> bCmmd("-l","average length");
  commandArg<double> eCmmd("-c","coverage");
  commandLineParser P(argc,argv);
  P.SetDescription("Testing the file parser.");
  P.registerArg(fileCmmd);
  P.registerArg(bCmmd);
  P.registerArg(eCmmd);
 
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  int l = P.GetIntValueFor(bCmmd);
  double cov = P.GetDoubleValueFor(eCmmd);
 
  vecDNAVector dna;
  dna.Read(fileName);

  svec<int> len;
  len.resize(dna.isize(), 0);

  int i, j;

  int total = 0;
  for (i=0; i<dna.isize(); i++) {
    total += dna[i].isize();
    len[i] = total;
    //cout << len[i] << endl;
  }

  double wanted = (double)total*cov;
  double have = 0.;
  int k = 0;
  while (have < wanted) {
    int s = RandomInt(total);
    for (j=0; j<len.isize(); j++) {
      //cout << s << "   " << len[j] << endl;
      if (s > len[j])
	continue;
      //j--;
      break;
    }
    if (j > 0)
      s -= len[j-1];
    int howmany = l + RandomInt(l/2)-l/4;
    const DNAVector & d = dna[j];
    cout << dna.Name(j) << "_" << k;
    k++;
    //cout << "  " << s << " " << j << " " << len[j] << " " << d.isize() << " " << howmany << endl;
    if (RandomFloat(1.) < 0.5) {
      cout << "+" << endl;
      for (j=s; j<s+howmany; j++) {
	if (j >= d.isize())
	  break;
	cout << d[j];
	have += 1.;
	if (j % 80 == 0)
	  cout << endl;
      }
      cout << endl;
    } else {
      cout << "-" << endl;
      for (j=s; j>s-howmany; j--) {
	if (j < 0)
	  break;
	//cout << d[j];
	if (d[j] == 'A')
	  cout << "T";
	if (d[j] == 'C')
	  cout << "G";
	if (d[j] == 'G')
	  cout << "C";
	if (d[j] == 'T')
	  cout << "A";
	if (d[j] == 'N')
	  cout << "N";
  	have += 1.;
 	if (j % 80 == 0)
	  cout << endl;
      }
      cout << endl;
    }
  }
  

  return 0;
}
