#include <string>
#include "base/CommandLineParser.h"
#include "src/DNAVector.h"
#include "visual/Histogram.h"


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
  
  int i, j;

  
  svec<double> win;
  int winsize = 1000;
  svec<double> seq;
  svec<double> size;

  for (i=0; i<dna.isize(); i++) {
    int n = 0;
    int gc = 0;
    int n_w = 0;
    int gc_w = 0;
    const DNAVector & d = dna[i];
    size.push_back(d.isize());

    for (j=0; j<d.isize(); j++) {
      if (d[j] == 'A' || d[j] == 'C' || d[j] == 'G' || d[j] == 'T') {
	n++;
	n_w++;
	if (d[j] == 'C' || d[j] == 'G') {
	  gc++;
	  gc_w++;
	}
	if (n_w == winsize) {
	  win.push_back((double)gc_w/(double)n_w);
	  n_w = 0;
	  gc_w = 0;
	}
      }
    }
    if (n > 50)
      seq.push_back((double)gc/(double)n);
  

  }


  Histogram hh;
  hh.Plot("seq.ps", seq, 100, 0., 1.);
  hh.Plot("win.ps", win, 100, 0., 1.);
  hh.Plot("size.ps", size, 100, color(0.99, 0.2, 0.2));
  hh.Scatter("scatter.ps", seq, size, color(0.3, 0.0, 0.4));


  return 0;
}
