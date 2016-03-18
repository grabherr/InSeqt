#include <string>
#include "base/CommandLineParser.h"
#include "src/DNAVector.h"
#include "visual/Histogram.h"
#include "base/FileParser.h"


void ReadFileNames(string & fileName)
{
  FlatFileParser parser;
  
  parser.Open(fileName);

  fileName = "";

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    for (int i=0; i<parser.GetItemCount(); i++) {
      if (fileName != "")
	fileName += ",";
      fileName += parser.AsString(i);
    }
  }
}
 
int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  //commandArg<bool> listCmmd("-f","", "");
  commandLineParser P(argc,argv);
  P.SetDescription("Testing the file parser.");
  P.registerArg(fileCmmd);
 
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
 
  FlatFileParser parser;
  
  parser.Open(fileName);
  parser.ParseLine();

  if (parser.AsString(0)[0] != '@') {
    ReadFileNames(fileName);
  }


  vecDNAVector dna;
  dna.Read(fileName);
  
  int i, j;

  
  svec<double> win;
  int winsize = 1000;
  svec<double> seq;
  svec<double> size;

  double avg = 0.;
  svec<int> lengths;

  string last;
  int same = 0;
  
  svec<int> cycles;

  for (i=0; i<dna.isize(); i++) {

    StringParser pp;
    pp.SetLine(dna.Name(i), "/");
    const string & id = pp.AsString(1);
    if (id == last) {
      same++;
    } else {
      if (same >= cycles.isize())
	cycles.resize(same+1, 0);
      cycles[same]++;
      same = 0;
      last = id;
    }


    int n = 0;
    int gc = 0;
    int n_w = 0;
    int gc_w = 0;
    const DNAVector & d = dna[i];
    size.push_back(d.isize());
    avg += d.isize();
    lengths.push_back(d.isize());

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
  avg /= (double)dna.isize();
  Sort(lengths);
  int median = lengths[lengths.isize()/2];

  cout << "Basic stats" << endl;
  cout << "Reads: " << dna.isize() << endl;
  cout << "Average: " << avg << endl;
  cout << "Median: " << median << endl;

  for (i=0; i<cycles.isize(); i++) {
    if (cycles[i] > 0)
      cout << "# " << i+1 << " " << cycles[i] << endl;
  }


  Histogram hh;
  hh.Plot("seq.ps", seq, 100, 0., 1.);
  hh.Plot("win.ps", win, 100, 0., 1.);
  hh.Plot("size.ps", size, 100, color(0.99, 0.2, 0.2));
  hh.Scatter("scatter.ps", seq, size, color(0.3, 0.0, 0.4));


  return 0;
}
