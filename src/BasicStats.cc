#include <string>
#include "base/CommandLineParser.h"
#include "src/DNAVector.h"
#include "visual/Histogram.h"
#include "base/FileParser.h"
#include "src/DataRead.h"

bool MakePNG(const string & fileName)
{
  char tmp[512];
  strcpy(tmp, fileName.c_str());
  tmp[strlen(tmp)-2] = 0;
  strcat(tmp, "png");
  string cmmd = "convert ";
  cmmd += fileName;
  cmmd += " ";
  cmmd += tmp;
  int r = system(cmmd.c_str());

  if (r != 0) {
    cout << "ERROR: imagemagick is not installed!!" << endl;
    exit(-1);
  }
  return true;
}


 
int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandArg<string> outCmmd("-o","output directory");
  //commandArg<bool> listCmmd("-f","", "");
  commandLineParser P(argc,argv);
  P.SetDescription("Testing the file parser.");
  P.registerArg(fileCmmd);
  P.registerArg(outCmmd);
 
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string outName = P.GetStringValueFor(outCmmd);
 
 
  string makeOut = "mkdir ";
  makeOut += outName;
  int rr = system(makeOut.c_str());

  vecDNAVector dna;
  ReadDNA(dna, fileName);
  
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
  double  total = 0;

  for (i=0; i<dna.isize(); i++) {

    StringParser pp;
    pp.SetLine(dna.Name(i), "/");
    if (pp.GetItemCount() > 2) {
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
    }

    int n = 0;
    int gc = 0;
    int n_w = 0;
    int gc_w = 0;
    const DNAVector & d = dna[i];
    size.push_back(d.isize());
    total += d.isize();
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
    //if (n > 50)
    seq.push_back((double)gc/(double)n);
  

  }
  avg /= (double)dna.isize();
  Sort(lengths);
  int median = lengths[lengths.isize()/2];

  double nn = 0;
  double n50 = 0;
  for (i=0; i<lengths.isize(); i++) {
    nn += lengths[i];
    if (nn >= total/2) {
      n50 = lengths[i];
      break;
    }
  }

  string summary = outName;
  summary += "/basicstats.out";
  FILE * pRep = fopen(summary.c_str(), "w");
  //cout << "Basic stats" << endl;
  cout << "Reads: " << dna.isize() << endl;
  cout << "Average: " << avg << endl;
  cout << "Median: " << median << endl;
  cout << "Total: " << total << endl;
  cout << "N50: " << n50 << endl;

  fprintf(pRep, "Reads: %d\n", dna.isize());
  fprintf(pRep, "Average: %f\n", avg);
  fprintf(pRep, "Median: %d\n", median);
  fprintf(pRep, "Total: %f\n", total);
  fprintf(pRep, "N50: %f\n", n50);


  for (i=0; i<cycles.isize(); i++) {
    if (cycles[i] > 0) {
      cout << "# " << i+1 << " " << cycles[i] << endl;
      fprintf(pRep, "# %d %d\n", i+1, cycles[i]);
    }
  }
  fclose(pRep);

  Histogram hh;

  string seq_n = outName;
  seq_n += "/seq.ps";
  string win_n = outName;
  win_n += "/win.ps";
  string size_n = outName;
  size_n += "/size.ps";
  string scatter_n = outName;
  scatter_n += "/scatter.ps";

  cout << "Plot Seq " << endl;
  hh.Plot(seq_n, seq, 100, 0., 1.);
  cout << "Plot Win " << endl;
  hh.Plot(win_n, win, 100, 0., 1.);
  cout << "Plot Size " << endl;
  hh.Plot(size_n, size, 100, color(0.99, 0.2, 0.2));
  cout << "Scatter " << seq.isize() << " " << size.isize() << endl;
  hh.Scatter(scatter_n, seq, size, color(0.3, 0.0, 0.4));

  MakePNG(seq_n);
  MakePNG(win_n);
  MakePNG(size_n);
  MakePNG(scatter_n);

  

  return 0;
}
