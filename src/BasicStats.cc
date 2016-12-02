#include <string>
#include "ryggrad/src/base/CommandLineParser.h"
#include "ryggrad/src/general/DNAVector.h"
#include "ryggrad/src/visual/Histogram.h"
#include "ryggrad/src/base/FileParser.h"
#include "DataRead.h"

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
  cout << "Calling " << cmmd << endl;
  int r = system(cmmd.c_str());

  if (r != 0) {
    cout << "ERROR: imagemagick could not convert image, returned: " << r << endl;
    exit(-1);
  }
  return true;
}

string GetShort(const string & full) {
  char tmp[2028];
  strcpy(tmp, full.c_str());

  int i;
  int k = 0;
  bool b = false;
  for (i=0; i<(int)full.size(); i++) {
    if (tmp[i] == '_') {
      k++;
      if (k == 3) {
	tmp[i] = 0;
	b = true;
	break;
      }
    }
  }
  string ret = tmp;
  if (!b)
    ret = "";
  return ret;
}

bool CheckIM()
{
  string cmmd = "convert > /dev/null ";
  int r = system(cmmd.c_str());
  if (r != 0 && r != 256) {
    cout << "ImageMagick convert not found, returned " << r << " Exiting." << endl;
    exit(-1);
  }
  cout << "ImageMagick installed, OK!" << endl;
  return true;
}
 
int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandArg<string> outCmmd("-o","output directory");
  
   
  commandLineParser P(argc,argv);
  P.SetDescription("Compute basic stats for PacBio reads.");
  P.registerArg(fileCmmd);
  P.registerArg(outCmmd);

  
  
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string outName = P.GetStringValueFor(outCmmd);
  CheckIM();

 
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
  svec<double> lengths_polym;

  string last;
  int same = 0;
  
  svec<int> cycles;
  double  total = 0;

  string summary = outName;
  summary += "/basicstats.out";
  FILE * pRep = fopen(summary.c_str(), "w");
  string cyclesummary = outName;
  cyclesummary += "/by_cycles.out";
  FILE * pCycle = fopen(cyclesummary.c_str(), "w");

  string chim = outName;
  chim += "/chimeras.out";
  FILE * pChim = fopen(chim.c_str(), "w");
  string trunc = outName;
  trunc += "/truncated.out";
  FILE * pTrunc = fopen(trunc.c_str(), "w");

  
  string lastName;
  int localReads = 0;
  double localTotal = 0;
  svec<int> localSize;

  svec<string> namecache;
  svec<int> lencache;

  int lastLen;
  int totalchim;
  int totaltrunc;
  
  for (i=0; i<dna.isize(); i++) {

    StringParser pp;
    pp.SetLine(dna.Name(i), "/");
    bool bCont = false;
    if (pp.GetItemCount() > 2) {
      const string & id = pp.AsString(1);
      if (id == last) {
	bCont = true;
	same++;
	namecache.push_back(dna.Name(i));
	lencache.push_back(dna[i].isize());
      } else {
	if (same >= cycles.isize())
	  cycles.resize(same+1, 0);
	cycles[same]++;
	int x;
	for (x=0; x<namecache.isize(); x++)
	  fprintf(pCycle, "%s\t%d\n", namecache[x].c_str(), same+1);

	bool bChim = false;
	for (x=1; x<lencache.isize(); x++) {
	  double rr = (double)lencache[x]/(double)lencache[x-1];
	  if (rr < 0.6) {
	    if (x == lencache.isize()-1) {
	      totaltrunc++;
	      fprintf(pTrunc, "%s\%d\n", namecache[x].c_str(), lencache[x]);
	    } else {
	      bChim = true;
	    }
	  }
	  if (rr > 1./0.6) {
	    bChim = true;
	  }
	}
	if (bChim) {
	  totalchim++;
	  for (x=0; x<namecache.isize(); x++)
	    fprintf(pChim, "%s\t%d\n", namecache[x].c_str(), same+1);
	}
	
	same = 0;
	last = id;
	namecache.clear();
	lencache.clear();
	namecache.push_back(dna.Name(i));
  	lencache.push_back(dna[i].isize());
    }
    }

    string shortName = GetShort(dna.Name(i));

    if (localReads > 0 && localSize.isize() > 0 && (shortName != lastName || i == dna.isize()-1)) {
      
      Sort(localSize);
      int localMedian = localSize[localSize.isize()/2];
      
      double local_nn = 0;
      double local_n50 = 0;
      for (j=0; j<localSize.isize(); j++) {
	local_nn += localSize[j];
	if (local_nn >= localTotal/2) {
	  local_n50 = localSize[j];
	  break;
	}
      }


      fprintf(pRep, "SMRT %s\n", lastName.c_str());
      fprintf(pRep, "Reads: %d\n", localReads);
      fprintf(pRep, "Average: %f\n", (double)localTotal/(double)localReads);
      fprintf(pRep, "Median: %d\n", localMedian);
      fprintf(pRep, "Total: %f\n", localTotal);
      fprintf(pRep, "N50: %f\n", local_n50);
      //lastName = shortName;
      localReads = 0;
      localTotal = 0;
      localSize.clear();
    }
    lastName = shortName;

 
    int n = 0;
    int gc = 0;
    int n_w = 0;
    int gc_w = 0;
    const DNAVector & d = dna[i];


    localReads++;
    localTotal += d.isize();
    localSize.push_back(d.isize());

    size.push_back(d.isize());
    total += d.isize();
    avg += d.isize();
    lengths.push_back(d.isize());

    if (!bCont) {
      lengths_polym.push_back(d.isize());
    } else {
      lengths_polym[lengths_polym.isize()-1] += d.isize();
    }

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

  for (int x=0; x<namecache.isize(); x++)
    fprintf(pCycle, "%s\t%d\n", namecache[x].c_str(), same+1);


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

  //cout << "Basic stats" << endl;
  cout << "smrt TOTAL " << endl;
  cout << "Reads: " << dna.isize() << endl;
  cout << "Average: " << avg << endl;
  cout << "Median: " << median << endl;
  cout << "Total: " << total << endl;
  cout << "N50: " << n50 << endl;

  fprintf(pRep, "SMRT total\n");
  fprintf(pRep, "Reads: %d\n", dna.isize());
  fprintf(pRep, "Average: %f\n", avg);
  fprintf(pRep, "Median: %d\n", median);
  fprintf(pRep, "Total: %f\n", total);
  fprintf(pRep, "N50: %f\n", n50);
  fprintf(pRep, "Chimeric: %d\n", totalchim);
  fprintf(pRep, "Truncated: %d\n", totaltrunc);

  cout << "Cycles." << endl;
  int max = 19;
  if (max > cycles.isize())
    max = cycles.isize();
  
  for (i=0; i<max; i++) {
    //if (cycles[i] > 0) {
    cout << "# " << i+1 << " " << cycles[i] << endl;
    fprintf(pRep, "# %d %d\n", i+1, cycles[i]);
      //}
  }
  int n20above = 0;
  for (i=max; i<cycles.isize(); i++) {
    n20above += cycles[i];
  }
  fprintf(pRep, "# >=%d %d\n", max+1, n20above);
  cout << endl;
  cout << "Wrote chimeras to " << chim << endl;
  cout << "Wrote truncated reads to " << trunc << endl;
  
  
  cout << "Done." << endl;
  fclose(pRep);
  fclose(pCycle);
  fclose(pChim);
  fclose(pTrunc);
  

  Histogram hh;

  string seq_n = outName;
  seq_n += "/seq.ps";
  string win_n = outName;
  win_n += "/win.ps";
  string size_n = outName;
  size_n += "/size.ps";
  string sizepolym_n = outName;
  sizepolym_n += "/sizepolym.ps";
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
  cout << "Plot Polym Size " << endl;
  hh.Plot(sizepolym_n, lengths_polym, 100, color(0.99, 0.2, 0.2));


  //MakePNG(seq_n);
  //MakePNG(win_n);
  //MakePNG(size_n);
  //MakePNG(sizepolym_n);
  //MakePNG(scatter_n);

  

  return 0;
}
