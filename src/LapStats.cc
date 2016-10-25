#include <string>
#include "ryggrad/src/base/CommandLineParser.h"
#include "ryggrad/src/visual/Histogram.h"
#include "ryggrad/src/base/FileParser.h"

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
    cout << "ERROR: imagemagick is not installed!!" << endl;
    exit(-1);
  }
  return true;
}

void LoadReads(svec<string> & reads, const string & fileName)
{
  FlatFileParser parser;
  
  parser.Open(fileName);
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    reads.push_back(parser.AsString(0));
  }
}
 
int main( int argc, char** argv )
{

  //commandArg<string> fileCmmd("-i","input overlap file");
  //commandArg<string> readCmmd("-r","input read file");
  commandArg<string> outCmmd("-i","data directory");
  //commandArg<bool> listCmmd("-f","", "");
  commandLineParser P(argc,argv);
  P.SetDescription("Computes statistics from overlaps.");
  //P.registerArg(fileCmmd);
  //P.registerArg(readCmmd);
  P.registerArg(outCmmd);
 
  P.parse();
  
  string fileName = P.GetStringValueFor(outCmmd);
  fileName += "/overlapcands.out";
  string readName = fileName + ".0.allreadnames.0";
  string outName = P.GetStringValueFor(outCmmd);
 
 
  //string makeOut = "mkdir ";
  //makeOut += outName;
  //int rr = system(makeOut.c_str());

  svec<string> reads;
  LoadReads(reads, readName);

  Sort(reads);

  svec<double> laps;
  laps.resize(reads.isize(), 0);
  svec<double> raw_len;
  raw_len.resize(reads.isize(), 0);
  svec<double> lap_len;
  lap_len.resize(reads.isize(), 0);
 
  FlatFileParser parser;
  
  parser.Open(fileName);

  //parser.ParseLine();
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    const string & s = parser.AsString(0);
    if (s == "Seq1")
      continue;
    int len = parser.AsInt(4);
    int len2 = parser.AsInt(5);

    //if (len < 10000)
    //   continue;
     

    int index = BinSearch(reads, s);
    if (index < 0)
      cout << "ERROR, read not found: " << s << endl;

    raw_len[index] = len;
    lap_len[index] += len2;
 
    laps[index] += 1.;
  }
 
  svec<double> laps_clean;

  svec<double> size;

  double totalLen = 0.;
  double totalLaps = 0.;
  double total = 0.;
  svec<int> for_med;
  for (int j=0; j<laps.isize(); j++) {
    if (laps[j] > 0 && laps[j] < 200) {
      total += laps[j];
      laps_clean.push_back(laps[j]);
      size.push_back(raw_len[j]);
      totalLaps += lap_len[j];
      totalLen += raw_len[j];
      for_med.push_back(laps[j]);
    }
  }

  Sort(for_med);
  double median = for_med[for_med.isize()/2];
  
  double cov = totalLaps / totalLen;
  cout << "Total length:  " << totalLen << endl;
  cout << "Total reads:   " << laps.isize() << endl;
  cout << "Passing reads: " << laps_clean.isize() << endl;
  cout << "Base coverage: " << cov << endl;
  cout << "Genome size:   " << totalLen / cov / 1000000 <<  " MB" << endl;
  cout << "NOTE: multiply the genome size estimate according to sub-sampling!" << endl; 

  double frac;
  FlatFileParser parser2;
  string fracName = P.GetStringValueFor(outCmmd);
  fracName += "/overlapcands.txt";
 
  parser2.Open(fracName);

  parser2.ParseLine();
  frac = parser2.AsFloat(1);
 

  string statsName = P.GetStringValueFor(outCmmd);
  statsName += "/overlapstats.txt";
  FILE * pStats = fopen(statsName.c_str(), "w");
  
 
  fprintf(pStats, "tested %d\n", laps.isize());
  fprintf(pStats, "testedper %f\n", 100*frac);
  fprintf(pStats, "perread %f\n", (double)total/(double)laps_clean.isize());
  fprintf(pStats, "perreadmed %f\n", median);
  fprintf(pStats, "passing %d\n", laps_clean.isize());
  fprintf(pStats, "total %f\n", total);
  fprintf(pStats, "totallength %f\n", totalLen);
  fprintf(pStats, "lapsplus 0\n");
  fprintf(pStats, "lapsminus 0\n");
  fprintf(pStats, "coverage %f\n", cov);
  fprintf(pStats, "genomesize %f\n", totalLen / cov / 1000000 / frac);

  fclose(pStats);
 
    

  Histogram hh;

  string dist_n = outName;
  dist_n += "/lap_dist.ps";
  string scatter_n = outName;
  scatter_n += "/lap_scatter.ps";

  cout << "Plot Dist " << endl;
  hh.Plot(dist_n, laps_clean, 100, 0, 200, color(0.4, 0.2, 0.9));

  MakePNG(dist_n);

  hh.Scatter(scatter_n, size, laps_clean, color(0.3, 0.0, 0.4));

  MakePNG(scatter_n);


  return 0;
}
