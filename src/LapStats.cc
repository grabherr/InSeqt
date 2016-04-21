#include <string>
#include "base/CommandLineParser.h"
#include "visual/Histogram.h"
#include "base/FileParser.h"

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

  commandArg<string> fileCmmd("-i","input overlap file");
  commandArg<string> readCmmd("-r","input read file");
  commandArg<string> outCmmd("-o","output directory");
  //commandArg<bool> listCmmd("-f","", "");
  commandLineParser P(argc,argv);
  P.SetDescription("Testing the file parser.");
  P.registerArg(fileCmmd);
  P.registerArg(readCmmd);
  P.registerArg(outCmmd);
 
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string readName = P.GetStringValueFor(readCmmd);
  string outName = P.GetStringValueFor(outCmmd);
 
 
  string makeOut = "mkdir ";
  makeOut += outName;
  int rr = system(makeOut.c_str());

  svec<string> reads;
  LoadReads(reads, readName);

  Sort(reads);

  svec<double> laps;
  laps.resize(reads.isize(), 0);

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

    // if (len > 10000)
    //  continue;

    int index = BinSearch(reads, s);
    if (index < 0)
      cout << "ERROR, read not found: " << s << endl;
    laps[index] += 1.;
  }
 
  svec<double> laps_clean;

  for (int j=0; j<laps.isize(); j++) {
    if (laps[j] > 0 && laps[j] < 200) {
      laps_clean.push_back(laps[j]);
    }
  }

  Histogram hh;

  string dist_n = outName;
  dist_n += "/lap_dist.ps";
 
  cout << "Plot Dist " << endl;
  hh.Plot(dist_n, laps_clean, 100, 0, 200, color(0.4, 0.2, 0.9));

  MakePNG(dist_n);
  
  

  return 0;
}
