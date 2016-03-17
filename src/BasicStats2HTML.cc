#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"



int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandLineParser P(argc,argv);
  P.SetDescription("Makes an HTML file from basic stats.");
  P.registerArg(fileCmmd);
 
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
 

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);
  int nLib = 1;
  int nReads = -1;
  double avg = -1.;
  int med = -1;
  int n50 = -1;
  int total = -1;

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    if (parser.AsString(0) == "Reads:")
      nReads = parser.AsInt(1);
    if (parser.AsString(0) == "Average:")
      avg = parser.AsFloat(1);
    if (parser.AsString(0) == "Median:")
      med = parser.AsInt(1);
    if (parser.AsString(0) == "Total:")
      total = parser.AsInt(1);
    if (parser.AsString(0) == "N50:")
      n50 = parser.AsInt(1);
  }



  cout << "  <!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">" << endl;
  cout << "<html>" << endl;
  cout << "  <head>" << endl;
  cout << "" << endl;
  cout << "    <meta http-equiv=\"content-type\" content=\"text/html; charset=UTF-8\">" << endl;
  cout << "    <title>InSeqt Report - Basic Stats</title>" << endl;
  cout << "  </head>" << endl;
  cout << "  <body>" << endl;
  cout << "    <big><big><b><i>InSeqt</i></b><b> Basic Statistics. Data set:</b></big></big><br>" << endl;
  cout << "    <br>" << endl;
  cout << "    <img alt=\"Size distribution\" src=\"size.png\" align=\"left\"" << endl;
  cout << "      height=\"340\" width=\"340\">Number of reads (total): " << nReads << "<br>" << endl;
  cout << "    Number of libraries: " << nLib << " <br>" << endl;
  cout << "    Total sequence: " << total << "<br>" << endl;
  cout << "    Average length: " << avg << "<br>" << endl;
  cout << "    Median length: " << med << "<br>" << endl;
  cout << "    N50: " << n50 << " <br>" << endl;
  cout << "    <br>" << endl;
  cout << "    <br>" << endl;
  cout << "    <br>" << endl;
  cout << "    <br>" << endl;
  cout << "    <br>" << endl;
  cout << "    <br>" << endl;
  cout << "    <br>" << endl;
  cout << "    <br>" << endl;
  cout << "    <br>" << endl;
  cout << "    <br>" << endl;
  cout << "    <br>" << endl;
  cout << "    <br>" << endl;
  cout << "    <br>" << endl;
  cout << "    <br>" << endl;
  cout << "    Average GC content:<br>" << endl;
  cout << "    <img alt=\"GC Content overall\" src=\"seq.png\" height=\"340\" width=\"340\"><img" << endl;
  cout << "      alt=\"GC content in 1kb windows\" src=\"win.png\" height=\"340\"" << endl;
  cout << "      width=\"340\"><br>" << endl;
  cout << "    <br>" << endl;
  cout << "    GC content versus read lengths:<br>" << endl;
  cout << "    <img alt=\"GC versus read lengths\" src=\"scatter.png\" height=\"340\"" << endl;
  cout << "      width=\"340\"><br>" << endl;
  cout << "    <br>" << endl;
  cout << "    <br>" << endl;
  cout << "    <br>" << endl;
  cout << "  </body>" << endl;
  cout << "</html>" << endl;
  

  return 0;
}
