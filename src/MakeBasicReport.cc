#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "src/Database.h"
#include "src/HTMLTemp.h"


int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandArg<string> tmpCmmd("-t","template file");
  commandArg<string> outCmmd("-o","HTML output");
  commandLineParser P(argc,argv);
  P.SetDescription("Makes an HTML file from basic stats.");
  P.registerArg(fileCmmd);
  P.registerArg(tmpCmmd);
  P.registerArg(outCmmd);
 
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string temp = P.GetStringValueFor(tmpCmmd);
  string out = P.GetStringValueFor(outCmmd);

  HTMLRead html;
  
  html.Read(temp, "##InSeqt##");

  Database db;
  ReadDB(db, fileName);

  html.FillWrite(db, out);

  return 0;
}
