#include <string>
#include "ryggrad/src/base/CommandLineParser.h"
#include "ryggrad/src/base/FileParser.h"
#include "Database.h"
#include "HTMLTemp.h"

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

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input directory");
  commandArg<string> tmpCmmd("-t","template file", "");
  commandArg<string> outCmmd("-o","HTML output file", "");
  commandLineParser P(argc,argv);
  P.SetDescription("Makes an HTML file from overlap stats.");
  P.registerArg(fileCmmd);
  P.registerArg(tmpCmmd);
  P.registerArg(outCmmd);
 
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string temp = P.GetStringValueFor(tmpCmmd);
  string out = P.GetStringValueFor(outCmmd);

  if (out == "")
    out = fileName + "/InSeqtOverlapReport.html";

  fileName += "/overlapstats.txt";

  string exe;
   
  int i;

  for (i=strlen(argv[0])-1; i>=0; i--) {
    if (argv[0][i] == '/') {
      char tmp[1024];
      strcpy(tmp, argv[0]);
      tmp[i+1] = 0;
      exe = tmp;
      break;
    }
  }
  if (temp == "") {
    temp = exe;
    temp += "../templates/lap_template.html";
  }

  HTMLRead html;
  html.SetRelativePath(P.GetStringValueFor(fileCmmd));
  html.Read(temp, "##InSeqt##");

  Database db;
  ReadDBPlain(db, fileName);

  html.FillWrite(db, out);

  return 0;
}
