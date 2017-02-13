#define FORCE_DEBUG
#include <string>
#include "ryggrad/src/base/CommandLineParser.h"
#include "ryggrad/src/base/FileParser.h"



int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandArg<string> seqCmmd("-s","sequence", "CATG");
  commandLineParser P(argc,argv);
  P.SetDescription("Build an in silico restriction map.");
  P.registerArg(fileCmmd);
  P.registerArg(seqCmmd);
 
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string site = P.GetStringValueFor(seqCmmd);

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);

  string l;
 
  int i, j;
  
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    if (parser.Line()[0] == '>') {
      int n = -1;
      for (i=0; i<(int)l.length()-(int)site.length(); i++) {
	//cout << i << " " << l.length() << endl;
	for (j=0; j<site.length(); j++) {
	  if (site[j] != l[i+j])
	    break;
	}
	if (j == site.length()) {
	  if (n >= 0) {
	    cout << " " << i-n;
	  } 
	  n = i;
	}
      }
      if (l != "")
	cout << endl;
      l = "";
      cout << parser.Line() << endl;
      continue;
    }
    l += parser.Line();
  }
  return 0;
}
