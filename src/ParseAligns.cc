#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"



int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input/output directory");
  commandLineParser P(argc,argv);
  P.SetDescription("Computes mismatch stats from alignments.");
  P.registerArg(fileCmmd);
 
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd) + "/alignments.out";
 

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);

  double n = 0;
  int i;

  int indelT = 0;
  int indelQ = 0;
  int mis = 0;
  double match = 0;
  
  svec<int> allQ;
  svec<int> allT;
  svec<int> allM;
  
  while (parser.ParseLine()) {
    if (parser.GetItemCount() < 4)
      continue;
   
    
    if (parser.AsString(0) != "Query:")
      continue;
    
    string q = parser.AsString(2);
    parser.ParseLine();
    parser.ParseLine();
    string t = parser.AsString(2);

    for (i=0; i<q.length(); i++) {
      if (q[i] != '-' || t[i] != '-')
	n += 1.;

      if (q[i] == '-') {
	indelQ++;
	continue;
      } else {
	if (indelQ > 0) {
	  allQ.push_back(indelQ);
	}
	indelQ = 0;	
      }
      if (t[i] == '-') {
	indelT++;
	continue;
      } else {
	if (indelT > 0) {
	  allT.push_back(indelT);
	}
	indelT = 0;
      }

      if (q[i] == t[i]) {
	match += 1.;
	if (mis > 0) {
	  allM.push_back(mis);
	}
	mis = 0;
      } else {
	mis++;
      }
    }
    
  }
    
  string out = P.GetStringValueFor(fileCmmd) + "/alignstats.txt";

  FILE * p = fopen(out.c_str(), "w");
  fprintf(p, "total %f\n", n);
  fprintf(p, "matches %f\n", match);
  fprintf(p, "identity %f\n", 100*match/n);
  fprintf(p, "insertions %d\n", allT.isize());
  fprintf(p, "deletions %d\n", allQ.isize());
  fprintf(p, "mismatches %d\n", allM.isize());

  
  fclose(p);
  
  return 0;
}
