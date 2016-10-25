#include <string>
#include "ryggrad/src/base/CommandLineParser.h"
#include "ryggrad/src/base/StringUtil.h"
#include <unistd.h>




int main( int argc, char** argv )
{

  int i;
  string dir;
  //cout << "Running " << m_exe << endl;
  for (i=strlen(argv[0])-1; i>=0; i--) {
    if (argv[0][i] == '/') {
      char tmp[1024];
      strcpy(tmp, argv[0]);
      tmp[i+1] = 0;
      dir = tmp;
      break;
    }
  }


  commandArg<string> aStringCmmd("-i","input sequence file or file list");
  commandArg<string> outCmmd("-o","output directory");
  commandArg<int> distCmmd("-distance","distance between seeds", 1);
  commandArg<int> numCmmd("-w","width of the filter", 2);
  commandArg<double> fracCmmd("-f","fraction of reads to examine (<=1.)", 0.01);
  commandArg<int> nCmmd("-n","number of total processes", 1);
  commandArg<int> mCmmd("-m","number of parallel processes", 1);

  commandLineParser P(argc,argv);

  P.SetDescription("Runs the overlap finding processes.");
  P.registerArg(aStringCmmd);
  P.registerArg(outCmmd);
  P.registerArg(nCmmd);
  P.registerArg(mCmmd);
  P.registerArg(distCmmd);
  P.registerArg(numCmmd);
  P.registerArg(fracCmmd);

  P.parse();

  string fileName = P.GetStringValueFor(aStringCmmd);
  string out = P.GetStringValueFor(outCmmd);
  int n = P.GetIntValueFor(nCmmd);
  int m = P.GetIntValueFor(mCmmd);
  int dist = P.GetIntValueFor(distCmmd);
  int num = P.GetIntValueFor(numCmmd);
  double frac = P.GetDoubleValueFor(fracCmmd);

  string mkdir = "mkdir " + out;
  int rr2 = system(mkdir.c_str());

  int k = 0;
  for (i=0; i<n; i++) {
    string cmmd = dir + "/FindOverlapCands -i ";
    cmmd += fileName;
    cmmd += " -o ";
    cmmd += out;
    cmmd += " -n " + Stringify(i);
    cmmd += " -c " + Stringify(n);
    cmmd += " -w " + Stringify(num);
    char tmp[256];
    sprintf(tmp, "%f", frac);
    cmmd += " -f ";
    cmmd += tmp;
    k++;
    if (k == m) {
      k = 0;
    } else {
      cmmd += " &";
    }
    int rr = system(cmmd.c_str());
  }

  
  for (i=0; i<n; i++) {
    string outName = out + "/overlapcands.out.";
    outName += Stringify(i);
    outName += ".done";
    while (true) {
      FILE * p = fopen(outName.c_str(), "r");
      if (p != NULL) {
	fclose(p);
	cout << outName << endl;
	break;
      }
      sleep(1);
    }
  }
  cout << "Concatenating output." << endl;
  string cat = "cat " + out + "/overlapcands.out.*.cand > ";
  cat += out + "/overlapcands.out";
  int rrr = system(cat.c_str());
  
  string stats = out + "/overlapcands.txt";
  FILE * pStats = fopen(stats.c_str(), "w");
  fprintf(pStats, "fraction %f\n", frac);
  fclose(pStats);

  cout << "All done!!" << endl;
  return 0;
}
