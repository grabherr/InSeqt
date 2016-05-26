#include <string>
#include "base/CommandLineParser.h"
#include "src/DNAVector.h"
#include "visual/Histogram.h"
#include "base/FileParser.h"
#include "src/DataRead.h"

/*
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
  }*/

void ReadLaps(svec<string> & a, svec<string> & b, const string & fileName)
{
  FlatFileParser parser;
  
  parser.Open(fileName);
  parser.ParseLine();
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    string one = parser.AsString(0);
    one += " ";
    one += parser.AsString(1);

    string two = parser.AsString(2);
    two += " ";
    two += parser.AsString(3);
    
    a.push_back(one);
    b.push_back(two);
    
  }

}
 
int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-f","input fasta file");
  commandArg<string> lapCmmd("-l","overlap file");
  commandArg<string> outCmmd("-o","output directory");
  //commandArg<bool> listCmmd("-f","", "");
  commandLineParser P(argc,argv);
  P.SetDescription("Estimates error rates.");
  P.registerArg(fileCmmd);
  P.registerArg(lapCmmd);
  P.registerArg(outCmmd);
 
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string lapName = P.GetStringValueFor(lapCmmd);
  string outName = P.GetStringValueFor(outCmmd);
 
 
  string makeOut = "mkdir ";
  makeOut += outName;
  int rr = system(makeOut.c_str());

  vecDNAVector dna;
  ReadDNA(dna, fileName);
  
  svec<string> a, b;
  ReadLaps(a, b, lapName);

  int i, j;

  int n = 20;
  int k = 0;

  //cmmd = "rm ";
  //cmmd += outName;
  //cmmd += "/*";

  //rr = system(cmmd.c_str());

  string fasta = outName;
  fasta += "/first.fasta";
  
  FILE * pOne = fopen(fasta.c_str(), "w");

  fasta = outName;
  fasta += "/second.fasta";   
  FILE * pTwo = fopen(fasta.c_str(), "w");
 
  string last;
  int counter = 0;
  

  while (counter < n) {
    const string & one = a[k];
    const string & two = b[k];
    k++;


    if (one == last)
      continue;
    last = one;

    counter++;

    const DNAVector & d1 = dna(one);
    const DNAVector & d2 = dna(two);

    fprintf(pOne, ">first_%d\n", counter);
    for (j=0; j<d1.isize(); j++)
      fprintf(pOne, "%c", d1[j]);
    fprintf(pOne, "\n");
   

    fprintf(pTwo, ">second_%d\n", counter);
    for (j=0; j<d2.isize(); j++)
      fprintf(pTwo, "%c", d2[j]);
    fprintf(pTwo, "\n");
 

  }
  
  fclose(pOne);
  fclose(pTwo);


  

  return 0;
}
