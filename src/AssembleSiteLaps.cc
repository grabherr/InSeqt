#define FORCE_DEBUG
#include <string>
#include "ryggrad/src/base/CommandLineParser.h"
#include "ryggrad/src/base/FileParser.h"
#include "ryggrad/src/general/DNAVector.h"


class LongestLap
{
public:
  
  LongestLap() {sites = size = preA = preB = postA = postB = startA = stopA = startB = stopB = 0;}
  
  string name;
  int size;
  int preA;
  int postA;
  int preB;
  int postB;
  string ori;
  int sites;
  
  int startA;
  int stopA;
  int startB;
  int stopB;
  
  string aN;
  string bN;

  void Check(const LongestLap & l) {
    if (l.sites > sites)
      *this = l;
    
  }
};


void Merge(DNAVector & a, DNAVector & b, const LongestLap & l)
{
  if (l.ori == "-")
    b.ReverseComplement();

  // Right extension
  DNAVector m;
  if (l.startA > l.preA) {
    int endA = l.startA;
    int startB = l.preB;
    m.SetToSubOf(a, 0, endA);
    DNAVector mm;
    mm.SetToSubOf(b, startB, b.isize()-startB);
    cout << "(->) Before a=" << a.isize() << " b=" << b.isize() << endl;
    m += mm;
    a = m;
    b.resize(0);
    cout << "(->) After a=" << a.isize() << endl;
  } else { // Left extension;
    int startA = l.preA;
    int endB = l.startB;
    m.SetToSubOf(b, 0, endB);
    DNAVector mm;
    mm.SetToSubOf(a, startA, a.isize()-startA);
    cout << "(<-) Before a=" << a.isize() << " b=" << b.isize() << endl;
    m += mm;
    a = m;
    b.resize(0);
    cout << "(<-) After a=" << a.isize() << endl;

  }

  
}

bool Equal(int a, int b)
{
  int c = a-b;
  if (c<0)
    c = -c;
  if (c < 100)
    return true;
  return false;
}

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","overlap file");
  commandArg<string> inCmmd("-f","fasta file");
  commandArg<string> outCmmd("-o","output file");
  commandLineParser P(argc,argv);
  P.SetDescription("Assembles overlaps from restriction maps.");
  P.registerArg(fileCmmd);
  P.registerArg(inCmmd);
  P.registerArg(outCmmd);
 
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string inName = P.GetStringValueFor(inCmmd);
  string outName = P.GetStringValueFor(outCmmd);

  vecDNAVector dna;
  dna.Read(inName);
  
  svec<LongestLap> laps;
  laps.resize(dna.isize());
    
  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);

  cout << "Finding longest laps" << endl;
  
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;

    LongestLap l;
    l.preA = parser.AsInt(4);
    l.postA = parser.AsInt(5);
    l.preB = parser.AsInt(9);
    l.postB = parser.AsInt(10);
    l.ori = parser.AsString(11);
    l.sites = parser.AsInt(12);

    l.startA = parser.AsInt(2);
    l.stopA = parser.AsInt(3);
    l.startB = parser.AsInt(7);
    l.stopB = parser.AsInt(8);

    l.aN = parser.AsString(1);
    l.bN = parser.AsString(6);

    int index = dna.NameIndex(l.aN); 
    laps[index].Check(l);
   
    
  }
  int i, j;
  cout << "Merging/erasing" << endl;

  for (i=0; i<laps.isize(); i++) {
    const LongestLap & l = laps[i];
    if (l.bN == "")
      continue;
    DNAVector & a = dna(l.aN);
    DNAVector & b = dna(l.bN);

    
    if (a.isize() == 0 || b.isize() == 0)
      continue;
    
    int index = dna.NameIndex(l.bN);
    const LongestLap & g = laps[index];
    //cout << "Longest " << l.aN << " " << l.bN << " <- " << g.aN << " " << g.bN << endl;

    if (Equal(l.startA, l.preA) && Equal(l.stopA, a.isize()-l.postA)) {
      cout << "Erasing " << l.aN << " " << l.startA << " " << l.preA << " <-> " << l.stopA << " " << a.isize()-l.postA << endl;
      a.resize(0);
      continue;
    }
    if (Equal(l.startB, l.preB) && Equal(l.stopB, b.isize()-l.postB)) {
      cout << "Erasing " << l.bN << endl;
      b.resize(0);
      continue;
    }

    if (g.bN == l.aN || g.sites == l.sites) {
      cout << "Found maximal mutal pair, merging " << l.aN << " " << l.bN << endl;
      Merge(a, b, l);
      string name = l.aN + l.bN;
      cout << "New contig: " << name << endl;
      dna.SetName(i, name);
    }

  }



  
  dna.Write(outName, true);
  
  return 0;
}
