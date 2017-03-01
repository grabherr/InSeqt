//#define FORCE_DEBUG
#include <string>
#include "ryggrad/src/base/CommandLineParser.h"
#include "ryggrad/src/base/FileParser.h"




string Clean(const string &s)
{
  int l = (int)s.length();
  if (s[l-2] == 'R' && s[l-1] == 'C') {
    char tmp[1024];
    strcpy(tmp, s.c_str());
    tmp[l-3] = 0;
    return tmp;
  }
  return s;
}


class Pool
{
public:
  Pool() {}

  void push_back(const string & s) {
    m_data.push_back(Clean(s));
  }
  int isize() const {return m_data.isize();}
  string & operator[] (int i) {return m_data[i];}
  const string & operator[] (int i) const {return m_data[i];}
  void clear() {m_data.clear();}
  
  void Absorb(Pool & p ) {
    int i;
    for (i=0; i<p.isize(); i++)
      m_data.push_back(p[i]);
    UniqueSort(m_data);
    p.clear();
  }
  
private:
  svec<string> m_data;
};

class PoolData
{
public:
  PoolData() {}

  void push_back(int s) {
    m_data.push_back(s);
  }
  int isize() const {return m_data.isize();}
  int & operator[] (int i) {return m_data[i];}  
  const int & operator[] (int i) const {return m_data[i];}
  void resize(int n) {
    m_data.resize(n);
  }
  
  void clear() {m_data.clear();}
  
  void Absorb(PoolData & p ) {
    int i;
    for (i=0; i<p.isize(); i++)
      m_data.push_back(p[i]);
    UniqueSort(m_data);
    p.clear();
  }

  void USort() {
    UniqueSort(m_data);
  }
private:
  svec<int> m_data;
};

class SeqLap
{
public:
  SeqLap() {}

  string seq1;
  string seq2;

  bool operator < (const SeqLap & s) const {
    return (seq1 < s.seq1);
  }
  
private:
};


int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandLineParser P(argc,argv);
  P.SetDescription("Find overlaps in restriction maps.");
  P.registerArg(fileCmmd);
 
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);

  
  int i, j;

  svec<SeqLap> laps;
  svec<string> all;
 
  cout << "LOG Reading data." << endl;
  int l = 0;
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;

    SeqLap tmp;
    tmp.seq1 = parser.AsString(0);
    tmp.seq2 = parser.AsString(1);
    laps.push_back(tmp);
    tmp.seq2 = parser.AsString(0);
    tmp.seq1 = parser.AsString(1);
    laps.push_back(tmp);
    all.push_back(parser.AsString(0));
    all.push_back(parser.AsString(1));
 
    if (l % 100000 == 0)
      cout << "LOG Read " << l/1000 << " k" << endl;
    l++;
  }

  cout << "LOG Sorting laps: " << laps.isize() << endl;
  //UniqueSort(laps);
  Sort(laps);

  //for (i=0; i<laps.isize(); i++)
  //cout << laps[i].seq1 << " - " << laps[i].seq2 << endl;
  
  cout << "LOG Sorting reads: " << all.isize() << endl;
  UniqueSort(all);
  
  cout << "LOG Done, have: " << all.isize() << endl;
  bool bOK = false;

  int k = 0;
  svec<int> pool;
  pool.resize(all.isize(), -1);

  for (i=0; i<all.isize(); i++) {
    if (pool[i] != -1)
      continue;

    int curr = i;

    svec<int> members;
    int counter = 0;
    members.push_back(curr);
    pool[curr] = k;
    int depth = 0;
    do {
      curr = members[counter];
      //cout << "Start w/ " << all[curr] << endl;
      SeqLap tmp;
      tmp.seq1 = all[curr];      
      int index = BinSearch(laps, tmp);
      //cout << "Index " << index << endl;
      int added = 0;

      int countCands = 0;
      for (j=index; j<laps.isize(); j++) {
	countCands++;
	if (laps[j].seq1 != all[curr]) {
	  break;
	}
      }

      if (countCands < 10) {
	for (j=index; j<laps.isize(); j++) {
	  //cout << "Try...  " << laps[j].seq1 <<  " " << laps[j].seq2 << " " << all[curr] << endl;
	  if (laps[j].seq1 != all[curr]) {
	    //cout << "BREAK" << endl;
	    break;
	  }
	  //cout << "Found lap w/ " << laps[j].seq2 << endl;
	  int mm = BinSearch(all, laps[j].seq2);
	  if (pool[mm] != -1)
	    continue;
	  pool[mm] = k;
	  members.push_back(mm);
	  added++;
	  //cout << "  members: " << members.isize() << " counter " << counter << endl;
	}
      }
      if (added > 0)
	depth++;
      counter++;
      //cout << "Cont, " << counter << endl;
    } while (counter < members.isize());
    cout << "POOL " << k << " members " << members.isize() << " depth " << depth << endl;
    for (j=0; j<members.isize(); j++) {
      cout << all[members[j]] << endl;
    }
    cout << "END" << endl;
    k++;
  }
  
  
  return 0;
}
