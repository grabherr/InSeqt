#define FORCE_DEBUG
#include <string>
#include "ryggrad/src/base/CommandLineParser.h"
#include "ryggrad/src/base/FileParser.h"


class OptiMer
{
public:
  OptiMer() {
    m_seq = -1;
    m_pos = -1;
  }

  bool operator < (const OptiMer & m) const {
    for (int i=0; i<m_data.isize(); i++) {
      if (m_data[i] != m.m_data[i])
	return m_data[i] < m.m_data[i];
    }
    return false;
  }
  
  bool operator != (const OptiMer & m) const {
    for (int i=0; i<m_data.isize(); i++) {
      if (m_data[i] != m.m_data[i])
	return true;
    }
    return false;
  }

  int & Seq() {return m_seq;}
  int & Pos() {return m_pos;}
 
  const int & Seq() const {return m_seq;}
  const int & Pos() const {return m_pos;}
 
  svec<int> & Data() {return m_data;}
  const svec<int> & Data() const {return m_data;}

  void Print() const {
    for (int i=0; i<m_data.isize(); i++)
      cout << " " << m_data[i];
    cout << " seq: " << m_seq << " pos: " << m_pos << endl;
  }
  
private:
  svec<int> m_data;
  int m_seq;
  int m_pos;
};


class OptiRead
{
public:
  OptiRead() {
    m_ori = 1;
    m_pre = m_post = 0;
  }
  
  svec<int> & Dist() {return m_dist;}
  string & Name() {return m_name;}
  const svec<int> & Dist() const {return m_dist;}
  const string & Name() const {return m_name;}
  int & Ori() {return m_ori;}
  const int & Ori() const {return m_ori;}
  
  int & Pre() {return m_pre;}
  const int & Pre() const {return m_pre;}
  int & Post() {return m_post;}
  const int & Post() const {return m_post;}

  void Flip() {
    m_ori = -m_ori;
    svec<int> tmp;
    tmp.resize(m_dist.isize());
    int k = m_dist.isize()-1;
    for (int i=0; i<m_dist.isize(); i++) {
      tmp[k] = m_dist[i];
      k--;
    }
    m_dist = tmp;
    int swap = m_pre;
    m_pre = m_post;
    m_post = swap;
  }

  bool FromList(const svec<int> & l) {
    int i;
    if (l.isize() > 2) {   
      m_dist.resize(l.isize()-2);
      for (i=1; i<l.isize()-1; i++)
	m_dist[i-1] = l[i];
      m_pre = l[0];
      m_post = l[l.isize()-1];
      return true;
    } else {
      m_dist.clear();
      m_pre = m_post = 0;
      return false;
    }
  }

  
  void AddOptimers(svec<OptiMer> & om, int k, int index) {
    int i, j;
    OptiMer mm;
    mm.Seq() = index;
    mm.Data().resize(k);
    svec<int> red;
    red.resize(k);
    //cout << "check... " << red.isize() << " k=" << k << endl;
    for (i=0; i<=m_dist.isize()-k; i++) {
      mm.Pos() = i;
      red.resize(k);
      for (j=0; j<k; j++) {
	//cout << i << " " << j << " " << m_dist.isize() << " ... " << red.isize() << " k=" << k << endl;
	//cout << " -> " << m_dist[i+j] << endl;
	mm.Data()[j] = m_dist[i+j];
	red[j] = m_dist[i+j];
      }
      //cout << "Done" << endl;
      UniqueSort(red);
      if (red.isize() >= k-1)
	om.push_back(mm);
    }

  }

  
private:
  svec<int> m_dist;
  
  string m_name;
  int m_ori;
  int m_pre;
  int m_post;
  
};






int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandArg<int> kCmmd("-k","seed size", 6);
  commandArg<int> wCmmd("-w","wiggle", 5);
  commandLineParser P(argc,argv);
  P.SetDescription("Find overlaps in restriction maps.");
  P.registerArg(fileCmmd);
  P.registerArg(kCmmd);
  P.registerArg(wCmmd);
 
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  int k = P.GetIntValueFor(kCmmd);
  int w = P.GetIntValueFor(wCmmd);

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);

  string l;
 
  int i, j;

  svec<int> mm;
  string name;


  svec<OptiRead> reads;
  cout << "Read data..." << endl;

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    if (parser.Line()[0] == '>') {
      OptiRead rr;
      rr.FromList(mm);
      rr.Name() = name;
      if (name != "" && mm.isize() >= k) {
	reads.push_back(rr);
	rr.Flip();
	rr.Name() += " RC";
	reads.push_back(rr);
      }
      
      name = parser.Line();
      mm.clear();
      continue;
    }
    for (i=0; i< parser.GetItemCount(); i++)
      mm.push_back(parser.AsFloat(i));

    
  }
  if (mm.isize() > k) {
    OptiRead rr2;
    rr2.FromList(mm);
    rr2.Name() = name;
    reads.push_back(rr2);
    rr2.Name() += " RC";
    rr2.Flip();
    reads.push_back(rr2);
  }
  
  cout << "Build mer list..." << endl;
  svec<OptiMer> mers;
  for (i=0; i<reads.isize(); i++) {
    reads[i].AddOptimers(mers, k, i);
  }

  cout << "Sort mers..." << endl;
  Sort(mers);


  //for (i=0; i<mers.isize(); i++)
  //mers[i].Print();

  
  cout << "Lookup mers..." << endl;
  svec<int> busy;
  busy.resize(reads.isize(), 0);

  int maxMis = 2;
  
  for (i=0; i<reads.isize(); i+=2) {
    svec<OptiMer> tmp;
    reads[i].AddOptimers(tmp, k, 0);
    cout << "Test read " << i << " " << reads[i].Name() << endl;
    
    svec<int> cands;
    for (j=0; j<tmp.isize(); j+= k) {
      int index = BinSearch(mers, tmp[j]);
       if (index < 0)
	continue;
      for (int x=index; x<mers.isize(); x++) {
	/*cout << "Found " << mers[x].Seq() << " " << mers[x].Pos() << endl;
	for (int z=0; z<mers[x].Data().isize(); z++)
	  cout << " " << tmp[j].Data()[z];
	cout << endl;
	for (int z=0; z<mers[x].Data().isize(); z++)
	  cout << " " << mers[x].Data()[z];
	  cout << endl;*/
	if (mers[x].Seq() == i)
	  continue;
 	if (mers[x] != mers[index]) {
	  //cout << "BREAK!" << endl;
	  break;
	}

	bool bNew = true;
	int y;

	if (busy[mers[x].Seq()] == 1)
	  continue;
	/*
	//cout << "size " << cands.isize() << endl;
	for (y=0; y<cands.isize(); y++) {
	  //cout << "??" << cands[y] << " " << mers[x].Seq() << endl;
	  if (cands[y] == mers[x].Seq()) {
	    bNew = false;
	    break;
	  }
	  }*/

	//cout << "NEW " << bNew << endl;
	//	if (bNew) {
	int shift = mers[x].Pos() - j;
	cout << "Overlap " << reads[i].Name() << " vs " << reads[mers[x].Seq()].Name() <<  " " << shift;
	cout << " " << mers[x].Pos() << " " << j << endl;

	cands.push_back(mers[x].Seq());
	busy[mers[x].Seq()] = 1;

	int lap = 0;
	int same = 0;
	
	int nnn = reads[i].Dist().isize()+reads[mers[x].Seq()].Dist().isize();
	int offA = reads[i].Pre();
	int offB = reads[mers[x].Seq()].Pre();

	int startA = -1;
	int startB = -1;
	int stopA = -1;
	int stopB = -1;
	
	int firstA = -1;
	int firstB = -1;
	
	for (y=-nnn; y<nnn; y++) {
	  int a = y;
	  int b = y+shift;
	  //int b = y;
	    //int a = y+shift;
	  bool bDo = false;
	  string toPrint;
	  int da = -1;
	  int db = -2;
	  if (a >=0 && a < reads[i].Dist().isize()) {
	    toPrint += Stringify(reads[i].Dist()[a]);
	    da = reads[i].Dist()[a];
	    offA += da;
	    bDo = true;
	    if (firstA < 0)
	      firstA = da;
	  } else {
	    toPrint += "---";
	  }
	  toPrint += "  ";
	  if (b >=0 && b < reads[mers[x].Seq()].Dist().isize()) {
	    toPrint += Stringify(reads[mers[x].Seq()].Dist()[b]);
	    db = reads[mers[x].Seq()].Dist()[b];
	    offB += db;
	    if (firstB < 0)
	      firstB = da;
	    bDo = true;
	  } else {
	    toPrint += "---";
	  }
	  if (bDo) {
	    if (da == db) {
	      same++;
	      if (startA == -1) {
		startA = offA;
		startB = offB;
	      }
	      stopA = offA;
	      stopB = offB;
	    }
	    if (da >=0 && db >= 0)
	      lap++;
		
	    cout << toPrint << endl;
	  }
	
	}
	if (lap-same < maxMis) {
	  cout << "PASSED" << endl;
	  cout << "OVERLAP: " << reads[i].Name() << " " << startA << " " << stopA << " ";
	  cout << reads[i].Pre() + firstA << " " << reads[i].Post() << " ";
	  StringParser ppp;
	  ppp.SetLine(reads[mers[x].Seq()].Name());
	  cout << ppp.AsString(0);
	  cout << " " << startB << " " << stopB << " ";
	  cout << reads[mers[x].Seq()].Pre() + firstB << " " << reads[mers[x].Seq()].Post() << " ";
	  if (strstr(reads[mers[x].Seq()].Name().c_str(), " RC") == NULL)
	    cout << "+";
	  else
	    cout << "-";
	  cout << " " << same << endl;
	  
	} else {
	  cout << "FAILED" << endl;
	}
      }
    }
    for (j=0; j<cands.isize(); j++)
      busy[cands[j]] = 0;
  }


  
  return 0;
}
