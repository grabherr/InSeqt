//#define FORCE_DEBUG
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
  bool operator == (const OptiMer & m) const {
    for (int i=0; i<m_data.isize(); i++) {
      if (m_data[i] != m.m_data[i])
	return false;
    }
    return true;
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
  }
  
  svec<int> & Dist() {return m_dist;}
  string & Name() {return m_name;}
  const svec<int> & Dist() const {return m_dist;}
  const string & Name() const {return m_name;}
  int & Ori() {return m_ori;}
  const int & Ori() const {return m_ori;}

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
  }

  void AddOptimers(svec<OptiMer> & om, int k, int index) {
    int i, j;
    OptiMer mm;
    mm.Seq() = index;
    mm.Data().resize(k);
    for (i=0; i<=m_dist.isize()-k; i++) {
      mm.Pos() = i;
      for (j=0; j<k; j++) {
	mm.Data()[j] = m_dist[i+j];
      }
      om.push_back(mm);
    }

  }

  
private:
  svec<int> m_dist;
  string m_name;
  int m_ori;
};



class Link
{
public:
  Link() {
    forward = -1;
    back = -1;
    last = -1;
    first = -1;
    id = -1;
  }

  int forward;
  int back;
  int last;
  int first;
  int id;
};


void DoLink(svec<Link> & l, int i, int j, int & id)
{
  Link & a = l[i];
  Link & b = l[j];

  /*
  if (a.forward == -1 && b.back == -1) {
    a.forward = j;
    b.back = i;
    a.last = b.last;
    b.first = a.first;
    return;
  }
  if (a.forward == -1) {
    a.forward = b.first;
    l[b.first].back = i;
    l[b.first].first = a.first;
    a.last = b.last;
    return;
  }
  if (b.back == -1) {
    b.back = a.last;
    l[a.last].forward = j;
    l[a.last].last = b.last;
    b.first = a.first;
    return;
    }*/

  if (a.id == b.id && a.id != -1) {
    //cout << "Nothing to do." << endl;
    return;
  }
  
  //cout << "Link " << i << " to " << j << endl;
  //cout << "a " << a.last << " " <<  b.first << endl;
  
  l[a.last].forward = b.first;
  l[b.first].back = a.last;
  l[a.last].last = l[b.first].last;
  l[b.first].first = l[a.last].first;

  if (a.id == -1 && b.id == -1) {
    a.id = id;
    b.id = id;
    a.last = j;
    b.first = i;
    //cout << "Create cluster " << id << " first: " << a.first << " last: " << b.last << endl;
    id++;
  } else {
    int first = a.first;
    int last = b.last;
    int next = a.first;
    //cout << "First: " << a.first << " last: " << b.last << endl; 
    do {
      //cout << reads[next].Name() << endl;
      l[next].first = first;
      l[next].last = last;
      //cout << "Update " << ext << " to " << first << " " << last << endl;
      next = l[next].forward;
      
    } while (next != -1);
    

    if (l[i].id != -1) {
      next = j;
      do {
	//cout << reads[next].Name() << endl;
	l[next].id = l[i].id;
	next = l[next].forward;
      } while (next != -1);
    } else {
      next = i;
      do {
	//cout << reads[next].Name() << endl;
	l[next].id = l[j].id;
	next = l[next].back;
      } while (next != -1);

    }
  }
}


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
      rr.Dist() = mm;
      rr.Name() = name;
      if (name != "" && mm.isize() >= k) {
	reads.push_back(rr);
	rr.Flip();
	rr.Name() += " RC";
	//reads.push_back(rr);
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
    rr2.Dist() = mm;
    rr2.Name() = name;
    reads.push_back(rr2);
    rr2.Name() += " RC";
    rr2.Flip();
    //reads.push_back(rr2);
  }
  
  cout << "Build mer list..." << endl;
  svec<OptiMer> mers;
  for (i=0; i<reads.isize(); i++) {
    reads[i].AddOptimers(mers, k, i);
  }

  cout << "Sort mers..." << endl;
  Sort(mers);

  svec<Link> link;
  link.resize(reads.isize());
  for (i=0; i<link.isize(); i++) {
    link[i].last = i;
    link[i].first = i;
  }
  
  int kp = 0;

  int counter = 0;
  int cluster = 0;
  for (i=0; i<mers.isize(); i++) {
    counter++;
    if (counter % 10000 == 0)
      cout << "Progress: " << 100*(double)i/(double)mers.isize() << "%" << endl;
    for (j=i+1; j<mers.isize(); j++) {
      if (mers[j] != mers[i]) {
	break;
      }
      //cout << "same " << i << " " << j << endl;
    }
    for (int x = i+1; x<j; x++) {
      //cout << "Loop " << reads[mers[x].Seq()].Name() << " " << reads[mers[x-1].Seq()].Name() << " " << i << " " << j << " " << x << endl;
      DoLink(link, mers[x-1].Seq(), mers[x].Seq(), cluster);

	/*
      //cout << "Pooling" << endl;
      Link & t = pool[mers[x].Seq()];
      if (t != -1) {
	//cout << "Adjust from " << t << " to " << kp << endl;
	for (int y=0; y<pool.isize(); y++) {
	  if (pool[y] == t) {
	    pool[y] = kp;
	  }
	}
      }
      //cout << "Assign " 
      pool[mers[x].Seq()] = kp;
	*/
    }
    kp++;
    i = j-1;
    
  }

  svec<int> done;
  done.resize(reads.isize(), 0);

  /*
  cout << "DEBUG" << endl;
  for (int i=0; i<reads.isize(); i++) {
    cout << "Read " << i << " " << reads[i].Name() << " ";
    cout << " <- " << link[i].back << " " << link[i].forward << " -> ";
    cout << " f " << link[i].first << " r " << link[i].last << endl;
   }
  */
  
  int pp = 0;
  for (i=0; i<reads.isize(); i++) {
    if (done[i])
      continue;
    Link & l = link[i];
    if (l.first != i)
      continue;
    int next = i;
    cout << "Pool " << pp << endl;
    pp++;
    do {
      cout << reads[next].Name() << endl;
      //<< " " << l.first << " " << link[next].id << endl;
      //cout << " -> " << link[next].forward << " <- " << link[next].back << endl;
      next = link[next].forward;
    } while (next != -1);
  }

  
  return 0;
}
