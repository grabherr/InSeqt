#ifndef FORCE_DEBUG
#define NDEBUG
#endif


#include "OptiMapAlignUnit.h"


void OptiRead::Flip() {
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

void OptiReads::LoadReads(const string& fileName, int seedSize) {
  FlatFileParser parser;
  parser.Open(fileName);

  string l;
  int i, j;

  svec<int> mm;
  string name;

  cout << "LOG Read data..." << endl;

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    if (parser.Line()[0] == '>') {
      OptiRead rr;
      rr.Dist() = mm;
      rr.Name() = name;
      if (name != "" && mm.isize() >= seedSize) {
        m_oReads.push_back(rr);
        rr.Flip();
        rr.Name() += "_RC";
        m_oReads.push_back(rr);
      }  
      
      name = parser.Line();
      mm.clear();
      continue;
    }
    for (i=0; i< parser.GetItemCount(); i++)
      mm.push_back(parser.AsFloat(i));
  }
  if (mm.isize() > seedSize) {
    OptiRead rr2;
    rr2.Dist() = mm;
    rr2.Name() = name;
    m_oReads.push_back(rr2);
    rr2.Name() += " RC";
    rr2.Flip();
    m_oReads.push_back(rr2);
  }
}  

bool Optimer::operator < (const Optimer & m) const {
  for (int i=0; i<m_data.isize(); i++) {
    if (m_data[i] != m.m_data[i])
      return m_data[i] < m.m_data[i];
  }
  return false;
}
  
bool Optimer::operator != (const Optimer & m) const {
  for (int i=0; i<m_data.isize(); i++) {
    if (m_data[i] != m.m_data[i])
      return true;
  }
  return false;
}
bool Optimer::operator == (const Optimer & m) const {
  for (int i=0; i<m_data.isize(); i++) {
    if (m_data[i] != m.m_data[i])
      return false;
  }
  return true;
}

void Optimer::Print() const {
  for (int i=0; i<m_data.isize(); i++)
    cout << " " << m_data[i];
  cout << " seq: " << m_seq << " pos: " << m_pos << endl;
}

void Optimers::AddSingleReadOptimers(const OptiReads& optiReads , int seedSize, int rIdx) {
  Optimer mm;
  mm.Seq() = rIdx;
  mm.Data().resize(seedSize);
  for (int i=0; i<=optiReads[rIdx].Dist().isize()-seedSize; i++) {
    mm.Pos() = i;
    for (int j=0; j<seedSize; j++) {
      mm.Data()[j] = optiReads[rIdx].Dist()[i+j];
    }
    m_mers.push_back(mm);
  }
}

void Optimers::BuildOptimers(const OptiReads& optiReads , int seedSize) {
  cout << "LOG Build mer list..." << endl;
  for (int rIdx=0; rIdx<optiReads.NumReads(); rIdx++) {
    AddSingleReadOptimers(optiReads, seedSize, rIdx);
  }
  cout << "LOG Sort mers... " << m_mers.isize() << endl;
  Sort(m_mers);
}


void OptiMapAlignUnit::PoolReads() {
  int kp      = 0;
  int counter = 0;
  int cluster = 0;
  int i,j     = 0;
  ORLinks link(m_reads.NumReads());
  cout << "Start going through mers..." << endl;
  for (i=0; i<m_optimers.NumMers(); i++) {
    counter++;
    if (counter % 1000 == 0)
      cout << "LOG Progress: " << 100*(double)i/(double)m_optimers.NumMers() << "%" << endl;
    for (j=i+1; j<m_optimers.NumMers(); j++) {
      if (m_optimers[j] != m_optimers[i]) {
        break;
      }
      //cout << "same " << i << " " << j << endl;
    }
    if (j-i < 25) {
      for (int x = i; x<j; x++) {
        cout << m_reads[m_optimers[x].Seq()].Name() << " ";
      //cout << "Loop " << m_reads[m_optimers[x].Seq()].Name() << " " << m_reads[m_optimers[x-1].Seq()].Name() << " " << i << " " << j << " " << x << endl;
      //DoLink(link, m_optimers[x-1].Seq(), m_optimers[x].Seq(), cluster);
	
      /*
      //cout << "Pooling" << endl; 
      int t = pool[m_optimers[x].Seq()];
      if (t != -1) {
      //cout << "Adjust from " << t << " to " << kp << endl;
      for (int y=0; y<pool.isize(); y++) {
        if (pool[y] == t) {
          pool[y] = kp;
	}
      }
    }
    //cout << "Assign " 
    pool[m_optimers[x].Seq()] = kp;
    */
      }
      kp++;
    }
    if (j > i)
      cout << endl;
    i = j-1;
  }

  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  return;
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  svec<int> done;
  done.resize(m_reads.NumReads(), 0);

  /*
  cout << "DEBUG" << endl;
  for (int i=0; i<m_reads.isize(); i++) {
    cout << "Read " << i << " " << m_reads[i].Name() << " ";
    cout << " <- " << link[i].back << " " << link[i].GetForward() << " -> ";
    cout << " f " << link[i].GetFirst() << " r " << link[i].last << endl;
   }
  */
  
  int pp = 0;
  for (int i=0; i<m_reads.NumReads(); i++) {
    if (done[i])
      continue;
    const ORLink & l = link[i];
    if (l.GetFirst() != i)
      continue;
    int next = i;
    cout << "Pool " << pp << endl;
    pp++;
    do {
      cout << m_reads[next].Name() << endl;
      //<< " " << l.GetFirst() << " " << link[next].id << endl;
      //cout << " -> " << link[next].GetForward() << " <- " << link[next].back << endl;
      next = link[next].GetForward();
    } while (next != -1);
  }  
}
