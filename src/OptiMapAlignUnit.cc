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
    rr2.Name() += "_RC";
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

void OverlapCandids::AddCandid(int rIdx1, int rIdx2, int offsetDelta) {
  // Make sure that rIdx1, rIdx2 are in increasing order (so that sorting will bring all relevant pairs together)
  if(rIdx1>rIdx2) {
    int temp = rIdx1;
    rIdx1 = rIdx2;
    rIdx2 = temp; // Swap rIdx1 & rIdx2
    offsetDelta *= -1;
  }
  m_candids.push_back(OverlapCandid(rIdx1, rIdx2, offsetDelta));
}
 
void OptiMapAlignUnit::WriteLapCandids(const OverlapCandids& candids) {
  for (int i = 0; i<candids.NumCandids(); i++) {
    cout << m_reads[candids[i].GetFirstReadIndex()].Name() << " " << m_reads[candids[i].GetSecondReadIndex()].Name()
         << " " << candids[i].GetOffsetDelta() << endl;
  }
}

void OptiMapAlignUnit::FindCandidLaps(int seedSize, OverlapCandids& lapCandids) {
  Optimers  optimers;  // To build optimers from optical reads
  optimers.BuildOptimers(m_reads, seedSize); 
  int counter = 0;
  int i,j     = 0;
  cout << "Start going through mers..." << endl;
  lapCandids.ReserveInit(optimers.NumMers());
  for (i=0; i<optimers.NumMers(); i++) {
    counter++;
    if (counter % 10000 == 0)
      cout << "LOG Progress: " << 100*(double)i/(double)optimers.NumMers() << "%" << endl;
    for (j=i+1; j<optimers.NumMers(); j++) {
      if (optimers[j] != optimers[i]) {
        break;
      }
    }
    if (j-i < 25) {
      for (int x = i; x<j; x++) {
        for(int y=x+1; y<j; y++) {
          lapCandids.AddCandid(optimers[x].Seq(), optimers[y].Seq(), optimers[y].Pos()-optimers[x].Pos());
        }
      }
    }
    i = j-1;
  }
  cout << "LOG Sort overlap candidates... " << lapCandids.NumCandids() << endl;
  lapCandids.SortAll();
  WriteLapCandids(lapCandids);
}
