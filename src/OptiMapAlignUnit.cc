#ifndef FORCE_DEBUG
#define NDEBUG
#endif


#include "OptiMapAlignUnit.h"


void RSiteRead::Flip() {
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

int RSiteReads::AddRead(const RSiteRead& rr) {
  m_oReads.push_back(rr);
  return m_oReads.isize()-1;
}

void RSiteReads::LoadReads(const string& fileName, int seedSize) {
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
      name = parser.Line();
      RSiteRead rr;
      rr.Name() = name;

      // Obtain the pre/post values and  distmer values
      parser.ParseLine();
      rr.PreDist()  = parser.AsFloat(0);
      for (i=1; i< parser.GetItemCount()-1; i++) // index starts at 1 and ends at length-1 since pre/postDist are loaded separately
        mm.push_back(parser.AsFloat(i));
      rr.PostDist() = parser.AsFloat(i); // Load postDist (last item in distmer sequence)
      rr.Dist() = mm;

      if (mm.isize() >= seedSize) {
        m_oReads.push_back(rr);
        rr.Flip();
        rr.Name() += "_RC";
        m_oReads.push_back(rr);
      }  
      mm.clear();
    }
  }
}  

bool Dmer::operator < (const Dmer & m) const {
  for (int i=0; i<m_data.isize(); i++) {
    if (m_data[i] != m.m_data[i])
      return m_data[i] < m.m_data[i];
  }
  return false;
}
  
bool Dmer::operator != (const Dmer & m) const {
  for (int i=0; i<m_data.isize(); i++) {
    if (m_data[i] != m.m_data[i])
      return true;
  }
  return false;
}

bool Dmer::operator == (const Dmer & m) const {
  for (int i=0; i<m_data.isize(); i++) {
    if (m_data[i] != m.m_data[i])
      return false;
  }
  return true;
}

void Dmer::Print() const {
  for (int i=0; i<m_data.isize(); i++)
    cout << " " << m_data[i];
  cout << " seq: " << m_seq << " pos: " << m_pos << endl;
}

void Dmers::AddSingleReadDmers(const RSiteReads& optiReads , int seedSize, int rIdx) {
  Dmer mm;
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

void Dmers::BuildDmers(const RSiteReads& optiReads , int seedSize) {
  cout << "LOG Build mer list..." << endl;
  for (int rIdx=0; rIdx<optiReads.NumReads(); rIdx++) {
    AddSingleReadDmers(optiReads, seedSize, rIdx);
  }
  cout << "LOG Sort mers... " << m_mers.isize() << endl;
  __gnu_parallel::sort(m_mers.begin(), m_mers.end());
}

void OverlapCandids::AddCandidSort(int rIdx1, int rIdx2, int offsetDelta) {
  // Make sure that rIdx1, rIdx2 are in increasing order (so that sorting will bring all relevant pairs together)
  if(rIdx1>rIdx2) {
    int temp = rIdx1;
    rIdx1 = rIdx2;
    rIdx2 = temp; // Swap rIdx1 & rIdx2
    offsetDelta *= -1;
  }
  m_candids.push_back(OverlapCandid(rIdx1, rIdx2, offsetDelta));
}
 
void OverlapCandids::AddCandid(const OverlapCandid& lapCandid) {
  m_candids.push_back(lapCandid);
}

void OptiMapAlignUnit::WriteLapCandids(const OverlapCandids& candids) {
  for (int i = 0; i<candids.NumCandids(); i++) {
    cout << m_reads[candids[i].GetFirstReadIndex()].Name() << " " << m_reads[candids[i].GetSecondReadIndex()].Name()
         << " " << candids[i].GetOffsetDelta() << endl;
  }
}

void OptiMapAlignUnit::MakeRSites(const string& fileName, const string& motif, int seedSize) {
  FlatFileParser parser;
  parser.Open(fileName);
  string l;
  int i, j;
  svec<int> mm;
  string name;
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    if (parser.Line()[0] == '>') {
      RSiteRead rr;
      rr.Name() = name;
      bool wrotePrefix = false;
      int n = -1;
      for (i=0; i<(int)l.length()-(int)motif.length(); i++) {
	for (j=0; j<motif.length(); j++) {
	  if (motif[j] != toupper(l[i+j]))
	    break;
	}
	if (j == motif.length()) {
	  if (n >= 0) {
            // Obtain the pre/post & dmer values
            if (!wrotePrefix) {
              rr.PreDist() = n; // prefix (number of trailing bits before the first motif location)
              wrotePrefix = true; 
            }
            mm.push_back(i-n);
	  } 
	  n = i;
	}
      }
      if (l != "") {
        if(wrotePrefix) { 
          rr.PostDist() = parser.AsFloat(i); // postfix (number of leading bits after last motif location, last item in distmer sequence)
        } 
      }
      if (mm.isize() >= seedSize) {
        rr.Dist() = mm;
        m_reads.AddRead(rr);
        rr.Flip();
        rr.Name() += "_RC";
        m_reads.AddRead(rr);
      }  
      mm.clear();
      l = "";
      name = parser.Line();
      continue;
    }
    l += parser.Line();
  }
}

void OptiMapAlignUnit::FindLapCandids(int seedSize, OverlapCandids& lapCandids) {
  Dmers  dmers;  // To build dmers from optical reads
  dmers.BuildDmers(m_reads, seedSize); 
  int counter = 0;
  int i,j     = 0;
  cout << "Start going through mers..." << endl;
  lapCandids.ReserveInit(dmers.NumMers());
  for (i=0; i<dmers.NumMers(); i++) {
    counter++;
    if (counter % 10000 == 0)
      cout << "LOG Progress: " << 100*(double)i/(double)dmers.NumMers() << "%" << endl;
    for (j=i+1; j<dmers.NumMers(); j++) {
      if (dmers[j] != dmers[i]) {
        break;
      }
    }
    if (j-i < 25) {
      for (int x = i; x<j; x++) {
        for(int y=x+1; y<j; y++) {
          lapCandids.AddCandidSort(dmers[x].Seq(), dmers[y].Seq(), dmers[y].Pos()-dmers[x].Pos());
        }
      }
    }
    i = j-1;
  }
  cout << "LOG Sort overlap candidates... " << lapCandids.NumCandids() << endl;
  lapCandids.SortAll();
}

void OptiMapAlignUnit::FinalOverlaps(const OverlapCandids& lapCandids, int tolerance,  OverlapCandids& finalOverlaps) {
  cout << "LOG Refine overlap candidates... " << lapCandids.NumCandids() << endl;
  finalOverlaps.ReserveInit(lapCandids.NumCandids()/2); // Rough estimate 
  OverlapCandid currCandid;
  int rejectCnt = 0;
  for(int i=0; i<lapCandids.NumCandids(); i++) {
    if(lapCandids[i]==currCandid) { continue; }
    currCandid = lapCandids[i];  
    const RSiteRead& read1    = m_reads[currCandid.GetFirstReadIndex()]; 
    const RSiteRead& read2    = m_reads[currCandid.GetSecondReadIndex()]; 
    int delta     = currCandid.GetOffsetDelta();
    int readPos1  = delta>=0 ? 0 : -delta;
    int readPos2  = delta>=0 ? delta : 0;
    bool overlaps = true;
    while(readPos1 < read1.Size() && readPos2 < read2.Size()) {
      if(abs(read1[readPos1] - read2[readPos2]) <= tolerance) {
        readPos1++;
        readPos2++;
      } else {
        overlaps = false;
        rejectCnt++;
        break;
      }
    }
    if(overlaps) { finalOverlaps.AddCandid(currCandid); }
  }
  cout << "LOG Total number of overlap candidates rejected: " << rejectCnt << endl;
  cout << "LOG Final number of overlaps: " << finalOverlaps.NumCandids() << endl;
  WriteLapCandids(finalOverlaps);
}
