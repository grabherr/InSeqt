#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include <cmath>
#include <sstream>
#include "RestSiteAlignUnit.h"

string RSiteRead::ToString() const {
  stringstream ss;
  ss << PreDist() << ", ";
  for (int i=0; i<m_dist.isize(); i++) {
    ss << m_dist[i] << ", ";
  }
  ss << PostDist();
  return ss.str();
}

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
  int tmp_pp; //Swap pre/postfix
  tmp_pp     = m_preDist;
  m_preDist  = m_postDist;
  m_postDist = tmp_pp;
}

int RSiteReads::AddRead(const RSiteRead& rr) {
  m_rReads[m_readCount] = rr;
  m_readCount++;
  return m_readCount-1;
}

string RSiteReads::ToString() const {
  string strOut;
  for(int i=0; i<m_readCount; i++) {
    strOut += m_rReads[i].ToString();
    strOut += "\n";
  }
  return strOut;
}

/*
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
        m_rReads.push_back(rr);
        rr.Flip();
        rr.Name() += "_RC";
        m_rReads.push_back(rr);
      }  
      mm.clear();
    }
  }
}  
*/

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

void Dmers::BuildDmers(const RSiteReads& rReads , int dmerLength, int motifLength, int countPerDimension) { 
  m_dmerLength = dmerLength;
  m_dimCount   = countPerDimension;
  m_mers.resize(pow(m_dimCount, m_dmerLength)); // TODO Check to be within memory limit
  cout << "LOG Build mer list..." << endl;
  for (int rIdx=0; rIdx<rReads.NumReads(); rIdx++) {
    AddSingleReadDmers(rReads, rIdx);
  }
}

void Dmers::AddSingleReadDmers(const RSiteReads& rReads, int rIdx) {
  Dmer mm;
  mm.Seq() = rIdx;
  mm.Data().resize(m_dmerLength);
  for (int i=0; i<=rReads[rIdx].Dist().isize()-m_dmerLength; i++) {
    mm.Pos() = i;
    for (int j=0; j<m_dmerLength; j++) {
      mm.Data()[j] = rReads[rIdx].Dist()[i+j];
    }
    //Find where this dmer should be placed in m_mers
    int merLoc = MapNToOneDim(mm.Data());
    m_mers[merLoc].push_back(mm);
    m_dmerCount++;
  }
  cout<<m_dmerCount<<endl;
}

int Dmers::MapNToOneDim(const svec<int>& nDims) {
  // This function does not do bound checking and assumes that nDims size is m_dmerLength and values are between 0 and m_dimCount
  int mapVal = 0;
  int coeff  = pow(m_dimCount, m_dmerLength-1);
  for(int i=0; i<m_dmerLength; i++) {
    int qVal = nDims[i]/30; //TODO fix exp allocation in separate function
    if(qVal>m_dimCount-1) { qVal = m_dimCount-1; }
    mapVal += qVal * coeff;
    coeff  /= m_dimCount;
  }
  return mapVal;
}

svec<int> Dmers::MapOneToNDim(int oneDMappedVal) {
  svec<int> nDims;
  nDims.resize(m_dmerLength);
  int coeff  = m_dimCount;
  for(int i=m_dmerLength-1; i>=0; i--) {
    nDims[i] = oneDMappedVal % coeff;
    oneDMappedVal /= coeff;
  }
  return nDims;
}

void Dmers::FindNeighbourCells(int initVal, svec<int>& result) {
  FindNeighbourCells(initVal, m_dmerLength-1, result);
}

void Dmers::FindNeighbourCells(int initVal, int depth, svec<int>& result) {
  if(depth == -1) { 
    result.push_back(initVal);
    return;
  }
  FindNeighbourCells(initVal, depth-1, result);
  svec<int> tempResult = result;
  int currDigit = (initVal % (int)pow(m_dimCount, depth+1)) / pow(m_dimCount, depth);
  if(currDigit < m_dimCount-1) { //only add one to the current digit if it has room to be increased 
    for(int elem:tempResult) {
      int newElem = elem + pow(m_dimCount, depth);
      result.push_back(newElem);
    }
  }
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
/*
void RestSiteAlignUnit::WriteLapCandids(const OverlapCandids& candids) {
  for (int i = 0; i<candids.NumCandids(); i++) {
    cout << m_rReads[candids[i].GetFirstReadIndex()].Name() << " " << m_rReads[candids[i].GetSecondReadIndex()].Name()
         << " " << candids[i].GetOffsetDelta() << endl;
  }
}
*/
void RestSiteAlignUnit::GenerateMotifs(int motifLength, int numOfMotifs) {
  m_motifs.reserve(numOfMotifs);
  vector<char> alphabet = {'A', 'C', 'G', 'T'}; //Should be in lexographic order
  vector<vector<char>> tempMotifs; 
  CartesianPower(alphabet, motifLength, tempMotifs);
  for(vector<char> sElem:tempMotifs){
    string motif = "";
    for(char cElem:sElem) {
      motif += cElem;
    }
    if(motif.length() == motifLength) { 
      //TODO only add each with RC once
      m_motifs.push_back(motif); 
      if(m_motifs.isize() == numOfMotifs) {
        break;
      }
    } 
  }
  m_rReads.resize(m_motifs.isize());
}
  
void RestSiteAlignUnit::CartesianPower(const vector<char>& input, unsigned k, vector<vector<char>>& result) {
  if (k == 1) {
    for (int value: input) {
      result.push_back( {value} );
    }
    return;
  } else {
    CartesianPower(input, k - 1, result);
    vector<vector<char>> smallerPower = result;
    for (int elem: input) {
      for (vector<char> sublist: smallerPower) {
        sublist.push_back(elem);
        result.push_back(sublist);
      }
    }
    return;
  }
} 

void RestSiteAlignUnit::MakeRSites(const string& fileName, int numOfReads) {
  FlatFileParser parser;
  parser.Open(fileName);
  string l;
  string name;
  //Initialize memory for reads
  for(int mi=0; mi<m_motifs.isize(); mi++) {
    m_rReads[mi].Resize(numOfReads*2); //Twice the number of reads to allow for recording reverse complements
  }
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    if (parser.Line()[0] == '>') {
      CreateRSitesPerString(l, name);
      l = "";
      name = parser.Line();
    }
    l += parser.Line();
  }
  if( l != "") {
    CreateRSitesPerString(l, name);
  }
  for(int i=0; i<m_rReads.isize(); i++) {
      FILE_LOG(logDEBUG3) << m_rReads[i].ToString();
  }
}

void RestSiteAlignUnit:: CreateRSitesPerString(const string& origString, const string& origName) {
  if (origString == "" && origName == "") {
    return;
  }
  svec<int> mm;
  for(int mi=0; mi<m_motifs.isize(); mi++) {
    RSiteRead rr;
    rr.Name() = origName;
    bool wrotePrefix = false;
    int n = -1;
    for (int i=0; i<(int)origString.length()-(int)m_motifs[mi].length(); i++) {
      int j = 0;
      for (j=0; j<m_motifs[mi].length(); j++) {
        if (m_motifs[mi][j] != toupper(origString[i+j]))
          break;
      }
      if (j == m_motifs[mi].length()) {
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
      if (origString != "") {
        if(wrotePrefix) { 
          rr.PostDist() = origString.length() - n - 1; // postfix (number of leading bits after last motif location, last item in dmer sequence)
        } 
      }
    } 
    rr.Dist() = mm;
    m_rReads[mi].AddRead(rr);
    rr.Flip();
    rr.Name() += "_RC";
    m_rReads[mi].AddRead(rr);
    mm.clear();
  }
}

void RestSiteAlignUnit::FindLapCandids(int seedSize, OverlapCandids& lapCandids) {
  Dmers  dmers;  // To build dmers from restriction site reads
  int dmerLen   = 6;
  int dimCount = 30;
  dmers.BuildDmers(m_rReads[0], seedSize, dmerLen, dimCount); //TODO parameterize
  cout << "Start iterating through dmers..." << endl;
  int counter = 0;
  int loopLim = pow(dimCount, dmerLen);
  lapCandids.ReserveInit(dmers.NumMers());
  for (int iterIndex=0; iterIndex<loopLim; iterIndex++) {
    counter++;
    if (counter % 10000 == 0)
      cout << "LOG Progress: " << 100*(double)iterIndex/(double)loopLim << "%" << endl;
    svec<int> neighbourCells;
    neighbourCells.reserve(pow(2, dmerLen));
    dmers.FindNeighbourCells(iterIndex, neighbourCells); 
    for (int nCell:neighbourCells) {
    }
    //lapCandids.AddCandidSort(dmers[x].Seq(), dmers[y].Seq(), dmers[y].Pos()-dmers[x].Pos());
  }
  cout << "LOG Sort overlap candidates... " << lapCandids.NumCandids() << endl;
  lapCandids.SortAll();
}

/*
void RestSiteAlignUnit::FinalOverlaps(const OverlapCandids& lapCandids, int tolerance,  OverlapCandids& finalOverlaps) {
  cout << "LOG Refine overlap candidates... " << lapCandids.NumCandids() << endl;
  finalOverlaps.ReserveInit(lapCandids.NumCandids()/2); // Rough estimate 
  OverlapCandid currCandid;
  int rejectCnt = 0;
  for(int i=0; i<lapCandids.NumCandids(); i++) {
    if(lapCandids[i]==currCandid) { continue; }
    currCandid = lapCandids[i];  
    const RSiteRead& read1    = m_rReads[currCandid.GetFirstReadIndex()]; 
    const RSiteRead& read2    = m_rReads[currCandid.GetSecondReadIndex()]; 
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
*/
