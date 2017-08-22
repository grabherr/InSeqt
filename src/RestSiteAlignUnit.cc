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

void Dmer::CalcDeviations(svec<int>& deviations, float indelVariance, float deviationCoeff) const {
  if(deviations.isize() < m_data.isize()) { deviations.resize(m_data.isize()); }
  for(int i=0; i<m_data.isize(); i++) {
    int deviation = sqrt(m_data[i]*indelVariance)*deviationCoeff;
    deviations[i] = deviation;
  }
}
 

string Dmer::ToString() const {
  stringstream ss;
  for (int i=0; i<m_data.isize(); i++)
    ss << " " << m_data[i];
  ss << " seq: " << m_seq << " pos: " << m_pos;
  return ss.str();
}

void Dmers::BuildDmers(const RSiteReads& rReads , int dmerLength, int motifLength, int countPerDimension) { 
  m_dmerLength = dmerLength;
  m_dimCount   = countPerDimension;
  m_mers.resize(pow(m_dimCount, m_dmerLength)); // TODO Check to be within memory limit
  SetRangeBounds(motifLength);
  cout << "LOG Build mer list..." << endl;
  FILE_LOG(logDEBUG2) << "LOG Build mer list...";
  for (int rIdx=0; rIdx<rReads.NumReads(); rIdx++) {
    AddSingleReadDmers(rReads, rIdx);
  }
}

void Dmers::SetRangeBounds(int motifSize) {
  double p1     = 1.0 / pow(4, motifSize); //for example for a motif size of 4 this will be 1/256
  double p2     = 1.0 - p1; //for example for a motif size of 4 this will be 255/256
  double pC     = 0;  //Cumulative probability
  int rangeLim  = 0;

  for(int dim=0; dim<m_dimCount-1; dim++) {
    while(pC < (double)(dim+1)/m_dimCount) {
      double pi = pow(p2, rangeLim+1) * p1;
      pC += pi;
      rangeLim++;
      m_dmerCellMap[rangeLim] = dim;
    }
    m_dimRangeBounds.push_back(rangeLim); 
    FILE_LOG(logDEBUG1) << "Dimension Range: " << m_dimRangeBounds.isize()-1 << "  " << rangeLim; 
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
  FILE_LOG(logDEBUG3) << "Read Index: " << rIdx << " total dmers so far: " << m_dmerCount << endl;
}

int Dmers::MapNToOneDim(const svec<int>& nDims) {
  // This function does not do bound checking and assumes that nDims size is m_dmerLength and values are between 0 and m_dimCount
  int mapVal = 0;
  int coeff  = pow(m_dimCount, m_dmerLength-1);
  for(int i=0; i<m_dmerLength; i++) {
    int qVal = m_dimCount - 1; // First set it to the highest possible digit and then check if it belongs in another cell 
    if(nDims[i] < m_dimRangeBounds[m_dimCount-2])  { qVal = m_dmerCellMap[nDims[i]]; } // the last digit range bound is in m_dimCount-2
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

void Dmers::FindNeighbourCells(int initVal, const Dmer& dmer, const svec<int>& deviations, svec<int>& result) {
  FindNeighbourCells(initVal, dmer, deviations, m_dmerLength-1, result);
}

void Dmers::FindNeighbourCells(int initVal, const Dmer& dmer, const svec<int>& deviations, int depth, svec<int>& result) {
  if(depth == -1) { 
    result.push_back(initVal);
    return;
  }
  FindNeighbourCells(initVal, dmer, deviations, depth-1, result);
  svec<int> tempResult = result;
  int currDigit = m_dimCount - 1; // First set it to the highest possible digit and then check if it belongs in another cell 
  if(dmer[depth] < m_dimRangeBounds[m_dimCount-2])  { currDigit = m_dmerCellMap[dmer[depth]]; } // the last digit range bound is in m_dimCount-2
  if((currDigit < m_dimCount-1)  //only add one to the current digit if it has room to be increased 
    && (dmer[depth]+deviations[depth] > m_dimRangeBounds[currDigit])) { // only try one cell up if the deviation limits don't fall within the same cell
    for(int elem:tempResult) {
      int newElem = elem + pow(m_dimCount, (m_dmerLength-1-depth));
      result.push_back(newElem);
    }
  }
}

string MatchCandid::ToString() const {
  stringstream ss;
  ss << GetFirstReadIndex() << "\t" <<  GetFirstMatchPos()  << "\t" 
     << GetSecondReadIndex() << "\t" << GetSecondMatchPos(); 
  return ss.str();
}

void MatchCandids::AddCandidSort(int rIdx1, int rIdx2, int rPos1, int rPos2) {
/*
  // Make sure that rIdx1, rIdx2 are in increasing order (so that sorting will bring all relevant pairs together)
  if(rIdx1>rIdx2) {
    int temp = rIdx1;
    rIdx1 = rIdx2;
    rIdx2 = temp; // Swap rIdx1 & rIdx2
    temp  = rPos1;
    rPos1 = rPos2;
    rPos2 = temp;
  }
  m_candids.push_back(MatchCandid(rIdx1, rIdx2, rPos1, rPos2));
*/
  m_candids[rIdx1][rIdx2]++;
  FILE_LOG(logDEBUG3) << "Adding overlap candidate: " << MatchCandid(rIdx1, rIdx2, rPos1, rPos2).ToString();
}
 
string MatchCandids::ToString() const {
  stringstream ss;
/*
  for(MatchCandid mc:m_candids) {
    ss << mc.ToString() << endl;
  }
*/
  return ss.str();
}

/*
void RestSiteAlignCore::WriteLapCandids(const MatchCandids& candids) {
  for (int i = 0; i<candids.NumCandids(); i++) {
    cout << m_rReads[candids[i].GetFirstReadIndex()].Name() << " " << m_rReads[candids[i].GetSecondReadIndex()].Name()
         << " " << candids[i].GetOffsetDelta() << endl;
  }
}
*/

void RestSiteAlignCore::MakeRSites(const string& fileName, int numOfReads) {
  FlatFileParser parser;
  parser.Open(fileName);
  string l;
  string name;
  //Initialize memory for reads
  m_rReads.Resize(numOfReads*2); //Twice the number of reads to allow for recording reverse complements
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
  FILE_LOG(logDEBUG3) << m_rReads.ToString();
}

void RestSiteAlignCore:: CreateRSitesPerString(const string& origString, const string& origName) {
  if (origString == "" && origName == "") {
    return;
  }
  svec<int> mm;
  RSiteRead rr;
  rr.Name() = origName;
  bool wrotePrefix = false;
  int n = -1;
  for (int i=0; i<(int)origString.length()-(int)m_motif.length(); i++) {
    int j = 0;
    for (j=0; j<m_motif.length(); j++) {
      if (m_motif[j] != toupper(origString[i+j]))
        break;
    }
    if (j == m_motif.length()) {
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
  int readIdx = m_rReads.AddRead(rr);
  m_totalSiteCnt += mm.isize();
  FILE_LOG(logDEBUG3) << "Adding Read: " << readIdx << "  " << rr.Name();
  
  rr.Flip();
  rr.Name() += "_RC";
  readIdx = m_rReads.AddRead(rr);
  m_totalSiteCnt += mm.isize();
  FILE_LOG(logDEBUG3) << "Adding Read: " << readIdx << "  " << rr.Name();
  mm.clear();
}

void RestSiteAlignCore::FindLapCandids(int dmerLength, int motifLength, float indelVariance, float deviationCoeff, MatchCandids& lapCandids) {
  Dmers  dmers;       // To build dmers from restriction site reads
  int dimCount = pow(TotalSiteCount()*3, 1.0/dmerLength);   // Number of bins per dimension
  if(dimCount>35) { 
    FILE_LOG(logWARNING) << "WARN: Input data size is too large ";
    dimCount = 35; 
  }
  FILE_LOG(logINFO) << "Estimate number of Dmers and dimension size for dmer storage: " << TotalSiteCount() << "  " << dimCount; 
  dmers.BuildDmers(m_rReads, dmerLength, motifLength, dimCount);
  cout << "Start iterating through " << dmers.NumMers() <<  "  dmers..." << endl;
  int counter    = 0;
  int matchCount = 0;
  int loopLim = pow(dimCount, dmerLength);
//  lapCandids.ReserveInit(dmers.NumMers());
  svec<int> neighbourCells;
  neighbourCells.reserve(pow(2, dmerLength));
  svec<int> deviations;
  deviations.resize(dmerLength);
  for (int iterIndex=0; iterIndex<loopLim; iterIndex++) {
    counter++;
    if (counter % 100000 == 0) {
      cout << "\rLOG Progress: " << 100*(double)iterIndex/(double)loopLim << "%" << flush;
    }
    if(!dmers[iterIndex].empty()) {
      FILE_LOG(logDEBUG2) << "Number of dmers in cell " << iterIndex << " " << dmers[iterIndex].isize(); 
      for (Dmer dm1:dmers[iterIndex]) {
        neighbourCells.clear();
        dm1.CalcDeviations(deviations, indelVariance, deviationCoeff); //TODO this does not need to be redone every time!
        dmers.FindNeighbourCells(iterIndex, dm1, deviations, neighbourCells); 
        for (int nCell:neighbourCells) {
          FILE_LOG(logDEBUG3) << "Checking Neighbour cell: " << nCell; 
          for (Dmer dm2:dmers[nCell]) {
            FILE_LOG(logDEBUG4) << endl << "Checking dmer match: dmer1 - " << dm1.ToString() << " dmer2 - " << dm2.ToString();
            if(dm1.IsMatch(dm2, deviations)) {
              FILE_LOG(logDEBUG4) << "Found Dmer Match";
                matchCount++;
//              lapCandids.AddCandidSort(dm1.Seq(), dm2.Seq(), dm1.Pos(), dm2.Pos());
            }
          }
        }
      }
    }
  }
//  cout << "LOG Sort overlap candidates... " << lapCandids.NumCandids() << endl;
//  lapCandids.SortAll();
cout << "Total number of matches recorded: " << matchCount << endl;
}

/*
void RestSiteAlignCore::FinalOverlaps(const MatchCandids& lapCandids, int tolerance,  MatchCandids& finalOverlaps) {
  cout << "LOG Refine overlap candidates... " << lapCandids.NumCandids() << endl;
  finalOverlaps.ReserveInit(lapCandids.NumCandids()/2); // Rough estimate 
  MatchCandid currCandid;
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

void RestSiteMapper::GenerateMotifs(int motifLength, int numOfMotifs) {
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
      if(ValidateMotif(motif)) {
        FILE_LOG(logDEBUG3) << "Motif: "  << m_motifs.isize() << " " << motif;
        m_motifs.push_back(motif); 
        if(m_motifs.isize() == numOfMotifs) {
          break;
        }
      }
    } 
  }
}
  
bool RestSiteMapper::ValidateMotif(const string& motif) const {
// 1. Simplicity Filter
  bool simple = true;
  for(int i=0; i<motif.length()-1; i++) {
    if(motif[i]!=motif[i+1]) {
      simple = false; 
    }
  }
  if(simple) { return false; }

// 2. RC Filter
// 3. 
  return true;
}

void RestSiteMapper::CartesianPower(const vector<char>& input, unsigned k, vector<vector<char>>& result) const {
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

void RestSiteMapper::FindMatches(const string& fileName, int readCnt, int motifIndex, MatchCandids& finalOverlaps) const {
  FILE_LOG(logDEBUG2) << "Finding sites/maps for motif: " << m_motifs[motifIndex] << endl;
  RestSiteAlignCore rsaCore(m_motifs[motifIndex]);
  rsaCore.MakeRSites(fileName, readCnt);
  MatchCandids lapCandids;
  rsaCore.FindLapCandids(m_modelParams.DmerLength(), m_modelParams.MotifLength(), 0.08, 1, lapCandids);
}
