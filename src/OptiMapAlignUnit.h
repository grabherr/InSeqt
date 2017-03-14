#ifndef OPTIMAPALIGNUNIT_H
#define OPTIMAPALIGNUNIT_H

#include <string>
#include "ryggrad/src/base/CommandLineParser.h"
#include "ryggrad/src/base/FileParser.h"

class OptiRead
{
public:
  OptiRead() {
    m_ori = 1;
  }
  
  int operator[](int idx) const  { return m_dist[idx];    }
  svec<int> & Dist()             { return m_dist;         }
  string & Name()                { return m_name;         }
  int & Ori()                    { return m_ori;          }
  const svec<int> & Dist() const { return m_dist;         }
  const string & Name() const    { return m_name;         }
  const int & Ori() const        { return m_ori;          }
  int Size() const               { return m_dist.isize(); }

  void Flip();

private:
  svec<int> m_dist;
  string m_name;
  int m_ori;
};

class OptiReads 
{
public:
  OptiReads(): m_oReads() {}
  const OptiRead& operator[](int idx) const { return m_oReads[idx];     }
  int NumReads() const                      { return m_oReads.isize();  }  

  void LoadReads(const string& fileName, int seedSize);

private:
  svec<OptiRead> m_oReads;  /// List of optical reads
};

class Optimer
{
public:
  Optimer() {
    m_seq = -1;
    m_pos = -1;
  }

  bool operator < (const Optimer & m) const;
  bool operator != (const Optimer & m) const; 
  bool operator == (const Optimer & m) const; 

  int & Seq() {return m_seq;}
  int & Pos() {return m_pos;}
 
  const int & Seq() const {return m_seq;}
  const int & Pos() const {return m_pos;}
 
  svec<int> & Data() {return m_data;}
  const svec<int> & Data() const {return m_data;}

  void Print() const;

private:
  svec<int> m_data;
  int m_seq;
  int m_pos;
};

class Optimers 
{
public:
  Optimers(): m_mers() {}

  const Optimer& operator[](int idx) const { return m_mers[idx];     }
  int NumMers()                            { return m_mers.isize();  }  

  void AddSingleReadOptimers(const OptiReads& optiReads , int seedSize, int rIdx);
  void BuildOptimers(const OptiReads& optiReads , int seedSize); 

private:
  svec<Optimer> m_mers;  /// List of optimes
};

class OverlapCandid
{
public:
  OverlapCandid(): m_rIdx1(-1), m_rIdx2(-1), m_offsetDelta(-1) {}
  OverlapCandid(int idx1, int idx2, int delta): m_rIdx1(idx1), m_rIdx2(idx2), m_offsetDelta(delta) {}

  int GetFirstReadIndex() const  { return m_rIdx1;       }
  int GetSecondReadIndex() const { return m_rIdx2;       }
  int GetOffsetDelta() const     { return m_offsetDelta; }

  inline bool operator < (const OverlapCandid& rhs) const {
    return(tie(m_rIdx1, m_rIdx2, m_offsetDelta)
       < tie(rhs.m_rIdx1, rhs.m_rIdx2, rhs.m_offsetDelta)); // keep the same order
  }

  inline bool operator == (const OverlapCandid& rhs) const {
    return(tie(m_rIdx1, m_rIdx2, m_offsetDelta)
       == tie(rhs.m_rIdx1, rhs.m_rIdx2, rhs.m_offsetDelta)); // keep the same order
  }


private:
  int m_rIdx1;        /// The index of the first read
  int m_rIdx2;        /// The index of the second read 
  int m_offsetDelta;  /// offset2 - offset1  Where the offsets are where the seed match starts from
};

class OverlapCandids
{
public:
  OverlapCandids(): m_candids() {}
  
  const OverlapCandid& operator[](int idx) const { return m_candids[idx]; }

  int NumCandids() const { return m_candids.isize(); }

  void ReserveInit(int initialCapacity) { m_candids.reserve(initialCapacity); } 

  void AddCandidSort(int rIdx1, int rIdx2, int offsetDelta);
  void AddCandid(const OverlapCandid& lapCandid);
  
  void SortAll() { Sort(m_candids); }

private:
  svec<OverlapCandid> m_candids; /// unordered set of overlap candidates (for uniqueness)
};

class OptiMapAlignUnit 
{
public:
  OptiMapAlignUnit(): m_reads() {}
  
  void LoadReads(const string& fileName, int seedSize)  { m_reads.LoadReads(fileName, seedSize); }

  void FindLapCandids(int seedSize, OverlapCandids& lapCandids);

  void FinalOverlaps(const OverlapCandids& lapCandids, int tolerance, OverlapCandids& finalOverlaps);

  void WriteLapCandids(const OverlapCandids& candids);

private:
  OptiReads m_reads;     /// Optical Reads
};
#endif //OPTIMAPALIGNUNIT_H
