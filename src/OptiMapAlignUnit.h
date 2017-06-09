#ifndef OPTIMAPALIGNUNIT_H
#define OPTIMAPALIGNUNIT_H

#include <string>
#include <parallel/algorithm>
#include "ryggrad/src/base/CommandLineParser.h"
#include "ryggrad/src/base/FileParser.h"

class RSiteRead
{
public:
  RSiteRead() {
    m_ori = 1;
  }
  
  int operator[](int idx) const  { return m_dist[idx];    }
  svec<int> & Dist()             { return m_dist;         }
  string & Name()                { return m_name;         }
  int& PreDist()                 { return m_preDist;      }
  int& PostDist()                { return m_postDist;     }
  int & Ori()                    { return m_ori;          }
  const svec<int> & Dist() const { return m_dist;         }
  const int& PreDist() const     { return m_preDist;      }
  const int& PostDist() const    { return m_postDist;     }
  const string & Name() const    { return m_name;         }
  const int & Ori() const        { return m_ori;          }
  int Size() const               { return m_dist.isize(); }

  string ToString() const;
  void Flip();

private:
  svec<int> m_dist;       /// Distmer values
  int m_preDist;          /// Number of bits prior to the start of the first distmer value
  int m_postDist;         /// Number of bits left over after the last distmer value
  string m_name;          /// Name of optiRead
  int m_ori;              /// Orientation
};

class RSiteReads 
{
public:
  // Default Ctor
  RSiteReads(): m_readCount(0), m_rReads() {}
  // Ctor 1
  RSiteReads(int numOfReads): m_readCount(0), m_rReads() {
    Resize(numOfReads);   
  }

  void Resize(int size)                      { m_rReads.resize(size);    }

  const RSiteRead& operator[](int idx) const { return m_rReads[idx];     }
  RSiteRead& operator[](int idx)             { return m_rReads[idx];     }
  int NumReads() const                       { return m_readCount;       }

  int AddRead(const RSiteRead& rr); 
  string ToString() const;

  //void LoadReads(const string& fileName, int seedSize);

private:
  int m_readCount;
  svec<RSiteRead> m_rReads;  /// List of site reads
};

class Dmer
{
public:
  Dmer() {
    m_seq = -1;
    m_pos = -1;
  }

  bool operator < (const Dmer & m) const;
  bool operator != (const Dmer & m) const; 
  bool operator == (const Dmer & m) const; 

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

class Dmers 
{
public:
  Dmers(): m_mers() {}

  const Dmer& operator[](int idx) const { return m_mers[idx];     }
  int NumMers()                            { return m_mers.isize();  }  

  void AddSingleReadDmers(const RSiteReads& rReads , int seedSize, int rIdx);
  void BuildDmers(const RSiteReads& rReads , int seedSize); 

private:
  svec<Dmer> m_mers;  /// List of optimes
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
  
  void SortAll() { 
    __gnu_parallel::sort(m_candids.begin(), m_candids.end());
  }

private:
  svec<OverlapCandid> m_candids; /// unordered set of overlap candidates (for uniqueness)
};

class OptiMapAlignUnit 
{
public:
  OptiMapAlignUnit(): m_rReads(), m_motifs() {}

  /* Generate Permutation of the given alphabet to reach number of motifs required */
  void GenerateMotifs(int motifLength, int numOfMotifs);  

  /* Function to obtain Permutations of string characters */
  void Permutation(string alphabet, int startLen, int len);

  /* Function to swap two characters */
  void Swap(char& a, char& b);

  void MakeRSites(const string& fileName, int numOfReads); 
  void CreateRSitesPerString(const string& origString, const string& origName); 

  //void LoadReads(const string& fileName, int seedSize)  { m_rReads.LoadReads(fileName, seedSize); }

  void FindLapCandids(int seedSize, OverlapCandids& lapCandids);

  void FinalOverlaps(const OverlapCandids& lapCandids, int tolerance, OverlapCandids& finalOverlaps);

  void WriteLapCandids(const OverlapCandids& candids);

private:
  svec<RSiteReads> m_rReads;        /// Restriction Site reads per motif
  svec<string> m_motifs;            /// Vector of all motifs for which restriction site reads have been generated
};
#endif //OPTIMAPALIGNUNIT_H
