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

  bool operator <  (const Dmer & m) const;
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
  Dmers(): m_mers(), m_dimCount(0), m_dmerLength(0), m_dmerCount(0) {}

  int NumMers()                            { return m_dmerCount;     }
  void BuildDmers(const RSiteReads& rReads , int dmerLength, int motifLength, int countPerDimension); 

private:
  void AddSingleReadDmers(const RSiteReads& rReads , int rIdx);
  inline int MapNToOneDim(const svec<int>& nDims);
  inline svec<int> MapOneToNDim(int oneDMappedVal);

  svec<svec<Dmer> > m_mers;  /// Multi-dimensional matrix representation of dmers projected on to dimensions
  int m_dimCount;            /// Number of cells in each dimension (this is dependent on the site values and the reduction coefficient)
  int m_dmerLength;          /// Number of dimensions in the matrix (i.e. dmer length)
  int m_dmerCount;           /// Total number of dmers
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

class RestSiteDataParams 
{
public:
  RestSiteDataParams( int totalNumReads=10000000, int meanReadLength=10000, int minAlignLength=4000, 
                      float deletionErr=0.03, float insertionErr=0.03, float substitutionErr=0.03)
                     :m_totalNumReads(totalNumReads), m_meanReadLength(meanReadLength), m_minAlignLength(minAlignLength),
                      m_deletionErr(deletionErr), m_insertionErr(insertionErr), m_substitutionErr(substitutionErr) { }

    bool   TotalNumReads() const     { return m_totalNumReads;   }
    int    MeanReadLength() const    { return m_meanReadLength;  }  
    int    MeanAlignLength() const   { return m_minAlignLength;  } 
    float  DeletionErr() const       { return m_deletionErr;     }
    float  InsertionErr() const      { return m_insertionErr;    }
    float  SubstitutionErr() const   { return m_substitutionErr; }

private: 
  bool    m_totalNumReads;    /// Flag specifying whether the reads are single or double strand
  int     m_meanReadLength;   /// Length of each motif
  int     m_minAlignLength;   /// Number of motifs to generate/use
  float   m_deletionErr;      /// The length of distmers to use for seed finding
  float   m_insertionErr;     /// The length of distmers to use for seed finding
  float   m_substitutionErr;  /// The length of distmers to use for seed finding
};

class RestSiteModelParams 
{
public:
  RestSiteModelParams(bool singleStrand=false, int motifSize=4, int numOfMotifs=20, 
                      int dmerLength=6, float cndfCoef=2.2, string alphabet="ACGT") 
                     :m_singleStrand(singleStrand), m_motifSize(motifSize), m_numOfMotifs(numOfMotifs),
                      m_dmerLength(dmerLength), m_cndfCoef(cndfCoef), m_alphabet(alphabet) { }

    bool   IsSingleStrand() const    { return m_singleStrand;   }
    int    MotifSize() const         { return m_motifSize;      }  
    int    NumOfMotifs() const       { return m_numOfMotifs;    } 
    int    DmerLength() const        { return m_dmerLength;     }
    float  CNDFCoef() const          { return m_cndfCoef;       }
    string Alphabet() const          { return m_alphabet;       }

private: 
  bool    m_singleStrand;   /// Flag specifying whether the reads are single or double strand
  int     m_motifSize;      /// Length of each motif
  int     m_numOfMotifs;    /// Number of motifs to generate/use
  int     m_dmerLength;     /// The length of distmers to use for seed finding
  float   m_cndfCoef;       /// Cumulative Normal Distribution Function coefficeint used for estimating similarity 
  string  m_alphabet;       /// Alphabet containing base letters used in the reads/motifs 
};

class RestSiteAlignUnit 
{
public:
  RestSiteAlignUnit(): m_rReads(), m_motifs() {}

  /* Generate Permutation of the given alphabet to reach number of motifs required */
  void GenerateMotifs(int motifLength, int numOfMotifs);  
  void CartesianPower(const vector<char>& input, unsigned k, vector<vector<char>>& result); 

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
