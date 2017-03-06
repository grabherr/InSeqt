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
  
  svec<int> & Dist() {return m_dist;}
  string & Name() {return m_name;}
  const svec<int> & Dist() const {return m_dist;}
  const string & Name() const {return m_name;}
  int & Ori() {return m_ori;}
  const int & Ori() const {return m_ori;}

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

class ORLink
{
public:
  ORLink(): m_forward(-1), m_back(-1), m_last(-1), m_first(-1), m_id(-1) {}

  int GetForward() const { return m_forward; }
  int GetBack() const    { return m_back;    }
  int GetLast() const    { return m_last;    }
  int GetFirst() const   { return m_first;   }
  int GetId() const      { return m_id;      }
  
  void SetForward(int fw) { m_forward = fw; }
  void SetBack(int bk)    { m_back    = bk; }
  void SetLast(int ls)    { m_last    = ls; }
  void SetFirst(int fs)   { m_first   = fs; }
  void SetId(int id)      { m_id      = id; }

private:
  int m_forward;
  int m_back;
  int m_last;
  int m_first;
  int m_id;
};

class ORLinks
{
public:
  ORLinks(): m_links() {}

  ORLinks(int capacity): m_links() {
    m_links.resize(capacity);
    // Reset all links
    for (int i=0; i<m_links.isize(); i++) {
      m_links[i].SetLast(i);
      m_links[i].SetFirst(i);
    }
  }

  const ORLink& operator[](int idx) const   { return m_links[idx];      }
  int NumLinks()                            { return m_links.isize();   }
private:
  svec<ORLink> m_links;  /// List containing all links
};

class OptiMapAlignUnit 
{
public:
  OptiMapAlignUnit(): m_reads(), m_optimers() {}
  
  void LoadReads(const string& fileName, int seedSize)  { m_reads.LoadReads(fileName, seedSize);       }
  void BuildOptimers(int seedSize)                      { m_optimers.BuildOptimers(m_reads, seedSize); }
  void PoolReads();

private:
  OptiReads m_reads;     /// Optical Reads
  Optimers  m_optimers;  /// Optimers built from optical reads
};
#endif //OPTIMAPALIGNUNIT_H
