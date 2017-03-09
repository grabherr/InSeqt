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

class OptiMapAlignUnit 
{
public:
  OptiMapAlignUnit(): m_reads() {}
  
  void LoadReads(const string& fileName, int seedSize)  { m_reads.LoadReads(fileName, seedSize);       }
  void FindCandidLaps(int seedSize);

private:
  OptiReads m_reads;     /// Optical Reads
};
#endif //OPTIMAPALIGNUNIT_H
