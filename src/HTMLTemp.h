#ifndef HTMLTEMP_H
#define HTMLTEMP_H

#include "base/SVector.h"

class HTMLPart
{
 public:
  HTMLPart() {
    m_br = false;
  }

  const string & Token() const {return m_token;}
  const string & Data() const {return m_data;}
  bool HasBR() const {return m_br;}

  void SetData(const string & d, bool br) {
    m_data = d;
    m_br = br;
  }
  void SetToken(const string & d) {
    m_token = d;
  }

 private:
  string m_token;
  string m_data;
  bool m_br;
};

class Database;

class HTMLRead
{
 public:
  void SetDelim(const string & d) {
    m_delim = d;
  }

  void Read(const string & fileName, const string & delim = "");
  void FillWrite(const Database & db, const string & fileName);

 
 private:
  string m_delim;
  svec<HTMLPart> m_parts;
  svec<HTMLPart> m_table;

};




#endif //HTMLTEMP_H



