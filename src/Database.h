#ifndef DATABASE_H
#define DATABASE_H

#include "ryggrad/src/base/SVector.h"

class KeyValue
{
 public:
  void Set(const string & key, const string & val) {
    m_key = key;
    m_val = val;
  }
  void Set(const string & key, double val) {
    m_key = key;
    char tmp[256];
    sprintf(tmp, "%f", val);
    m_val = tmp;
  }
  void Set(const string & key, int val) {
    m_key = key;
    char tmp[256];
    sprintf(tmp, "%d", val);
    m_val = tmp;
  }
  
  const string & Key() const {return m_key;}
  const string & Val() const {return m_val;}

 private:
  string m_val;
  string m_key;
};


class KeyValueSet
{
 public:
  string Get(const string & key) const {
    int i;
    for (i=0; i<m_data.isize(); i++) {
      if (key == m_data[i].Key())
	return m_data[i].Val();
    }
    string empty = "n/a";
    return empty;
  }

  void Add(const KeyValue & k) {
    m_data.push_back(k);
  }
  void Add(const string & key, const string & val) {
    KeyValue v;
    v.Set(key, val);
    m_data.push_back(v);
  }

  void Add(const string & key, int val) {
    KeyValue v;
    v.Set(key, val);
    m_data.push_back(v);
  }

  void Add(const string & key, double val) {
    KeyValue v;
    v.Set(key, val);
    m_data.push_back(v);
  }
  svec<KeyValue> & Data() {return m_data;}

 private:
  svec<KeyValue> m_data;
};


class Database
{
 public:
  

  int GetLibCount() const {return m_lib.isize();}
  const KeyValueSet & GetLib(int i) const {
    return m_lib[i];
  }
  void Add(const KeyValueSet & s) {
    m_lib.push_back(s);
  }


  string Get(const string & key) const {
    return m_global.Get(key);
  }
  KeyValueSet & Global() {return m_global;}


  void Add(const string & key, const string & val) {
    m_global.Add(key, val);
  }
  
  void Add(const string & key, int val) {
    m_global.Add(key, val);
  }
  
  void Add(const string & key, double val) {
    m_global.Add(key, val);
  }


 private:
  KeyValueSet m_global;
  svec<KeyValueSet> m_lib;

};

void ReadDB(Database & d, const string & fileName);

void ReadDBPlain(Database & d, const string & fileName);


#endif //DATABASE_H
