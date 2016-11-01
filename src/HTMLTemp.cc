#include "HTMLTemp.h"
#include "Database.h"
#include "ryggrad/src/base/FileParser.h"


string Clean(const string & in)
{
  if (in.size() < 2)
    return in;
  if (in[0] == '.' && in[1] == '/') {
    string in2 = &in[2];
    return in2;
  }
  return in;  
}

string HTMLRead::GetRelativePath(const string & out_raw)
{
  if (m_relPath =="")
    return "";

  m_relPath = Clean(m_relPath);
  string out = Clean(out_raw);
  
  
  int i;

  cout << "Rel: " << m_relPath << endl;
  cout << "HTM: " << out << endl;
  
  StringParser p1;
  p1.SetLine(m_relPath, "/");
  StringParser p2;
  p2.SetLine(out, "/");

  int k = 0;
  int n = p1.GetItemCount();
  if (p2.GetItemCount() > n)
    n = p2.GetItemCount() ;

  string plus;
  
  for (i=0; i<n; i++) {
    if (i == p2.GetItemCount()-1)
      break;
    if (i == p1.GetItemCount())
      break;

    if (p1.AsString(i) == p2.AsString(i))
      continue;

    cout << p1.AsString(i) << " -- " <<  p2.AsString(i) << endl;
    plus += p1.AsString(i) + "/";
    plus = "../" + plus;
    cout << "Added " << plus << endl;
    
  }
  int j;
  for (j=i; j<p1.GetItemCount(); j++)
    plus += p1.AsString(j) + "/";
  for (j=i; j<p2.GetItemCount()-1; j++)
    plus = "../" + plus;

  return plus;
}

void HTMLRead::Read(const string & fileName, const string & delim)
{
  if (delim != "")
    m_delim = delim;

  FlatFileParser parser;
  
  parser.Open(fileName);
  int i;
  bool bTable = false;
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;

    //cout << "Reading template " << parser.Line() << endl;

    StringParser pp;
    pp.SetLine(parser.Line(), m_delim);
    //cout << "Parts " << pp.GetItemCount() << endl;
    if (pp.GetItemCount() > 1 && pp.AsString(1) == "tablerow_begin") {
      //cout << "############ Prepare table!!" << endl;
      bTable = true;
      HTMLPart tabins;
      tabins.SetToken("internal_table");
      m_parts.push_back(tabins);
      continue;

    }
    if (pp.GetItemCount() > 1 && pp.AsString(1) == "tablerow_end") {
      bTable = false;
      continue;
    }

    for (i=0; i<pp.GetItemCount(); i++) {
      bool br = false;
      if (i+1 == pp.GetItemCount())
	br = true;
      HTMLPart t;
      if (i % 2 == 0) {
        t.SetData(pp.AsString(i), br);       	
      } else {
	t.SetToken(pp.AsString(i));
      }
      if (bTable) 
	m_table.push_back(t);
      else 
	m_parts.push_back(t);
	
    }

  }


}
void HTMLRead::FillWrite(const Database & db, const string & fileName)
{
 
  
  int i, j, k;
  FILE * pOut = fopen(fileName.c_str(), "w");

  string imgPath = GetRelativePath(fileName);

  cout << "Image path: " << imgPath << endl;

  //cout << "FillWrite" << endl;
  for (i=0; i<m_parts.isize(); i++) {
    //cout << "Getting part " << i << endl;
    const HTMLPart & p = m_parts[i];
    //cout << "Got part " << i << endl;
    if (p.Token() == "internal_table") {
      //cout << "Inserting table!" << endl;
      for (j=0; j<db.GetLibCount(); j++) {
	const KeyValueSet & t_db = db.GetLib(j);
	//cout << "Process " << j << endl;
	for (k=0; k<m_table.isize(); k++) {
	  const HTMLPart & t = m_table[k];	
	  if (t.Token() == "") {
	    //cout << "Print " << t.Data() << endl;
	    fprintf(pOut, "%s", t.Data().c_str());
	    if (t.HasBR())
	      fprintf(pOut, "\n");
	  } else {
	    string rep = t_db.Get(t.Token());
	    //cout << "Fill " << t.Token() << " -> " << rep << endl;
	    fprintf(pOut, "%s", rep.c_str());      
	  }
	}
      }
      continue;
    }

    if (p.Token() == "") {
      //cout << "Write " << p.Data() << endl;

      string parsed = p.Data();
      if (imgPath != "") {
	StringParser rep;
	rep.SetLine(parsed, "img src=\"");
	string abs;
	//cout << "Len " << rep.GetItemCount() << endl;
	for (int x=0; x<rep.GetItemCount(); x++) {
	  if (x > 0) {
	    abs += "img src=\"";
	    abs += imgPath;
	  }
	  abs += rep.AsString(x);
	  // if (x > 0) {
	  //cout << "Insert " << abs << endl;
	  // cout << "Before " << rep.AsString(x) << endl;
	  //}
	}
	parsed = abs;
      }

      
      fprintf(pOut, "%s", parsed.c_str());
      if (p.HasBR())
	fprintf(pOut, "\n");
    } else {
      string rep = db.Get(p.Token());
      //cout << "Replace " << p.Token() << " --> " << rep << endl;
      fprintf(pOut, "%s", rep.c_str());      
    }
  }
  fclose(pOut);
}

