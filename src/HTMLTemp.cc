#include "HTMLTemp.h"
#include "Database.h"
#include "ryggrad/src/base/FileParser.h"


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
      fprintf(pOut, "%s", p.Data().c_str());
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

