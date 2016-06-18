#include "src/Database.h"
#include "base/FileParser.h"


void Read(Database & d, const string & fileName)
{

  FlatFileParser parser;
  
  parser.Open(fileName);
  int i;
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;

    if (parser.AsString(0) == "SMRT" && parser.AsString(1) != "total") {      
      KeyValueSet n;
      n.Add("name", parser.AsString(1));
      parser.ParseLine();
      n.Add("readcount", parser.AsInt(1));
      parser.ParseLine();
      n.Add("average", parser.AsFloat(1));
      parser.ParseLine();
      n.Add("median", parser.AsInt(1));
      parser.ParseLine();
      n.Add("sequence", parser.AsFloat(1));
      parser.ParseLine();
      n.Add("n50", parser.AsInt(1));
      d.Add(n);
    }
    if (parser.AsString(0) == "SMRT" && parser.AsString(1) == "total") {      
      KeyValueSet & n = d.Global();
      n.Add("readcount_total", parser.AsInt(1));
      parser.ParseLine();
      n.Add("average_total", parser.AsFloat(1));
      parser.ParseLine();
      n.Add("median_total", parser.AsInt(1));
      parser.ParseLine();
      n.Add("sequence_total", parser.AsFloat(1));
      parser.ParseLine();
      n.Add("n50_total", parser.AsInt(1));      
    }
    if (parser.AsString(0) == "#") {      
      KeyValueSet & n = d.Global();
      n.Add(parser.AsString(1), parser.AsInt(2));
    }

  }
}
