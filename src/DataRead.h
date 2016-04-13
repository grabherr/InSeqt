#ifndef DATAREAD_H
#define DATAREAD_H

#include "src/DNAVector.h"
#include "base/FileParser.h"


inline void ReadFileNames(string & fileName)
{
  FlatFileParser parser;
  
  parser.Open(fileName);

  fileName = "";

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    for (int i=0; i<parser.GetItemCount(); i++) {
      if (fileName != "")
	fileName += ",";
      fileName += parser.AsString(i);
    }
  }
}


inline void ReadDNA(vecDNAVector & dna, const string & fileName)
{
  FlatFileParser parser;
  
  parser.Open(fileName);
  parser.ParseLine();

  string input = fileName;

  if (parser.AsString(0)[0] != '@' && parser.AsString(0)[0] != '>') {
    ReadFileNames(input);
  }

  dna.Read(input);
}


#endif

