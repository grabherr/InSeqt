#include <string>

#include "ryggrad/src/base/CommandLineParser.h"
#include "ryggrad/src/general/DNAVector.h"
#include "aligns/KmerAlignCore.h"
#include "SequenceMatch.h"
#include "DataRead.h"
#include "ryggrad/src/base/RandomStuff.h"
#include <math.h>
#include "ryggrad/src/base/StringUtil.h"

void MakeRandomList(svec<int> & l, int n)
{
  l.resize(n, 0);
  int i;
  for (i=0; i<l.isize(); i++) {
    l[i] = i;
  }
  for (i=0; i<l.isize(); i++) {
    int swap = RandomInt(n);
    int tmp = l[i];
    l[i] = l[swap];
    l[swap] = tmp;
  }
  
}


class SingleHitNoPos
{
public:
  SingleHitNoPos() {
    m_contig = -1;
    m_pos = -1;
    m_rc = false;
  }

  SingleHitNoPos(int contig, int pos, bool rc) {
    m_contig = contig;
    m_pos = pos;
    m_rc = rc;
  }

  bool operator < (const SingleHitNoPos & h) const {
    return m_contig < h.m_contig;
  }

  int Contig() const {return m_contig;}
  int Pos() const {return m_pos;}
  bool RC() const {return m_rc;}

private:
  int m_contig;
  int m_pos;
  bool m_rc;
};

class SingleHit
{
public:
  SingleHit() {
    m_contig = -1;
    m_pos = -1;
    m_qPos = -1;
    m_rc = false;
  }

  SingleHit(int contig, int pos, int qPos, bool rc) {
    m_contig = contig;
    m_pos = pos;
    m_qPos = qPos;
    m_rc = rc;
  }

  bool operator < (const SingleHit & h) const {
    if (m_rc != h.m_rc)
      return h.m_rc;
    if (m_contig != h.m_contig)
      return m_contig < h.m_contig;
    return m_pos < h.m_pos;
  }

  int Contig() const {return m_contig;}
  int Pos() const {return m_pos;}
  bool RC() const {return m_rc;}
  char RCC() const {
    if (m_rc)
      return '-';
    else
      return '+';
  }

private:
  int m_contig;
  int m_pos;
  int m_qPos;
  bool m_rc;
};


int main(int argc,char** argv)
{

  
  commandArg<string> aStringCmmd("-i","input sequence file or file list");
  commandArg<string> outCmmd("-o","output directory");
  commandArg<int> nCmmd("-n","chunk index (0-based)", 0);
  commandArg<int> chunksCmmd("-c","total chunks", 1);
  commandArg<int> distCmmd("-distance","distance between seeds", 1);
  commandArg<int> numCmmd("-w","width of the filter", 2);
  commandArg<double> fracCmmd("-f","fraction of reads to examine (<=1.)", 1.0);
  commandLineParser P(argc,argv);
  P.SetDescription("Finds overlaps between sequences.");
  P.registerArg(aStringCmmd);
  P.registerArg(outCmmd);
  P.registerArg(nCmmd);
  P.registerArg(chunksCmmd);
  P.registerArg(distCmmd);
  P.registerArg(numCmmd);
  P.registerArg(fracCmmd);

  P.parse();
 
  int chunks = P.GetIntValueFor(chunksCmmd);
  int thisone = P.GetIntValueFor(nCmmd);

  string fileName = P.GetStringValueFor(aStringCmmd);
  string outName = P.GetStringValueFor(outCmmd);
  outName += "/overlapcands.out.";
  outName += Stringify(thisone);
 
  int distance = P.GetIntValueFor(distCmmd);
  int num12 = P.GetIntValueFor(numCmmd);
  double fraction = P.GetDoubleValueFor(fracCmmd);
  
 
  vecDNAVector dna;
   
  ReadDNA(dna, fileName);
 

  int i, j, l;

  svec<int> order;
  MakeRandomList(order, dna.isize());

  int upto = (int)((double)dna.isize()*fraction);
  cout << "Will examine " << upto << " sequences." << endl;


  // Set up query
  vecDNAVector query;  
  for (i=0; i<upto; i++) {
    query.push_back(dna[order[i]], dna.Name(order[i]));
  }
  
  // There is a little bit of a round-off error here...
  int size = dna.isize()/chunks;
  int from = thisone*size;
  int to = from + size;

  cout << "Only keep sequences from " << from << " to " << to << " of " << dna.isize() << endl;
  for (i=0; i<from; i++)
    dna[i].clear();
  for (i=to; i<dna.isize(); i++)
    dna[i].clear();


  KmerAlignCore<KmerAlignCoreRecordWithScore> core;
  core.SetNumTables(num12);
  int k = 12;
  TranslateBasesToNumberExact trans;
  trans.SetSize(k); // Use 15-mers
  core.SetTranslator(&trans);
  
  cout << "Adding target." << endl;
  core.AddData(dna, false, 1);
  core.SortAll();
  cout << "done" << endl;

  string outNamePlus = outName + ".cand";
  FILE * pOut = fopen(outNamePlus.c_str(), "w");
  string done = outName + ".done";
  string remove = "rm " + done;
  int retRem = system(remove.c_str());
  outName += ".allreadnames.";
  outName += Stringify(thisone);
  FILE * pAll = fopen(outName.c_str(), "w");
  fprintf(pOut, "Seq1\tSeq2\tLen1\tLen2\tOri\t#Kmers\n");
  
  

  int maxSingle = 100;
 
  for (int i=0; i<query.isize(); i++) {   
    const DNAVector & d = query[i];
    //cout << query.Name(i) << endl;
    svec<SingleHit> hits;
    int min_k = d.isize()/k/20;
    fprintf(pAll, "%s\n", 
	    query.Name(i).c_str());
    fflush(pAll);

    //cout << "Start." << endl;
    int maxHits = d.isize()*100;
    for (j=0; j<=d.isize()-k*num12; j+=distance) {
      //cout << "Lookup " << j << endl;
      svec<SingleHitNoPos> tmp;
      DNAVector sub;
      sub.SetToSubOf(d, j, k*num12);
      svec<KmerAlignCoreRecordWithScore> matches;
      //cout << "Before." << endl;
      core.GetMatches(matches, sub);
      //cout << "After." << endl;
      int found = 0;

      int l;
      for (l=0; l<matches.isize(); l++) {
	//cout << "Match w/ " << matches[l].GetContig() << endl;
	if (matches[l].GetContig() != order[i])
	  tmp.push_back(SingleHitNoPos(matches[l].GetContig(), matches[l].GetPosition(), false));
      }

      
      svec<KmerAlignCoreRecordWithScore> matches_rc;     
      sub.ReverseComplement();
      core.GetMatches(matches_rc, sub);
      int n = matches.isize()+matches_rc.isize();
      //  if (n >1)
      // cout << "Matches: " << n << " " << log((double)n) << endl;
      //cout << "Matches: " << matches.isize() << " " << matches_rc.isize() << endl;
      int fw_count = tmp.isize();
      //cout << i << endl;
      //cout << order[i] << endl;
      for (l=0; l<matches_rc.isize(); l++) {
	if (matches_rc[l].GetContig() != order[i])
	  tmp.push_back(SingleHitNoPos(matches_rc[l].GetContig(), matches_rc[l].GetPosition(), true));
      }
      //cout << "Unique" << endl;
      //UniqueSort(tmp);
      //cout << j << " -> " << tmp.isize() << endl;
      //if (tmp.isize() > 0) {
      //cout << "    pos: " << tmp[0].Pos() << endl;
      //}
      if (tmp.isize() < maxSingle) {
	//cout << tmp.isize() << " " << fw_count << endl;
	for (l=0; l<fw_count; l++) {
	  hits.push_back(SingleHit(tmp[l].Contig(), tmp[l].Pos(), j, tmp[l].RC()));
	}
      
      
	for (; l<tmp.isize(); l++) {
	  hits.push_back(SingleHit(tmp[l].Contig(), tmp[l].Pos(), d.isize()-j, tmp[l].RC()));
	}
      }
      if (hits.isize() > maxHits)
	break;
     
    }
    //cout << "Done" << endl;
    if (hits.isize() > maxHits)
      continue;

    UniqueSort(hits);
    //cout << "After unique." << endl;
    int kmers = 0;
    //cout << "Hits: " << hits.isize() << endl;
    for (j=1; j<hits.isize(); j++) {
      //cout << j << " Contig " << hits[j].Contig() << " Pos " << hits[j].Pos() << endl;
      if (hits[j].Contig() == hits[j-1].Contig()) {
	if (hits[j].Pos() > hits[j-1].Pos())
	  kmers++;
      } else {
	if (kmers > min_k) {
	  cout << query.Name(i) << " overlaps w/ " << dna.Name(hits[j-1].Contig()) << " " << kmers << endl;
	  fprintf(pOut, "%s\t%s\t%d\t%d\t%c\t%d\n", 
		  query.Name(i).c_str(), dna.Name(hits[j-1].Contig()).c_str(), 
		  query[i].isize(), dna[hits[j-1].Contig()].isize(), hits[j-1].RCC(), kmers);
	  fflush(pOut);
	}
	kmers = 0;
      }
    }
    if (kmers > min_k) {
      cout << dna.Name(i) << " overlaps w/ " << dna.Name(hits[j-1].Contig()) << " " << kmers << endl;
      fprintf(pOut, "%s\t%s\t%d\t%d\t%c\t%d\n", 
	      query.Name(i).c_str(), dna.Name(hits[j-1].Contig()).c_str(), 
	      query[i].isize(), dna[hits[j-1].Contig()].isize(), hits[j-1].RCC(), kmers);
    }
    //cout << "Done." << endl;
  }

  fclose(pOut);
  fclose(pAll);

  FILE * pDone = fopen(done.c_str(), "w");
  fprintf(pDone, "done\n");
  fclose(pDone);

  return 0;

}
  
