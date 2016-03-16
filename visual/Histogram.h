#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include "base/SVector.h"
#include "visual/Whiteboard.h"



class Histogram
{
 public:
  Histogram() {}

  void Plot(const string & fileName,
	    const svec<double>&data, 
	    int bins, 	    
	    color c = color(0., 0., 0.7));
  
  void Plot(const string & fileName,
	    const svec<double>&data, 
	    int bins, 	    
	    double lo,
	    double hi,
	    color c = color(0., 0., 0.7));

  void Scatter(const string & fileName,
	    const svec<double>&data1, 
	    const svec<double>&data2, 	  
	    color c = color(0., 0., 0.7));
 private:
  

};


#endif

