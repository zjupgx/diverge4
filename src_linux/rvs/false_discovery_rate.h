#ifndef _FALSE_DISCOVERY_RATE_2_H_
#define _FALSE_DISCOVERY_RATE_2_H_

#include <string>
#include <vector>

#include "tree.h"
#include "sequence.h"


bool false_discovery_rate_compute(const std::vector<Tree> &trees, 
		  const std::vector<sequence_t> &sequences, 
          std::vector<std::vector<double> > &rets1, 
	      std::vector<std::vector<double> > &rets2, 
		  std::vector<std::vector<double> > &rets3, 
		  std::vector<std::vector<double> > &rets4); 


double GAMMA(double); 

int trans2(int x); 
const char * trans2Str(int x); 
const char * trans16(int x); 

char * strcat(char * s1, char * s2); 


#endif

