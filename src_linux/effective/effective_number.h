#ifndef _EFFECTIVE_NUMBER_SITES_H_
#define _EFFECTIVE_NUMBER_SITES_H_

#include <string>
#include <vector>

#include "tree.h"
#include "sequence.h"


bool effective_number_compute(const std::vector<Tree> &trees,
		  const std::vector<sequence_t> &sequences,
		  std::vector<std::vector<double> > &rets1,
		  std::vector<std::vector<double> > &rets2, 
		  int &EffectiveNumberType1, int &EffectiveNumberType2); 


double GAMMA(double); 

int trans2(int x); 
const char * trans2Str(int x); 
const char * trans16(int x); 

char * strcat(char * s1, char * s2); 


#endif

