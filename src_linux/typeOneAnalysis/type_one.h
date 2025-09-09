#ifndef _TYPE_ONE_H_
#define _TYPE_ONE_H_

#include <string>
#include <vector>

#include "tree.h"
#include "sequence.h"


//#pragma once


bool type_one_compute(const std::vector<Tree> &trees,
		  const std::vector<sequence_t> &sequences,
		  std::vector<std::vector<double> > &summary,
		  std::vector<std::vector<double> > &rets2); 


double GAMMA(double); 

int trans2(int x); 
const char * trans2Str(int x); 
const char * trans16(int x); 

char * strcat(char * s1, char * s2); 


#endif

