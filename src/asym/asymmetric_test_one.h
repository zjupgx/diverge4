#ifndef _ASYMMETRIC_TEST_ONE_H_
#define _ASYMMETRIC_TEST_ONE_H_

#include <string>
#include <vector>

#include "tree.h"
#include "sequence.h"



bool asymmetric_test_one_compute(const std::vector<Tree> &trees,
		  const std::vector<sequence_t> &sequences,
		  std::vector<std::vector<double> > &rets); 


double GAMMA(double); 

int trans2(int x); 
const char * trans2Str(int x); 
const char * trans16(int x); 

char * strcat(char * s1, char * s2); 


double mean(const std::vector<double> X, int n); 
double variance(const std::vector<double> X, int n); 

int solveNonlinearEquation_1(std::vector<double> thetaN1A, std::vector<double> &thetaN2A, 
	                         const std::vector<double> coefficientA, const int n); 


#endif

