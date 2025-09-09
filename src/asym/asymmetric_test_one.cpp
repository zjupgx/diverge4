#define _CRT_SECURE_NO_WARNINGS
#define WIN32

#include <string>
#include <vector>

#include <stdio.h>
#include <stdlib.h>

#include "asymmetric_test_one.h"
#include "gz97.h"
#include "GZf2.h"

#include "tree.h"
#include "matrices.h"
#include "matrix.h"
#include <math.h> 





#include <algorithm>

#include <iostream>
#ifdef WIN32
#include <minmax.h>
#endif

using namespace std;
//----------------------------------------------------------------------

extern double p0;
extern DVector site_muta;

bool asymmetric_test_one_compute(const vector<Tree> &trees, const vector<sequence_t> &sequences,
	     vector<vector<double> > &rets) 
{
  
//  Type One Computation. 
	int ntrees = trees.size();

	int seqLength = sequences[0].sequence.size();
	vector<DVector> retsV(ntrees);

	DVector freq(jtt_freq, 20);
	DMatrix2D prob(jtt_prob, 20, 20);

	vector< DVector > site_muta_groups;
	site_muta_groups.clear();


	for(int i = 0; i < ntrees; i++) 
	{
		int treeE;

		

		vector<string> taxa;
		trees[i].leaf_names(taxa);

		CMatrix2D alignment(taxa.size(), sequences[0].sequence.size());

		{
			// generate sub_sequences and sub_seqNames
			vector<string>::const_iterator i;
			vector<sequence_t>::const_iterator j;
			int i2;
			for(i2 = 0, i = taxa.begin(); i != taxa.end(); i++) {
				for(j = sequences.begin(); j != sequences.end(); j++) {
					if(j->label == *i) {
						for(int k=0; k<(int)j->sequence.size(); k++) {
							alignment(i2, k) = j->sequence[k];
						}
						i2++;
						break;
					}
				}

				if(j == sequences.end()) {
					abort();
				} 
			}
		}

		string tree_str;
		if(!trees[i].gen_str_wrt_seq(taxa, tree_str)) 
			return false;

		IMatrix2D gu_tree;

		if(!parse_tree(tree_str.c_str(), alignment.rows(), treeE, gu_tree)) {
			return false;
		}
		if(!gz97_compute(true, 2.4, alignment, treeE, gu_tree, freq, prob, retsV[i])) {
			return false;
		}

		site_muta_groups.push_back(site_muta);

	}

	double p000=0; 

	int iZero=0;
	for(int i=0;i<seqLength;i++)
	{
		if(site_muta_groups[0](i)<0.001 && site_muta_groups[1](i)<0.001 && site_muta_groups[2](i)<0.001) iZero++;
	}


	vector<DVector> &numSub = retsV; 
	vector<vector<double> > summary; 
	vector<vector<double> > rets1; 

    if(!GZf2_compute(numSub, summary, rets1)) 
    	return false; 


	int n = 3; 

	std::vector<double> theta0A; 
	std::vector<double> solvedThetaA; 

	std::vector<double> coefficientA; 

	for (int i = 0; i < n; i++)
	{
		theta0A.push_back(0.50); 
	}

	// sum[0] = theta,    sum[1] = se
	// sum[2] = r_X,      sum[3] = r_max,
	// sum[4] = z_score   (sum[0-4] are model-free
	// sum[5] = thetaML,  sum[6] = alphaML,
	// sum[7] = se_theta, sum[8] = LRT

  coefficientA.push_back(summary[0][0]); 
	coefficientA.push_back(summary[1][0]); 
	coefficientA.push_back(summary[2][0]); 

	double d12 = 0.0; 
	double d13 = 0.0; 
	double d23 = 0.0; 

	double d1 = 0.0; 
	double d2 = 0.0; 
	double d3 = 0.0; 

	d12 = log(1 - summary[0][0]); 
	d13 = log(1 - summary[1][0]); 
	d23 = log(1 - summary[2][0]); 

	d1 = (d12 + d13 - d23) / 2; 
	d2 = (d12 + d23 - d13) / 2;
	d3 = (d13 + d23 - d12) / 2; 

	solvedThetaA.push_back(1 - exp(d1)); 
	solvedThetaA.push_back(1 - exp(d2)); 
	solvedThetaA.push_back(1 - exp(d3)); 


	for (int i = 0; i < n; i++)
	{
		if (solvedThetaA[i] < 0)
		{
			solvedThetaA[i] = 0.0; 
		}
	}
	
	std::vector<double> varianceDeltaA; 

	// Pay attention to the following values of vector named by varianceDeltaA. 
	// Three values are respectively according to the corresponding outgroup chosen from Cluster 1, Cluster 2 or Cluster 3. 

	varianceDeltaA.push_back(summary[0][1] * summary[0][1] + summary[1][1] * summary[1][1] - 2 * solvedThetaA[2] * (summary[0][1] * summary[1][1])); 
	varianceDeltaA.push_back(summary[0][1] * summary[0][1] + summary[2][1] * summary[2][1] - 2 * solvedThetaA[1] * (summary[0][1] * summary[2][1])); 
	varianceDeltaA.push_back(summary[1][1] * summary[1][1] + summary[2][1] * summary[2][1] - 2 * solvedThetaA[0] * (summary[1][1] * summary[2][1])); 
	

	vector<double> retsLine; 

	retsLine.clear(); 
	retsLine.push_back(solvedThetaA[2])    ; 
	retsLine.push_back(varianceDeltaA[0]); 
	rets.push_back(retsLine); 
	retsLine.clear(); 
	retsLine.push_back(solvedThetaA[1]); 
	retsLine.push_back(varianceDeltaA[1]); 
	rets.push_back(retsLine); 
	retsLine.clear(); 
	retsLine.push_back(solvedThetaA[0]); 
	retsLine.push_back(varianceDeltaA[2]); 
	rets.push_back(retsLine); 
	
  
  return true; 

}



/*----------------------------------------------------------------------*/


double mean(const std::vector<double> X, int n) 
{
	double meanValue = 0.0; 

	for (int i = 0; i < n; i++)
	{
		meanValue += X[i]; 
	}
	meanValue /= n; 

	return meanValue; 
}


/*--------------------------------------------------------------------*/


double variance(const std::vector<double> X, int n) 
{
	double meanValue = 0.0; 
	double varianceValue= 0.0; 
	
	meanValue = mean(X, n); 

	for (int i = 0; i < n; i++)
	{
		varianceValue += (X[i] - meanValue) * (X[i] - meanValue); 
	}
	varianceValue /= n; 

	return varianceValue; 
}


/*--------------------------------------------------------------------*/


int solveNonlinearEquation_1(std::vector<double> thetaN1A, std::vector<double> &thetaN2A, 
	                         const std::vector<double> coefficientA, const int n) 
{
	int rtnCode = 0; 
	
	if (n != 3)
	{
		printf("The number of variables of equations is not Three !");
		return rtnCode; 
	}

	double v1 = 0.0; 
	double v2 = 0.0; 
	double v3 = 0.0; 
	double v4 = 0.0; 
	double v5 = 0.0; 
	double v6 = 0.0; 
	double v7 = 0.0; 
	double v8 = 0.0; 
	double v9 = 0.0; 
	double v10 = 0.0; 

	double solutionDist = 100.0; 

	int maxIterNum = 50000; 
	int i = 0; 

	for (int j = 0; j < n; j++)
	{
		thetaN2A.push_back(thetaN1A[j]); 
	}

	while (solutionDist > 1e-4 && i < maxIterNum) 
	{
		for (int j = 0; j < n; j++)
		{
			thetaN1A[j] = thetaN2A[j]; 
		}

		v1 = 1.0 / (2 * (thetaN1A[0] - 1) * (thetaN1A[1] - 1) * (thetaN1A[2] - 1)); 
		v2 = - (thetaN1A[0] - 1) * (thetaN1A[0] - 1); 
		v3 = (thetaN1A[0] - 1) * (thetaN1A[1] - 1); 
		v4 = - (thetaN1A[1] - 1) * (thetaN1A[1] - 1); 
		v5 = (thetaN1A[0] - 1) * (thetaN1A[2] - 1); 
		v6 = (thetaN1A[1] - 1) * (thetaN1A[2] - 1); 
		v7 = - (thetaN1A[2] - 1) * (thetaN1A[2] - 1); 
		v8 = (thetaN1A[0] - 1) * (thetaN1A[1] - 1) - (thetaN1A[0] - 1) - (thetaN1A[1] - 1) + coefficientA[0]; 
		v9 = (thetaN1A[0] - 1) * (thetaN1A[2] - 1) - (thetaN1A[0] - 1) - (thetaN1A[2] - 1) + coefficientA[1]; 
		v10 = (thetaN1A[1] - 1) * (thetaN1A[2] - 1) - (thetaN1A[1] - 1) - (thetaN1A[2] - 1) + coefficientA[2]; 

		thetaN2A[0] = thetaN1A[0] - (v5 * v8 + v3 * v9 + v2 * v10) / v1; 
		thetaN2A[1] = thetaN1A[1] - (v6 * v8 + v4 * v9 + v3 * v10) / v1; 
		thetaN2A[2] = thetaN1A[2] - (v7 * v8 + v6 * v9 + v5 * v10) / v1; 

		solutionDist = 0.0; 
		for (int j = 0; j < n; j++)
		{
			solutionDist += (thetaN1A[j] - thetaN2A[j]) * (thetaN1A[j] - thetaN2A[j]); 
		}
		solutionDist = sqrt(solutionDist); 

		i++; 
	}

	if (solutionDist <= 1e-4) 
	{
		rtnCode = 1; 
	}
	return rtnCode; 

}


/*--------------------------------------------------------------------*/







