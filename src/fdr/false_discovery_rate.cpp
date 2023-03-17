#define _CRT_SECURE_NO_WARNINGS
#define WIN32
#include <string>
#include <vector>

#include <stdio.h>
#include <stdlib.h>

#include "false_discovery_rate.h"
#include "gz97.h"
#include "GZf2.h"
#include "type_two.h"

#include "tree.h"
#include "matrices.h"
#include "matrix.h"


#include <math.h> 

#include <math.h>

#include <iostream>

//----------------------------------------------------------------------

using namespace std;

//----------------------------------------------------------------------

extern double p0;
DVector site_muta;

bool false_discovery_rate_compute(const vector<Tree> &trees, const vector<sequence_t> &sequences,
	     vector<vector<double> > &rets1, vector<vector<double> > &rets2, 
		 vector<vector<double> > &rets3, vector<vector<double> > &rets4) 
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
	for( int i = 0; i < seqLength; i++)
	{
		bool site_muta_groups_vbl = true; 

		for (int j = 0; j < ntrees; j++)
		{
			if (site_muta_groups_vbl && site_muta_groups[0](j) >= 0.001) 
			{
			    site_muta_groups_vbl = false; 
			}
		}

		if (site_muta_groups_vbl)
		{
			iZero++; 
		}
	}


	vector<DVector> &numSub = retsV; 
	vector<vector<double> > summary; 
	vector<vector<double> > rets; 

    if(!GZf2_compute(numSub, summary, rets)) 
    	return false; 


//  Sort the array of posterior probability in rets. 

	vector<vector<double> > sortedRets; 

	sortedRets = rets; 

	for( int i = 0; i < sortedRets.size(); i++)
	{

		double tempValue = 0.0; 

		for(int j = 0; j < sortedRets[i].size() - 1; j++)
		{
			for(int k = j + 1; k < sortedRets[i].size(); k++)
			{
				if (sortedRets[i][j] > sortedRets[i][k])
				{
					tempValue = sortedRets[i][j]; 
					sortedRets[i][j] = sortedRets[i][k]; 
					sortedRets[i][k] = tempValue; 
				}
			}
		}

	}


	double initialProb = 0.00; 
	double terminalProb = 1.00; 
	double stepProb = 0.02; 

	double cutoff = 0.0; 

	double sumProb = 0.0; 
	double FDR = 0.0; 

	int Lc = 0; 
	int n = 0; 


	cutoff = initialProb; 

	vector<double> cutoffV; 

	while (cutoff <= terminalProb + stepProb / 100)
	{
		cutoffV.push_back(cutoff); 
		cutoff += stepProb; 
	}
	rets1.push_back(cutoffV); 

	cutoffV.clear(); 
	

	for (int  i = 0; i < sortedRets.size(); i++)
	{
		
		cutoff = initialProb; 
		n = 0; 

		vector<double> FDRV; 

		while (cutoff <= terminalProb + stepProb / 100)
		{
			
			sumProb = 0.0; 
			FDR = 0.0; 

			Lc = 0; 

			for (int j = sortedRets[i].size() - 1; j >= 0; j--)
			{
				if (sortedRets[i][j] >= cutoff)
				{
					sumProb += sortedRets[i][j]; 
					Lc++; 
				}
				else
				{
				    break; 
				}
			}

			if (Lc != 0)
			{
				FDR = 1 - sumProb / Lc; 
			}

			FDRV.push_back(FDR); 

			cutoff += stepProb; 
			n++; 
		}

		rets1.push_back(FDRV); 


		FDRV.clear(); 

		sumProb = 0.0; 
		FDR = 0.0; 

		Lc = 0; 

		FDRV.push_back(FDR); 

		for (int j = sortedRets[i].size() - 1; j >= 0; j--)
		{
			sumProb += sortedRets[i][j]; 
			Lc++; 

			FDR = 1 - sumProb / Lc; 
			FDRV.push_back(FDR); 
		}

		for (  int j = sortedRets[i].size() - 2; j >= 0; j--)
		{
			if (sortedRets[i][j] == sortedRets[i][j + 1])
			{
				FDRV[sortedRets[i].size() - j - 1] = -1.0; 
			}
		}

		rets3.push_back(FDRV); 

	}
	


//  Type Two Computation. 

//	vector<vector<double> > summary; 
//	vector<vector<double> > rets; 
	rets.clear(); 

	if (!type_two_compute(trees, sequences, summary, rets)) 
	{
		return false;
	}


	//  Sort the array of posterior probability in rets. 

//	vector<vector<double> > sortedRets; 
	sortedRets.clear(); 

	sortedRets = rets; 

	for( int i = 0; i < sortedRets.size(); i++)
	{

		double tempValue = 0.0; 

		for(int j = 0; j < sortedRets[i].size() - 1; j++)
		{
			for(int k = j + 1; k < sortedRets[i].size(); k++)
			{
				if (sortedRets[i][j] > sortedRets[i][k])
				{
					tempValue = sortedRets[i][j]; 
					sortedRets[i][j] = sortedRets[i][k]; 
					sortedRets[i][k] = tempValue; 
				}
			}
		}

	}


	initialProb = 0.0; 
	terminalProb = 1.0; 
	stepProb = 0.02; 

	cutoff = 0.0; 

	sumProb = 0.0; 
	FDR = 0.0; 

	Lc = 0; 
	n = 0; 


	cutoff = initialProb; 

	while (cutoff <= terminalProb + stepProb / 100)
	{
	    cutoffV.push_back(cutoff); 
	    cutoff += stepProb; 
	}

	rets2.push_back(cutoffV); 

	cutoffV.clear(); 


	for ( int i = 0; i < sortedRets.size(); i++)
	{

		cutoff = initialProb; 
		n = 0; 

		vector<double> FDRV; 

		while (cutoff <= terminalProb + stepProb / 100)
		{

			sumProb = 0.0; 
			FDR = 0.0; 

			Lc = 0; 

			for (int j = sortedRets[i].size() - 1; j >= 0; j--)
			{
				if (sortedRets[i][j] >= cutoff)
				{
					sumProb += sortedRets[i][j]; 
					Lc++; 
				}
				else
				{
					break; 
				}
			}

			if (Lc != 0)
			{
				FDR = 1 - sumProb / Lc; 
			}

			FDRV.push_back(FDR); 

			cutoff += stepProb; 
			n++; 
		}

		rets2.push_back(FDRV); 


		FDRV.clear(); 

		sumProb = 0.0; 
		FDR = 0.0; 

		Lc = 0; 

		FDRV.push_back(FDR); 

		for (int j = sortedRets[i].size() - 1; j >= 0; j--)
		{
			sumProb += sortedRets[i][j]; 
			Lc++; 

			FDR = 1 - sumProb / Lc; 
			FDRV.push_back(FDR); 
		}

		for (int  j = sortedRets[i].size() - 2; j >= 0; j--)
		{
			if (sortedRets[i][j] == sortedRets[i][j + 1])
			{
				FDRV[sortedRets[i].size() - j - 1] = -1.0; 
			}
		}

		rets4.push_back(FDRV); 

	}
    
    return true; 

}



/*----------------------------------------------------------------------*/



/*--------------------------------------------------------------------*/




