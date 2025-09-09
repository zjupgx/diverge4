#define _CRT_SECURE_NO_WARNINGS
#define WIN32
#include <string>
#include <vector>

#include <stdio.h>
#include <stdlib.h>

#include "gu99.h"
#include "gz97.h"
#include "GZf2.h"
#include "common.h"
#include "tree.h"
#include "matrices.h"
#include "matrix.h"

//#define DB_GP

//----------------------------------------------------------------------

using namespace std;
extern Tree input_tree;
//----------------------------------------------------------------------

bool
gu99_compute(const vector<Tree> &trees, const vector<sequence_t> &sequences,
	     vector<vector<double> > &summary, vector<vector<double> > &rets2) {
  int ntrees = trees.size();
  
  vector<DVector> rets(ntrees);

  DVector freq(jtt_freq, 20);
  DMatrix2D prob(jtt_prob, 20, 20);

  for(int i = 0; i < ntrees; i++) {
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
#ifdef DB_GP
		  printf("i %s\n",(*i).c_str());
#endif
		  for(j = sequences.begin(); j != sequences.end(); j++) {
#ifdef DB_GP
		  //printf("j %s\n",(j->label).c_str());
#endif
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
#if 0
	string mytree_str;
	vector<string> mytaxa;
	input_tree.gen_str_wrt_seq(mytaxa, mytree_str);
	parse_tree(mytree_str.c_str(),alignment.rows()));
#endif
	//printf("trees[i].gen_str():%s\ntree[i].gen_str_wrt_seq():%s\n",(trees[i].gen_str()).c_str(),tree_str.c_str());


    IMatrix2D gu_tree;

    if(!parse_tree(tree_str.c_str(), alignment.rows(), treeE, gu_tree)) {
		return false;
    }
#ifdef DB_GP
	for(int mm=0;mm<2*(alignment.rows());++mm)
	{
		printf("gu_tree(%d,0)\ngu_tree(%i,1)\n", gu_tree(mm,0), gu_tree(mm,1));
	}
#endif
    if(!gz97_compute(true, 2.4, alignment, treeE, gu_tree, freq, prob, rets[i])) {
		return false;
    }
  }

  {
    vector<DVector> &numSub = rets;
    
    if(!GZf2_compute(numSub, summary, rets2)) return false;
  }

  return true;
}

//----------------------------------------------------------------------
