#define _CRT_SECURE_NO_WARNINGS

#include <string>
#include <vector>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "gu99.h"
#include "gz97.h"
#include "GZf2.h"
#include "common.h"
#include "tree.h"
#include "matrices.h"
#include "matrix.h"


//----------------------------------------------------------------------

using namespace std;

//----------------------------------------------------------------------
 
bool
rvs_compute(const vector<Tree> &trees, const vector<sequence_t> &sequences,
	     vector<vector<double> > &summary, vector<vector<double> > &rets2) {
  int ntrees = trees.size();

  rets2.resize(ntrees * 2);
  summary.resize(ntrees);
  
  DVector freq(jtt_freq, 20);
  DMatrix2D prob(jtt_prob, 20, 20);

  int i, j, k;
  
  for(i = 0; i < ntrees; i++) {
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
    if(!trees[i].gen_str_wrt_seq(taxa, tree_str)) return false;
    
    IMatrix2D gu_tree;

    if(!parse_tree(tree_str.c_str(), alignment.rows(), treeE, gu_tree)) {
      return false;
    }

    double alpha=2.4, se_alpha=0.0;
    
    double inferred_alpha = 0;
    DVector rets;
    if(!gz97_compute_adapter_for_rvs(true, alpha, alignment,
		     treeE, gu_tree, freq, prob, rets, inferred_alpha)) {
      return false;
    }

    int N = rets.size();
    double a = 0.0;
    double D = 0.0;
    for(k=0; k<N; k++) {
      double ak = 0.0;
      for(j=0; j<(int)(rets(j) + 0.5); j++) {
	ak -= 1.0 / ( (alpha + j) * (alpha + j) );
      }
      a += ak;
      D += rets(k);
    }
    D /= N;
    double b = - N * (D + alpha) / alpha;
    double c = - N * (D + alpha) * (D + alpha) * (1 / (D * D) + 1 / (alpha * alpha));
    
    rets2[i*2+0].resize(N);
    rets2[i*2+1].resize(N);
    for(k=0; k<N; k++) {
      rets2[i*2+0][k] = rets(k);
      rets2[i*2+1][k] = (rets(k) + alpha) / (D + alpha);
    }

    //cout << a << '\t' << b << '\t' << c << '\t' << D << '\t' << N << '\n';
    
    double denom = b*b - a*c;
    if(denom != 0.0) {
      double var_alpha = c / denom;
      if(var_alpha < 0.0) {
	se_alpha = -1.0;
      } else {
	se_alpha = sqrt(var_alpha);
      }
    } else {
      se_alpha = -1.0;
    }
	
    summary[i].resize(3);
    summary[i][0] = inferred_alpha;
  //summary[i][1] = se_alpha;
	summary[i][1] = D;
    summary[i][2] = int(N);
 
  }

#if 0
  {
    cout << "rets.size(): " << rets.size() << '\n';
    vector<vector<double> >::const_iterator i;
    vector<double>::const_iterator j;
    for(i = rets.begin(); i != rets.end(); i++) {
      for(j =i->begin(); j != i->end(); j++) {
	cout << *j << '\t';
      }
      cout << '\n';
    }
  }
#endif
  
  return true;
}

//----------------------------------------------------------------------
