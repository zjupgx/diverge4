#include <string>
#include <vector>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "type_two.h"
#include "gz97.h"
#include "gu99.h"
#include "GZf2.h"
#include "common.h"
#include "tree.h"
#include "matrices.h"
#include "matrix.h"

#include "progress.h"
#include "ancestral_compute.h"

#include <iostream>
#include <fstream>

//----------------------------------------------------------------------

using namespace std;

//----------------------------------------------------------------------
enum SITE_CLASSIFY {IDENTICAL, SAMEGROUP, DIFFGROUP};

//----------------------------------------------------------------------
bool type_two_compute(const vector<Tree> &trees, const vector<sequence_t> &sequences,
	     vector<vector<double> > &summary, vector<vector<double> > &rets2);
bool
type_two_divergence(const Tree &t1, const Tree &t2, const vector<sequence_t> &s1,
					const vector<sequence_t> &s2, const double alpha, const vector<double> &tree1_numsub,
					const vector<double> &tree2_numsub, vector<double> &summary, 
					vector<double> &rets);
bool GZ97_compute_adapter(const vector<Tree> &trees, const vector<sequence_t> &sequences,
					 vector<double> &alpha, vector< vector<double> > &numSub);
int amino_acid_group(char a);
enum SITE_CLASSIFY amino_acid_compare(char a, char b);
int num_site_with_amino_acid_classify(enum SITE_CLASSIFY cl, const sequence_t & seq1,
							 const sequence_t &seq2);
void site_profile();
bool ancestral_inference(const Tree &t, const vector<sequence_t> &seq, 
						 sequence_t &ans_seq);
int F(const vector<double> &changes1, const vector<double> &changes2,
			 const sequence_t &root1, const sequence_t &root2, enum SITE_CLASSIFY c);
int G(sequence_t &s1, sequence_t &s2, enum SITE_CLASSIFY c);

void construct_seq(const Tree &t, const vector<sequence_t> &seqs, vector<sequence_t> &s);

//----------------------------------------------------------------------

bool
type_two_compute(const vector<Tree> &trees, const vector<sequence_t> &sequences,
	     vector<vector<double> > &summary, vector<vector<double> > &rets2) 
{
  int ntrees = trees.size();
  int i, j, totalcount;

  //find out the alpha_ML and number of substition 
  vector<double> alpha;
  vector< vector<double> > numsubs;
  GZ97_compute_adapter(trees, sequences, alpha, numsubs);
  ///////////////////

  totalcount = 0;

  for(i = 0; i < ntrees; i++) {
	  for(j = i + 1; j < ntrees; j++, totalcount++) {
		vector<double> sum;
		vector<double> ret;
		vector<sequence_t> s1, s2;
		Tree t1 = trees[i];
		Tree t2 = trees[j];
		construct_seq(t1, sequences, s1);
		construct_seq(t2, sequences, s2);

		type_two_divergence(t1, t2, s1, s2, alpha[totalcount], numsubs[i], numsubs[j] , sum, ret);
		summary.push_back(sum);
		rets2.push_back(ret);
		Progress::increment(100);
	  }
  }

  return true;
}

bool
type_two_divergence(const Tree &t1, const Tree &t2, const vector<sequence_t> &s1,
					const vector<sequence_t> &s2, const double alpha, const vector<double> &tree1_numsub,
					const vector<double> &tree2_numsub, vector<double> &summary, 
					vector<double> &rets) {
	
	int site_num;
	int tree1_seqnum, tree2_seqnum;
	int i;
		
	site_num = s1.front().sequence.size();
	tree1_seqnum = s1.size();
	tree2_seqnum = s2.size();

	//---
	double Da, Db;


	Da = 0;
	Db = 0;
	for(i = 0; i < site_num; i++)
		Da += tree1_numsub[i];
	Da = Da/site_num;

	for(i = 0; i < site_num; i++)
		Db += tree2_numsub[i];
	Db = Db/site_num;
	//---
	sequence_t a, b;
	int N, C, R;

	ancestral_inference(t1, s1, a);
	ancestral_inference(t2, s2, b);

	N = num_site_with_amino_acid_classify(IDENTICAL, a, b);
	C = num_site_with_amino_acid_classify(SAMEGROUP, a, b);
	R = num_site_with_amino_acid_classify(DIFFGROUP, a, b);



	double P = 1 - ((double)N)/site_num;
	double F00_N = ((double)F(tree1_numsub, tree2_numsub, a, b, IDENTICAL))/site_num;
	double F00_C = ((double)F(tree1_numsub, tree2_numsub, a, b, SAMEGROUP))/site_num;
	double F00_R = ((double)F(tree1_numsub, tree2_numsub, a, b, DIFFGROUP))/site_num;

	double d = 0;
	double W = 0; 
	double theta = 0.0; 
	double tmp = 0.0;

	while(1) {
		d = -1 * log(1 - (P-theta)/(1-theta));
		W = pow(alpha/(Da+Db+d+alpha), alpha);
		theta = 1 - F00_N/W;
		if(theta <= 0.0)
			break;
		else {
			if(fabs(tmp - theta) < 0.001)
				break;
			else
				tmp = theta;
		}
	}

	double Z = pow(alpha/(Da+Db+alpha), alpha);

	double theta_se = sqrt( ((double)F00_N * (1-F00_N))/(N+C+R) )/W;

	double Gr = ((double)R)/(R+C);
	double Gc = ((double)C)/(R+C);

	double PIr=0.31213;
	double Ar = (P*Gr - (1-theta) * PIr * (P-theta))/theta;
	//double Ar = (F00_R - (Z-W) * Gr)/(theta * W);
	//double PIr = ( Gr * Z - F00_R )/((1 - theta) * W);
    

	double Ac = 1 - Ar;
	double PIc = 1 - PIr;

	double h = d/(Da + Db + d + alpha);
	double Q = ((double)R)/(1 + R);

	for(i = 0; i < site_num; i++) {
		double R;
		char tmp1 = a.sequence[i];
		char tmp2 = b.sequence[i];
		switch(amino_acid_compare(tmp1, tmp2)) {
		case IDENTICAL:
			R = 0.0;
			break;
		case SAMEGROUP:
			R = (theta * Ac)/((1-theta) * PIc * (1-pow((1-h), (tree1_numsub[i] + tree2_numsub[i] + alpha))));
			break;
		case DIFFGROUP:
			R = (theta * Ar)/((1-theta) * PIr * (1-pow((1-h), (tree1_numsub[i] + tree2_numsub[i] + alpha))));
			break;
		}
		rets.push_back(R);
	}

	summary.push_back(Da);
	summary.push_back(Db);
	summary.push_back(N);
	summary.push_back(C);
	summary.push_back(R);
	summary.push_back(alpha);
	summary.push_back(theta);
	summary.push_back(theta_se);
	summary.push_back(Ar);
	summary.push_back(PIr);
	summary.push_back(P);
	summary.push_back(d);
	summary.push_back(W);
	summary.push_back(Z);
	summary.push_back(Gr);
	summary.push_back(Gc);
	summary.push_back(h);
	summary.push_back(Q);
	summary.push_back(F00_N);
	summary.push_back(F00_R);
	summary.push_back(F00_C);



	return true;
}

//----------------------------------------------------------------------
/*This method is to calculate the number of site changes
and alpha value. It is adapted from Gu99 file.
 */


bool
GZ97_compute_adapter(const vector<Tree> &trees, const vector<sequence_t> &sequences,
					 vector<double> &alpha, vector< vector<double> > &numSub) {
    

  vector< vector<double> > rets2;
  vector< vector<double> > summary;
  int i, j;

  int ntrees = trees.size();
  
  vector<DVector> rets(ntrees);

  DVector freq(jtt_freq, 20);
  DMatrix2D prob(jtt_prob, 20, 20);

  for(i = 0; i < ntrees; i++) {
    int treeE;

    if(Progress::wasCancelled()) return false;
    
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

    if(!gz97_compute(true, 2.4, alignment,
		     treeE, gu_tree, freq, prob, rets[i])) {
      return false;
    }
  }
   
  vector<DVector> &ns = rets;

  if(!GZf2_compute(ns, summary, rets2)) 
	  return false;

  //post processed
  for(i = 0; i < ntrees; i++) {
	  vector<double> tmp;
	  for(j = 0; j < rets[0].size(); j++) {
		  tmp.push_back(rets[i](j));
	  }
	  numSub.push_back(tmp);
  }


  for(i = 0; i < summary.size(); i++) {
	alpha.push_back(summary[i][6]);
  }
  


  return true;
}

//----------------------------------------------------------------------

int amino_acid_group(char a) {
	char upper_a = toupper(a);
	switch(a) {
	case 'K':
	case 'R': 
	case 'H': 
		return 1;
		break;
	case 'D':
	case 'E':
		return 2;
		break;
	case 'S':
	case 'T': 
	case 'N': 
	case 'Q':
	case 'C': 
	case 'G': 
	case 'P': 
		return 3;
		break;
	case 'A':
	case 'I': 
	case 'L': 
	case 'M': 
	case 'F': 
	case 'W': 
	case 'V': 
	case 'Y': 
		return 4;
		break;
	default:
		return 0;
	}
}

//----------------------------------------------------------------------
enum SITE_CLASSIFY amino_acid_compare(char a, char b) {
	if(a == b)
		return IDENTICAL;
	else if(amino_acid_group(a) == amino_acid_group(b))
		return SAMEGROUP;
	else if(amino_acid_group(a) != amino_acid_group(b))
		return DIFFGROUP;
}
//----------------------------------------------------------------------

int num_site_with_amino_acid_classify(enum SITE_CLASSIFY cl, const sequence_t & seq1,
							 const sequence_t &seq2) {
	int result = 0, i;
	char c1, c2;
	if(seq1.sequence.size() != seq2.sequence.size())
		return 0;
	for(i = 0; i < seq1.sequence.length(); i++) {
		c1 = seq1.sequence.at(i);
		c2 = seq2.sequence.at(i);
		if(c1 == c2 && cl == IDENTICAL)
			result++;
		else if(c1 != c2 && amino_acid_group(c1) == amino_acid_group(c2) && 
			cl == SAMEGROUP)
			result++;
		else if(amino_acid_group(c1) != amino_acid_group(c2) && 
			(cl == DIFFGROUP))
			result++;
	}
	return result;
}
//----------------------------------------------------------------------

void site_profile() {

}

//----------------------------------------------------------------------

bool ancestral_inference(const Tree &t, const vector<sequence_t> &seq, 
						 sequence_t &ans_seq) {
	vector<sequence_t> seq_out;

	Tree tt = t;
	tt.reorder_id();

	QString warn;
	ancestral_compute(tt, seq, seq_out, warn);


	ans_seq = seq_out[seq.size()+1];

	
	
/*
	ofstream out("c:\\temp\\test.data");
  
	if(!out) {
		return false;
	}

	  string tree_str;
  {
    vector<string> str_seqNames;
    
    vector<sequence_t>::const_iterator i;
    for(i = seq.begin(); i != seq.end(); i++) {
      str_seqNames.push_back(i->label);
    }
    
    if(!t.gen_str_wrt_seq(str_seqNames, tree_str)) {
      return false;
    }
  }

	out << seq_out.size() << "   " << seq_out[0].sequence.size() << endl;

  for(int ii = 0; ii < seq_out.size(); ii++) {
		string s = seq_out[ii].sequence;
		out << seq_out[ii].label << endl;
		out << s << endl;
  }

  out << tree_str << endl;
  out.close();
	
	  */
	return true;
 }

//----------------------------------------------------------------------

//----------------------------------------------------------------------

int F(const vector<double> &changes1, const vector<double> &changes2,
			 const sequence_t &root1, const sequence_t &root2, enum SITE_CLASSIFY c) {
	int result = 0;
	int i, j;

	int num_site = root1.sequence.size();

	for(i = 0; i < num_site; i++) {
		if(changes1[i] == 0 && changes2[i] == 0 && 
			amino_acid_compare(root1.sequence.at(i), root2.sequence.at(i)) == c)
			result++;
	}
	return result;
}

int G(sequence_t &s1, sequence_t &s2, enum SITE_CLASSIFY c) {
	int num_site = s1.sequence.size();
	int i, result;

	for(i = 0; i < num_site; i++)
		if (amino_acid_compare(s1.sequence.at(i), s2.sequence.at(i)) == c)
			result++;

	return result;
}

void construct_seq(const Tree &t, const vector<sequence_t> &seqs, vector<sequence_t> &s) {
	vector<string> names;
	int i, j;

	t.leaf_names(names);
	string t1, t2;
	for(i = 0; i < names.size(); i++) {
		for(j = 0; j < seqs.size(); j++) {
			t1 = seqs[j].label;
			t2 = names[i];
			if(names[i] == seqs[j].label)
				s.push_back(seqs[j]);
		}
	}
}


