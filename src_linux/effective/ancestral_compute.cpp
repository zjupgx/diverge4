#define _CRT_SECURE_NO_WARNINGS

#include <vector>
#include <string>
#include <iostream>
#include <fstream>


#include <stdio.h>
#include <stdlib.h>


#include "gu99.h"
#include "gz97.h"
#include "common.h"
#include "tree.h"
#include "matrices.h"
#include "sequence.h"


//----------------------------------------------------------------------

using namespace std;

//----------------------------------------------------------------------

static bool
/*ClusteringTab::*/containsPolytomy(const node_t *node, bool root) {
  if(!node) return false;
  if(node->n_children() > 3 || !root && node->n_children() > 2) return true;
  const list<node_t*> &children = node->children();
  list<node_t*>::const_iterator i;
  for(i=children.begin(); i!=children.end(); i++) {
    if(containsPolytomy(*i, false)) return true;
  }
  return false;
}

bool
ancestral_compute(Tree &in_tree, const vector<sequence_t> &orig_sequences,
				  vector<sequence_t> &sequencesOut, std::string & warning) {
  warning = "";

  if(in_tree.breadth() <= 3)  {
	  warning = "The number of tree leaf is less than 4";
	  return false;
  }

  if(containsPolytomy(in_tree.root(), true)) {
    warning = "Selected tree contains an unsupported polytomy.";
    return false;
  }

  Tree tree;
  if(!in_tree.polyroot(tree)) 
	  abort();

  tree.assign_ids();

  std::vector<sequence_t> sequences;
  std::vector<std::string> leaf_names;
  in_tree.leaf_names(leaf_names);
  if(!reorder_sequences(leaf_names, orig_sequences, sequences)) {
	  warning = "Reorder error: please contact the authors for detailed information!";
	  return false;
  }
  

  DVector freq(jtt_freq, 20);
  DMatrix2D prob(jtt_prob, 20, 20);
  
  int treeE, k;

  CMatrix2D alignment(sequences.size(), sequences[0].sequence.size());
  
  {
    // generate sub_sequences and sub_seqNames
    vector<sequence_t>::const_iterator j;
    int i;
    for(i=0, j = sequences.begin(); j != sequences.end(); j++, i++) {
      for(k=0; k<(int)j->sequence.size(); k++) {
		  alignment(i, k) = j->sequence[k];
      }
    }
  }

  string tree_str;
  vector<string> str_seqNames;
  { 
    vector<sequence_t>::const_iterator i;
    for(i = sequences.begin(); i != sequences.end(); i++) {
      str_seqNames.push_back(i->label);
    }
    
    if(!tree.gen_str_wrt_seq(str_seqNames, tree_str)) {
		warning = "Tree sequence error: please contact the authors for more information!";
      return false;
    }
  }

  //cout << tree_str << '\n';
  ///////////////TEST ANCESTOR/////////////
  /*ofstream out("c:\\temp\\test.data");
  
  if(!out) {
		return false;
  }

  out << sequences.size() << "   " << sequences[0].sequence.size() << endl;

  for(int ii = 0; ii < sequences.size(); ii++) {
		string s = sequences[ii].sequence;
		out << sequences[ii].label << endl;
		out << s << endl;
  }

  out << tree_str << endl;
  out.close();
*/
  

  /////////////////////////////////////////
  IMatrix2D gu_tree;

  if(!parse_tree(tree_str.c_str(), alignment.rows(),
		 treeE, gu_tree)) {
	  warning = "Parse tree error: please contact the authors for more information!";
    return false;
  }

  /*ofstream out("c:\\temp\\test1.txt");
  
  if(!out) {
		return false;
  }

  out << tree_str << endl;

  for(int ii = 0; ii < gu_tree.rows(); ii++) {
	  for( int jj = 0; jj < gu_tree.cols(); jj++)
		  out << gu_tree(ii, jj) << ", ";
	  out << endl;
  }

  for(int iii = 0; iii < str_seqNames.size(); iii++) {
	out << str_seqNames[iii] << endl;
  }

  out.close();
  */


  CMatrix2D alignmentOut;
  if(!ancestral_inference(true, 2.4,
			  alignment,
			  treeE, gu_tree, freq, prob,
			  alignmentOut)) {
	  warning = "Ancestral inference failure!";
	  return false;
  }

  //the main reason for to do reorder is: the tree shown on the window is different
  //from the internal one (tree) used for infering the ancestral genes. Therefore,
  //we have to remap the internal node of tree to the internal nodes of in_tree. 
  //we need to do the reordering using the following procedure:
  //- if one of the child is a leaf, then ask tree to provide its parent id
  //- if both of them are internal nodes, then add their id by 1 and ask tree to provide its
  //  parent id
  vector<int> internal_node_id;
  int i, j;
  internal_node_id.clear();
  int parent;
  
  for(i = sequences.size() + 1; i < 2 * sequences.size() - 1; i++) {
	  if(gu_tree(i, 0) <= sequences.size()) {
		  int id = in_tree.id_lookup(str_seqNames[gu_tree(i, 0) - 1]);
		  if(!in_tree.parent_id(id, parent) || parent == -1) {
			  warning = "Node mapping error: please contact the authors for more information!";
			  return false;
		  }
		  internal_node_id.push_back(parent);
	  }
	  else if (gu_tree(i, 1) <= sequences.size()) {
		  int id = in_tree.id_lookup(str_seqNames[gu_tree(i, 1) - 1]);
		  if(!in_tree.parent_id(id, parent) || parent == -1) {
			  warning = "Node mapping error: please contact the authors for more information!";
			  return false;
		  }
		  internal_node_id.push_back(parent);
	  }
	  else {
		  //first, we collect the set of child nodes and their relative depth to current node
		  vector<int> tmp1, tmp2;
		  vector<int> children, depth;

		  tmp1.push_back(gu_tree(i,0)); tmp2.push_back(1);
		  tmp1.push_back(gu_tree(i,1)); tmp2.push_back(1);

		  while(!tmp1.empty()) {
			  int b = tmp1.back();
			  int d = tmp2.back();
			  tmp1.pop_back();
			  tmp2.pop_back();
			  if(b <= sequences.size()) {
				  children.push_back(b);
				  depth.push_back(d);
			  }
			  else {
				  	tmp1.push_back(gu_tree(b,0)); tmp2.push_back(d+1);
					tmp1.push_back(gu_tree(b,1)); tmp2.push_back(d+1);
			  }
		  }

		  //second, we use the child node string to find out its id on original
		  //tree and the depth in that tree
		  vector<int> orig_children, orig_depth;
		  orig_children.clear();
		  orig_depth.clear();

		  int j;

		  for(j = 0; j < children.size(); j++) {
			  orig_children.push_back(in_tree.id_lookup(str_seqNames[children[j]-1]));
			  orig_depth.push_back(in_tree.id_depth(orig_children[j]));
		  }

		  //third, find the node whose depth in orginal tree is greater than on 
		  //new tree and trace back
		  for(j = 0; j < children.size(); j++) {
			  if(orig_depth[j] > depth[j]) {
				  vector<int> parents;
				  if(!in_tree.parents(orig_children[j], parents)) {
					  warning = "Node mapping error: please contact the authors for more information!";
					  return false;
				  }
				  internal_node_id.push_back(parents[orig_depth[j] - depth[j]]);
				  break;
			  }
		  }
	  }
  }


  sequencesOut.resize(alignmentOut.rows());
  
  for(i=0; i<alignmentOut.rows()-1; i++) {
    if(i < sequences.size()) {
      sequencesOut[i] = sequences[i];
    } else if(i == sequences.size()) {
      // skip
    } else {
      char str[128];
      sprintf(str, "Internal Node %d", internal_node_id[i - sequences.size() - 1]);
      string tmp = "";
      for(j=0; j<alignmentOut.cols(); j++) {
		  tmp += alignmentOut(i, j);
      }
      sequencesOut[i].label = str;
      sequencesOut[i].sequence = tmp;
    }
  }

  /*ofstream out("c:\\temp\\test.data");
  
  if(!out) {
		return false;
  }

  out << sequencesOut.size() << "   " << sequencesOut[0].sequence.size() << endl;

  for(int ii = 0; ii < sequencesOut.size(); ii++) {
		string s = sequencesOut[ii].sequence;
		out << sequencesOut[ii].label << endl;
		out << s << endl;
  }

  out << tree_str << endl;
  out.close();*/

  return true;
}

//----------------------------------------------------------------------
