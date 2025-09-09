#define _CRT_SECURE_NO_WARNINGS

#include <string>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <iostream>
#include "tree.h"

//======================================================================

using namespace std;

//======================================================================

static char *
str_clean_rn(char *str) {
  char *p;

  if((p=strrchr(str, '\r'))) *p='\0';
  if((p=strrchr(str, '\n'))) *p='\0';

  return str;
}

//----------------------------------------------------------------------

static bool
get_entire_line(FILE *file, string &line) {
  line = "";
  
  int c;
  do {
    c = fgetc(file);
    if(c != EOF) {
      line += (char)c;
    }
  } while(c != EOF);

  return feof(file) != 0;
}

//----------------------------------------------------------------------

static int
getc_nws(string &str) {
  if(str.empty()) return EOF;
  int c;
  do {
    c = str[0];
    str.erase(0, 1);
  } while(isspace(c));
  return c;
}

//======================================================================

node_t::node_t()
  : _length(0.0),
    _x(0.0), _y(0.0),
    _min_x(1e6),_max_x(-1e6),
    _min_y(1e6), _max_y(-1e6),
    _selected(false),
    _highlight(false),
    _id(-1)
{}

//----------------------------------------------------------------------

node_t::node_t(const node_t &n)
  : _label(n._label),
    _length(n._length),
    _x(n._x), _y(n._y),
    _min_x(n._min_x),_max_x(n._max_x),
    _min_y(n._min_y), _max_y(n._max_y),
    _selected(n._selected),
    _highlight(n._highlight),
    _id(n._id)
{
  list<node_t *>::const_iterator i;
  for(i = n._children.begin(); i != n._children.end(); i++) {
    add_node(new node_t(**i));
  }
}

void node_t::reorder_id(int &id) {
	list<node_t *>::const_iterator j;

	if(n_children() > 0) {
		for(j = _children.begin(); j != _children.end(); j++) {
			(*j)->reorder_id(id); 
		}
		_id = id;
		id++;
	}
	
}

void node_t::generate_profile(int * profile, double * branch, int total_leaves) {
	list<node_t *>::const_iterator j;
	node_t * middle;
	int i = 0;
	for(j = _children.begin(); j != _children.end(); j++) {
		(*j)->generate_profile(profile, branch, total_leaves); 
		if(i == 1) 
		  middle = *j;
		i++;
	}


  //try

  try{

    if(n_children() == 0) {// no children, it is a leaf
	  profile[3 * (_id - 1) + 0] = 0;
	  profile[3 * (_id - 1) + 2] = 0;
	  profile[3 * (_id - 1) + 1] = _id;
	  branch[_id - 1] = _length;
	}
	else if(n_children() == 2) {//two children, normal internal node
	  profile[3 * (_id - 1) + 0] = (_children.front())->id();
	  profile[3 * (_id - 1) + 2] = (_children.back())->id();
	  profile[3 * (_id - 1) + 1] = _id;
	  branch[_id - 1] = _length;
	}
	else if(n_children() == 3) {//root node
	  //
	  profile[3 * (_id - 1) + 0] = (_children.front())->id();
	  profile[3 * (_id - 1) + 2] = (_children.back())->id();
	  profile[3 * (_id - 1) + 1] = _id;
	  branch[_id - 1] = (middle->length() / 2);
	  
	  
	  int total_nodes = 2 * total_leaves - 1;

	  //create a new node as root
	  profile[3 * (total_nodes - 1) + 0] = _id;
	  profile[3 * (total_nodes - 1) + 2] = middle->id();
	  profile[3 * (total_nodes - 1) + 1] = total_nodes;
	  branch[(total_nodes - 1)] =  (middle->length() / 2);
	  
	}
	else {
	  throw std::invalid_argument("error in treefile.");
	}
	}
	catch(const char* &e){
		cout<<e<<endl;
	}

  //end 


	return;
}

//----------------------------------------------------------------------

node_t::node_t(const string &label, double length)
  : _label(label),
    _length(length),
    _x(0.0), _y(0.0),
    _min_x(1e6),_max_x(-1e6),
    _min_y(1e6), _max_y(-1e6),
    _selected(false),
    _highlight(false),
    _id(-1)
{}

//----------------------------------------------------------------------
const node_t * Tree::getSelectedRoot(node_t * node) const {
  if(!node) return NULL;

  if(node->_selected)
	  return node;
  
  list<node_t *>::const_iterator i;
  for(i = node->_children.begin(); i != node->_children.end(); i++) {
    const node_t * k = getSelectedRoot(*i);
	if(k)
		return k;
  }

  return NULL;
}

//----------------------------------------------------------------------
const node_t * Tree::getSelectedRoot() const {
  return getSelectedRoot(root_);
}
//======================================================================
bool node_t::flipNode() {
	node_t * b, *f;
	if(_children.size() <= 0)
		return false;
	f = _children.front();
	b = _children.back();
	_children.clear();
	_children.push_back(b);
	_children.push_back(f);

	return true;
}
//----------------------------------------------------------------------

node_t::~node_t() {
  list<node_t *>::iterator i;
  for(i = _children.begin(); i != _children.end(); i++) {
    delete *i;
  }
}




int

node_t::id_depth(int id) {

	if(_id == id)

		return 0;

	else {

		list<node_t *>::iterator i;

		for(i = _children.begin(); i!= _children.end(); i++) {

			if((*i)->contain_id(id)) {

				return ((*i)->id_depth(id) + 1);

			}

		}

	}

	return -1;

}
//----------------------------------------------------------------------

node_t &
node_t::operator = (const node_t &n) {
  _label = n._label;
  _length = n._length;
  _x = n._x;
  _y = n._y;
  _min_x = n._min_x;
  _max_x = n._max_x;
  _min_y = n._min_y;
  _max_y = n._max_y;
  _selected = n._selected;
  _highlight = n._highlight;

  list<node_t *>::const_iterator i;
  for(i = n._children.begin(); i != n._children.end(); i++) {
    delete *i;
  }

  _children.clear();
  
  for(i = n._children.begin(); i != n._children.end(); i++) {
    add_node(new node_t(**i));
  }

  return *this;
}

//----------------------------------------------------------------------

void
node_t::add_node(node_t *n) {
  _children.push_back(n);
}



bool 

node_t::parents(int id, std::vector<int> &parents) {

	if(_id == id) {

		return true;

	}

	else {

		

		list<node_t *>::const_iterator i;

		for(i = _children.begin(); i != _children.end(); i++) {

			if((*i)->contain_id(id)) {

				parents.push_back(_id);

				(*i)->parents(id, parents);

				return true;

			}

		}



		return false;

	}

}



//----------------------------------------------------------------------


bool 

node_t::contain_id(int id) {

	if(_id == id)

		return true;

	

	if(n_children() <= 0)

		return false;



	list<node_t *>::const_iterator i;

	for(i = _children.begin(); i != _children.end(); i++) {

		if((*i)->contain_id(id))

			return true;

	}



	return false;

}


//----------------------------------------------------------------------

bool
node_t::remove_node(node_t *n) {
  list<node_t *>::iterator i = find(_children.begin(), _children.end(), n);
  if(i != _children.end()) {
    _children.erase(i);
    return true;
  }
  return false;
}

//----------------------------------------------------------------------

int
node_t::height() const {
  int ret = 1;

  list<node_t *>::const_iterator i;
  for(i = _children.begin(); i != _children.end(); i++) {
    ret = max(ret, (*i)->height()+1);
  }

  return ret;
}



//----------------------------------------------------------------------

//if not found, return false;

//if it is the root, return true with parent_id == -1;

//otherwise, return ture with the real parent id;

bool

node_t::parent_id(int child_id, int & parent_id) {

	if(this->_id == child_id) {

		parent_id = -1;

		return true;

	}

	else {

		list<node_t *>::const_iterator i;

		for(i = _children.begin(); i != _children.end(); i++) {

			int id;

			if((*i)->parent_id(child_id, id)) {

				if(id == -1)

					parent_id = _id;

				else 

					parent_id = id;

				return true;

			}

		}

	}



	return false;

}

//----------------------------------------------------------------------

double
node_t::max_branch_length() const {
  double ret = 0.0;

  list<node_t *>::const_iterator i;
  for(i = _children.begin(); i != _children.end(); i++) {
    ret = max(ret, (*i)->max_branch_length());
  }

  return ret + _length;
}

//----------------------------------------------------------------------

double
node_t::total_branch_length() const {
  double ret = _length;

  list<node_t *>::const_iterator i;
  for(i = _children.begin(); i != _children.end(); i++) {
    ret += (*i)->total_branch_length();
  }

  return ret;
}

//----------------------------------------------------------------------

int
node_t::n_nodes() const {
  int ret = 1;

  list<node_t *>::const_iterator i;
  for(i = _children.begin(); i != _children.end(); i++) {
    ret += (*i)->n_nodes();
  }

  return ret;
}

//----------------------------------------------------------------------

int
node_t::n_children() const {
  return _children.size();
}

//----------------------------------------------------------------------

int
node_t::breadth() const {
  int ret = 0;

  if(_children.empty()) {
    ret = 1;
  } else {
    list<node_t *>::const_iterator i;
    for(i = _children.begin(); i != _children.end(); i++) {
      ret += (*i)->breadth();
    }
  }

  return ret;
}

//----------------------------------------------------------------------

int
node_t::id() const {
  return _id;
}

//----------------------------------------------------------------------

#if 0
void
node_t::print(int indent) const {
  for(int i=0; i<indent; i++) {
    cout << ' ';
  }

  cout << '(' << _id << ") ";
  
  if(_label.empty()) {
    cout << '*';
  } else {
    cout << _label;
  }
  
  cout << " @ (" << _x << ", " << _y << ")";
  if(_children.size() > 0) {
    cout << " with " << _children.size() << " children\n";
  } else {
    cout << '\n';
  }
  
  list<node_t *>::const_iterator j;
  for(j = _children.begin(); j != _children.end(); j++) {
    (*j)->print(indent + 2);
  }
}
#endif

//----------------------------------------------------------------------

void
node_t::leaf_names(vector<string> &names) const {
  if(_children.empty()) {
    names.push_back(_label);
  } else {
    list<node_t *>::const_iterator i;
    for(i = _children.begin(); i != _children.end(); i++) {
      (*i)->leaf_names(names);
    }
  }
}

//----------------------------------------------------------------------

void
node_t::leaf_names(set<string> &names) const {
  if(_children.empty()) {
    names.insert(_label);
  } else {
    list<node_t *>::const_iterator i;
    for(i = _children.begin(); i != _children.end(); i++) {
      (*i)->leaf_names(names);
    }
  }
}

//======================================================================

void
Tree::leaf_names(vector<string> &names) const {
  names.clear();
  if(root_) root_->leaf_names(names);
}

//----------------------------------------------------------------------

void
Tree::leaf_names(set<string> &names) const {
  names.clear();
  if(root_) root_->leaf_names(names);
}

//----------------------------------------------------------------------

void
Tree::gen_str(string &str, bool include_lengths, const node_t *tree) const {
  if(!tree) return;
  
  str += '(';
  list<node_t *>::const_iterator i = tree->_children.begin();
  while(i != tree->_children.end()) {
    if(!(*i)->_children.empty()) {
      gen_str(str, include_lengths, *i);
    } else {
      string::const_iterator j;
      for(j = (*i)->_label.begin(); j != (*i)->_label.end(); j++) {
	if(isspace(*j)) {
	  str += '_';
	} else {
	  str += *j;
	}
      }
      if(include_lengths && has_lengths /*(*i)->_length > 0*/) {
	char s[64];
	sprintf(s, ":%f", (*i)->_length);
	str += s;
      }
    }
    i++;
    if(i != tree->_children.end()) {
      str += ',';
    }
  }
  str += ')';
  if(include_lengths && has_lengths /*tree->_length > 0*/) {
    char s[64];
    sprintf(s, ":%f", tree->_length);
    str += s;
  }
}

//----------------------------------------------------------------------

string
Tree::gen_str(bool include_lengths) const {
  string str;
  gen_str(str, include_lengths, root_);
  return str;
}

//----------------------------------------------------------------------

bool
Tree::gen_str_pruned(string &str, bool include_lengths, const node_t *node, const node_t *prune_node) const {
  if(!node) return false;

  bool rv = false;
  string tmp;
  int n = 0;
  
  list<node_t *>::const_iterator i = node->_children.begin();
  while(i != node->_children.end()) {
    if(*i == prune_node) {
      rv = true;
      i++;
    } else {
      if(!(*i)->_children.empty()) {
	rv = gen_str_pruned(tmp, include_lengths, *i, prune_node) || rv;
      } else {
	string::const_iterator j;
	for(j = (*i)->_label.begin(); j != (*i)->_label.end(); j++) {
	  if(isspace(*j)) {
	    tmp += '_';
	  } else {
	    tmp += *j;
	  }
	}
	if(include_lengths && has_lengths /*(*i)->_length > 0*/) {
	  char s[64];
	  sprintf(s, ":%f", (*i)->_length);
	  tmp += s;
	}
      }
      i++;
      if(i != node->_children.end() && *i != prune_node) {
	tmp += ',';
      }
      n++;
    }
  }

  if(n >= 2) {
    str += '(';
    str += tmp;
    str += ')';
    
    if(include_lengths && has_lengths /*tree->_length > 0*/) {
      char s[64];
      sprintf(s, ":%f", node->_length);
      str += s;
    }
  } else {
    str += tmp;
  }
  
  return rv;
}

//----------------------------------------------------------------------

bool
Tree::gen_str_pruned(const node_t *node, string &tree, bool include_lengths) const {
  tree = "";
  return gen_str_pruned(tree, include_lengths, root_, node);
}


//----------------------------------------------------------------------

bool
Tree::gen_str_wrt_seq(string &str, const vector<string> &names, const node_t *tree) const {
  if(!tree) return false;

  

  str += '(';

  list<node_t *>::const_iterator i = tree->_children.begin();

  while(i != tree->_children.end()) {

    if(!(*i)->_children.empty()) {

      if(!gen_str_wrt_seq(str, names, *i)) {

		  return false;

      }

    } else {

      int j;

      for(j=0; j<(int)names.size(); j++) {

		  if((*i)->_label == names[j]) {

			  char s[64];

			  sprintf(s, "%d", j+1);

			  str += s;

			  break;

		  }

      }

      if(j == (int)names.size()) {

        //cerr << "Unable to find " << (*i)->_label << " in sequences\n";

		  return false;

      }

    }

    i++;

    if(i != tree->_children.end()) {

      str += ',';

    }

  }

  str += ')';



  return true; 

  



  /*if(!tree) return false;

  



  if(tree->_children.size() == 3) { //root, so we reverse it

	  str += '(';

	  list<node_t *>::const_reverse_iterator i = tree->_children.rbegin();

	  i++;

	  while(i != tree->_children.rend()) {

		  if(!(*i)->_children.empty()) {

			  if(!gen_str_wrt_seq(str, names, *i)) {

				  return false;

			  }

		  } 

		  else {

			  int j;

			  for(j=0; j<(int)names.size(); j++) {

				  if((*i)->_label == names[j]) {

					  char s[64];

					  sprintf(s, "%d", j+1);

					  str += s;

					  break;

				  }

			  }

			  

			  if(j == (int)names.size()) {

				  //cerr << "Unable to find " << (*i)->_label << " in sequences\n";

				  return false;

			  }

		  }

		  i++;

		  if(i != tree->_children.rend()) {

			  str += ',';

		  }

	  }



	  str += ',';



	  i = tree->_children.rbegin();

	  if(!(*i)->_children.empty()) {

		  if(!gen_str_wrt_seq(str, names, *i)) {

			  return false;

		  }

	  } 

	  else {

		  int j;

		  for(j=0; j<(int)names.size(); j++) {

			  if((*i)->_label == names[j]) {

				  char s[64];

				  sprintf(s, "%d", j+1);

				  str += s;

				  break;

			  }

		  }

			  

		  if(j == (int)names.size()) {

			  //cerr << "Unable to find " << (*i)->_label << " in sequences\n";

			  return false;

		  }

	  }



	  str += ')';

  }

  else {

	  str += '(';

	  list<node_t *>::const_iterator i = tree->_children.begin();

	  while(i != tree->_children.end()) {

		  if(!(*i)->_children.empty()) {

			  if(!gen_str_wrt_seq(str, names, *i)) {

				  return false;

			  }

		  } 

		  else {

			  int j;

			  for(j=0; j<(int)names.size(); j++) {

				  if((*i)->_label == names[j]) {

					  char s[64];

					  sprintf(s, "%d", j+1);

					  str += s;

					  break;

				  }

			  }

			  

			  if(j == (int)names.size()) {

				  //cerr << "Unable to find " << (*i)->_label << " in sequences\n";

				  return false;

			  }

		  }

		  i++;

		  if(i != tree->_children.end()) {

			  str += ',';

		  }

	  }

	  str += ')';

  }

  

  return true;

  */
}

//----------------------------------------------------------------------

bool
Tree::gen_str_wrt_seq(const vector<string> &names, string &str) const {
  str = "";
  return gen_str_wrt_seq(str, names, root_);
}

//----------------------------------------------------------------------

int
Tree::height() const {
  if(!root_) return 0;
  return root_->height();
}

//----------------------------------------------------------------------

int
Tree::breadth() const {
  if(!root_) return 0;
  return root_->breadth();
}

//----------------------------------------------------------------------

double
Tree::max_branch_length() const {
  if(!root_) return 0.0;
  return root_->max_branch_length();
}

//----------------------------------------------------------------------

double
Tree::total_branch_length() const {
  if(!root_) return 0.0;
  return root_->total_branch_length();
}

//----------------------------------------------------------------------

void
Tree::req_size(int text_width, int text_height, int &req_width, int &req_height) const {
  req_width = height()*20+text_width;
  req_height = (text_height+4) * (breadth() - 1) + text_height*4;
}

//----------------------------------------------------------------------

void
Tree::calc_position(int va_width, int va_height,
		    int text_width, int text_height,
		    bool use_branch_lengths) {
  int height_ = height();
  int breadth_ = breadth();

#if 0
  cout << "max:   " << max_branch_length() << '\n'
       << "total: " << total_branch_length() << '\n';
#endif
  
  flip_hor_axis_ = false;
  flip_vert_axis_ = false;
  
  double margin = text_height * 3;

#if 0
  cout << "text_width:" << text_width << '\t' << "text_height:" << text_height << '\n'
       << "tree_height:" << height_ << '\t' << "tree_breadth:" << breadth_ << '\n'
       << "va_width:" << va_width << '\t' << "va_height:" << va_height << '\n'
       << "margin:" << margin << '\n';
#endif

  double cy = margin;
  double dx = (va_width - text_width - 2*margin) / double(height_-1);
  double dy = 0.0;
  if(breadth_ > 1) {
    dy = (va_height - 2*margin) / double(breadth_-1);
    cy = margin;
  } else {
    dy = (va_height - 2*margin) / 2.0;
    cy = margin + dy;
  }

#if 0
  cout << "dx:" << dx << '\t' << "dy:" << dy << '\n';
#endif

  calc_position(root_, use_branch_lengths, dx, dy, margin, (va_width - text_width - 2*margin), int(va_width - margin - text_width + .5), va_width, va_height, text_height, 0.0, max_branch_length(), margin, cy);
}

//----------------------------------------------------------------------

void
Tree::calc_position(node_t *node, bool use_branch_lengths, double dx, double dy, double margin, double max_tree_width, int va_width_m, int va_width, int va_height, int text_height, double c_branch_depth, double max_branch_length, double cx, double &cy) {
  if(!node) return;

  c_branch_depth += node->_length;
  
  if(node->_children.empty()) {  // A leaf node
    node->_y = cy;
    cy += dy;
    if(has_lengths && use_branch_lengths) {
      node->_x = c_branch_depth / max_branch_length * max_tree_width + margin;
    } else {
      node->_x = va_width_m;
    }
    node->_min_x = va_width_m;
    node->_max_x = va_width;
    node->_min_y = node->_y - text_height/2.0;
    node->_max_y = node->_y + text_height/2.0;
    //cout << node->_x << '\t' << node->_y << '\n';
  } else {                       // A node with children
    //double min_x = va_width_m;
    node->_min_x = 1e6;
    node->_max_x = -1e6;
    node->_min_y = 1e6;
    node->_max_y = -1e6;
    list<node_t *>::iterator i;
    for(i = node->_children.begin(); i != node->_children.end(); i++) {
      calc_position(*i, use_branch_lengths, dx, dy, margin, max_tree_width, va_width_m, va_width, va_height, text_height, c_branch_depth, max_branch_length, cx+dx, cy);
      //cout << "a: " << (*i)->x() << endl;
      //min_x = min(min_x, (*i)->_x);
      node->_min_x = min(node->_min_x, (*i)->_x);
      node->_max_x = max(node->_max_x, (*i)->_max_x);
      node->_min_y = min(node->_min_y, (*i)->_min_y);
      node->_max_y = max(node->_max_y, (*i)->_max_y);
    }
    //cout << '[' << node->_min_x << ',' << node->_max_x << "] [" << node->_min_y << ',' << node->_max_y << "]\n";
    //node->_x = cx;


    if(has_lengths && use_branch_lengths) {
      node->_x = c_branch_depth / max_branch_length * max_tree_width + margin;
    } else {
      node->_x = node->_min_x - dx;
    }

    
    //node->_y = (node->_children.front()->y() + node->_children.back()->y()) / 2.0;
    node->_y = (node->_min_y + node->_max_y) / 2.0;
  }
}

//----------------------------------------------------------------------



/*

find out the direct parent of a child with child_id

return false, if the tree does not contain this id;

return true & -1 if the child_id is the root

return true & real parent id fi the child_id is the not the root

*/



bool 

Tree::parent_id(int child_id, int & parent_id) {

	if(root_) {

		return root_->parent_id(child_id, parent_id);

	}

	else {

		return false;

	}

}



//----------------------------------------------------------------------


void
Tree::shift(int x_shift, int y_shift) {
  if(root_) {
    shift(x_shift, y_shift, root_);
  }
}

//----------------------------------------------------------------------
//the following method is to genrate the tree profile using the ID. The format
// has to complye with the Gu_Zhang 2001 demo code's format. 
//<LChild> <ID> <Rchild> <Branch Length>
//The array has to be initlized before passing to this method
//the parameter of branch is a 1D array to hold the branch length
//The profile parameter is an 1D array to hold the topology relation
void Tree::generate_tree_profile(int * profile, double * branch, int total_leaves) {
	if(root_) {
		reorder_id();
		root_->generate_profile(profile, branch, total_leaves);
	}
}

void Tree::reorder_id() {
	int id = 1;
	if(root_) {
		id = root_->breadth()+1;
		root_->reorder_id(id);
	}
}
//----------------------------------------------------------------------

void
Tree::shift(int x_shift, int y_shift, node_t *tree) {
  if(!tree) return;
  
  list<node_t *>::iterator i;
  for(i = tree->_children.begin(); i != tree->_children.end(); i++) {
    shift(x_shift, y_shift, *i);
  }

  tree->_x += x_shift;
  tree->_y += y_shift;
  tree->_min_x += x_shift;
  tree->_max_x += x_shift;
  tree->_min_y += y_shift;
  tree->_max_y += y_shift;  
}

//----------------------------------------------------------------------

void
Tree::flip(bool hor_axis, bool vert_axis) {
  if(root_) {
    flip(hor_axis, vert_axis, root_);
  }

  if(hor_axis) {
    flip_hor_axis_ = !flip_hor_axis_;
  }
  if(vert_axis) {
    flip_vert_axis_ = !flip_vert_axis_;
  }
}


//----------------------------------------------------------------------

bool
Tree::flipNode(int id) {
	return true;
}

//----------------------------------------------------------------------
bool
Tree::flipNode() {
	int id;
	const node_t * r = this->getSelectedRoot();
	if(r == NULL)
		return false;
	else
		id = r->id();

	((node_t *)r)->flipNode();
	return true;
}

//----------------------------------------------------------------------

void
Tree::flip(bool hor_axis, bool vert_axis, node_t *tree) {
  if(!tree) return;
  
  list<node_t *>::iterator i;
  for(i = tree->_children.begin(); i != tree->_children.end(); i++) {
    flip(hor_axis, vert_axis, *i);
  }

  if(hor_axis) {
    tree->_y = root_->_max_y - (tree->_y - root_->_min_y);
  }

  if(vert_axis) {
    tree->_x = root_->_max_x - (tree->_x - root_->_x);
  }
}

//----------------------------------------------------------------------

Tree::Tree()
  : root_(NULL),
    ch(0),
    has_lengths(false),
    draw_root_(true),
    draw_style(square),
    draw_node_id(false),
    flip_hor_axis_(false),
    flip_vert_axis_(false)
{}

//----------------------------------------------------------------------

Tree::Tree(const Tree &t)
  : root_(NULL),
    ch(t.ch),
    has_lengths(t.has_lengths),
    filename_(t.filename_),
    draw_root_(t.draw_root_),
    draw_style(t.draw_style),
    draw_node_id(t.draw_node_id),
    flip_hor_axis_(t.flip_hor_axis_),
    flip_vert_axis_(t.flip_vert_axis_)
{
  if(t.root_) {
    root_ = new node_t(*t.root_);
  }
}

//----------------------------------------------------------------------

Tree::Tree(const node_t *n)
  : root_(NULL),
    ch(0),
    has_lengths(false),  //XXX
    draw_root_(true),
    draw_style(square),
    draw_node_id(false),
    flip_hor_axis_(false),
    flip_vert_axis_(false)
{
  if(n) {
    root_ = new node_t(*n);
  }
  assign_ids();
}

//----------------------------------------------------------------------

void
Tree::set_draw_style(draw_style_t ds) {
  draw_style = ds;
}

//----------------------------------------------------------------------

void
Tree::set_draw_node_id(bool dni) {
  draw_node_id = dni;
}

//----------------------------------------------------------------------

Tree::~Tree() {
  if(root_) delete root_;
}

//----------------------------------------------------------------------

Tree &
Tree::operator = (const Tree &t) {
  if(this != &t) {
    if(root_) delete root_;
    root_ = NULL;
    if(t.root_) {
      root_ = new node_t(*t.root_);
    }
    ch = t.ch;
    has_lengths = t.has_lengths;
    filename_ = t.filename_;
    
    draw_root_ = t.draw_root_;
    
    draw_style = t.draw_style;
    draw_node_id = t.draw_node_id;
    
    flip_hor_axis_ = t.flip_hor_axis_;
    flip_vert_axis_ = t.flip_vert_axis_;
  }
  return *this;
}

//----------------------------------------------------------------------

bool
Tree::draw_root() const {
  return draw_root_;
}

//----------------------------------------------------------------------

void
Tree::set_draw_root(bool draw_root) {
  draw_root_ = draw_root;
}

//----------------------------------------------------------------------

double
Tree::parse_length(string &tree_str) {
  int numerator = 0;
  int denominator = 0;

  do {
    ch = getc_nws(tree_str);
    if(isdigit(ch)) {
      numerator = numerator * 10 + (ch - '0');
      if(denominator) denominator *= 10;
    } else if(ch == '.') denominator = 1;
    else if(ch == '-') numerator = -numerator;
  } while(isdigit(ch) || ch == '.' || ch == '-');

  double len;
  if(denominator) {
    has_lengths = true;
    len = double(numerator)/double(denominator);
  } else {
    len = double(numerator);
  }

  return len;
}



//----------------------------------------------------------------------
bool
Tree::import(string &tree_str, node_t *parent) {
  ch = getc_nws(tree_str);

  if(ch == '(') {
    node_t *p = new node_t;

    if(parent) {
      parent->add_node(p);
    } else {
      if(root_) delete root_;
      root_ = p;
    }

    do {
      if(!import(tree_str, p)) return false;
    } while(ch != ')' && ch != EOF && ch != ';');

    if(ch == ')') {
      ch = getc_nws(tree_str);
      if(ch == ':') {
	has_lengths = true;
	p->_length = parse_length(tree_str);
      }
    }
  } else {
    string label;
    double len = 0;

    if(!parent) return false;

    do {
      if(ch == '_') ch = ' ';
      label += ch;
      ch = getc_nws(tree_str);
    } while(ch != EOF && ch != ':' && ch != '(' && ch != ';' && ch != ',' && ch != ')');

    if(ch == ':') {
      has_lengths = true;
      len = parse_length(tree_str);
    }

    parent->add_node(new node_t(label, len));
  }

  return true;
}

//----------------------------------------------------------------------

bool
Tree::import(string tree_str) {
  if(!import(tree_str, NULL)) return false;
  assign_ids();
  return true;
}

//----------------------------------------------------------------------
//20240301:modify tree load rule
bool
Tree::load(const string str,const string filename) {

  if(!import(str)) {
    return false;
  }

  filename_ = filename;

  return true;
}

//----------------------------------------------------------------------

bool
Tree::load_augment(const string filename) {
  FILE *file = fopen(filename.c_str(), "r");
  if(!file) return false;
  
  map<string, string> aug_map;
  string path, label;
  int c;
  int line_num = 1;
  int tab_num = 0;
  do {
    c = getc(file);
    if(c == '\n') {
      // In the case of a sparse file there may not be any data
      if(tab_num == 1) {
	//cout << path << '\t' << label << endl;
	aug_map[path] = label;
      }
      
      line_num++;
      tab_num = 0;
      path.erase();
      label.erase();
    } else if(line_num >= 1 && tab_num < 2) {
      if(c == '\t') {
	tab_num++;
	if(tab_num == 1) {
	  //str += ' ';
	} else if(tab_num == 2) {
	  //cout << path << '\t' << label << endl;
	  aug_map[path] = label;
	}
      } else if(c == '(' || c == ')') {
      } else {
	if(tab_num == 0) {
	  path += (char)c;
	} else if(tab_num == 1) {
	  label += (char)c;
	}
      }
    }
  } while(c != EOF);
  
  fclose(file);

  // traverse the tree replacing the numerical label with the string label
  augment(root_, aug_map);
  
  return true;
}

//----------------------------------------------------------------------

void
Tree::augment(node_t *node, const map<string, string> &aug_map) {
  if(!node) return;
  
  if(!node->_label.empty()) {
    map<string, string>::const_iterator i;
    i = aug_map.find(node->_label);
    if(i != aug_map.end()) {
      node->_label = i->second;
    }
  }

  list<node_t *>::const_iterator i;
  for(i = node->_children.begin(); i != node->_children.end(); i++) {
    augment(*i, aug_map);
  }  
}

//----------------------------------------------------------------------

bool
Tree::load_highlights(const string filename) {
  FILE *file = fopen(filename.c_str(), "r");
  if(!file) return false;
  
  set<string> highlights;
  
  const int line_len = 1024;
  char line[line_len];

  while(fgets(line, line_len, file)) {
    if(line[0] == '#') {
      continue;
    }
    str_clean_rn(line);
    highlights.insert(line);
  }
  
  fclose(file);

  highlight(highlights);
  
  return true;
}

//----------------------------------------------------------------------

void
Tree::highlight(const set<string> &highlights) {
  highlight(root_, highlights);
}

//----------------------------------------------------------------------

void
Tree::highlight(node_t *node, const set<string> &highlights) {
  if(!node) return;
  
  if(!node->_label.empty()) {
    set<string>::const_iterator i;
    i = highlights.find(node->_label);
    node->_highlight = i != highlights.end();
  }

  list<node_t *>::const_iterator i;
  for(i = node->_children.begin(); i != node->_children.end(); i++) {
    highlight(*i, highlights);
  }  
}

//----------------------------------------------------------------------

bool
Tree::save(const string filename) {
  FILE *file = fopen(filename.c_str(), "w");
  if(!file) return false;

  string str = gen_str();

  fprintf(file, "%s;\n", str.c_str());

  fclose(file);

  return true;
}

//----------------------------------------------------------------------

const node_t *
Tree::select(int x, int y, bool exclusive, int radius) {
  node_t *n = find(x, y, root_, radius);

  if(n) {
    // The new state should also dependent on the parent...
    bool new_state = !n->_selected;
    if(exclusive/* && !n->_selected*/) {
      unselect_all();
    }
    select(n, new_state);
  }

  return n;
}


//----------------------------------------------------------------------



bool 

Tree::parents(int child_id, std::vector<int> & parents) {

	if(!root_)

		return false;

	else {

		return root_->parents(child_id, parents);

	}

}



//----------------------------------------------------------------------

//return the depth of node with id

int

Tree::id_depth(int id) {

	if(!root_)

		return -1;

	else

		return root_->id_depth(id);

}




//----------------------------------------------------------------------

void
Tree::select(node_t *node, bool state) {
  if(!node) return;
  
  node->_selected = state;

  list<node_t *>::const_iterator i;
  for(i = node->_children.begin(); i != node->_children.end(); i++) {
    select(*i, state);
  }
}

//----------------------------------------------------------------------

void
Tree::select_all() {
  select(root_, true);
}

//----------------------------------------------------------------------

void
Tree::unselect_all() {
  select(root_, false);
}

//----------------------------------------------------------------------

node_t *
Tree::find(int x, int y, node_t *node, int radius) {
  if(node) {
    if(abs(int(node->_x - x + 0.5)) <= radius && abs(int(node->_y - y + 0.5)) <= radius) {
      return node;
    }
    list<node_t *>::const_iterator i;
    for(i = node->_children.begin(); i != node->_children.end(); i++) {
      node_t *n;
      n = find(x, y, *i, radius);
      if(n) return n;
    }
  }
  return NULL;
}

//----------------------------------------------------------------------

void
Tree::build_clades(vector<set<string> > &clades) const {
  build_clades(clades, root_);
}

//----------------------------------------------------------------------

void
Tree::build_clades(vector<set<string> > &clades, const node_t *tree) {
  if(!tree) return;
  
  if(!tree->_children.empty()) {
    vector<string> vln;
    tree->leaf_names(vln);
    
    set<string> sln;
    vector<string>::const_iterator i;
    for(i = vln.begin(); i != vln.end(); i++) {
      sln.insert(*i);
    }
    
    clades.push_back(sln);
    
    list<node_t *>::const_iterator j;
    for(j = tree->_children.begin(); j != tree->_children.end(); j++) {
      build_clades(clades, *j);
    }
  }
}

//----------------------------------------------------------------------

void
Tree::compare_tree(const Tree &tree2) {
  vector<string> t2l;
  tree2.leaf_names(t2l);
  vector<set<string> > t2c;
  tree2.build_clades(t2c);
  compare_tree(t2l, t2c, root_);
}

//----------------------------------------------------------------------

void
Tree::compare_tree(const vector<string> &leaves,
		   const vector<set<string> > &clades, node_t *tree) {
  if(!tree) return;

  tree->_highlight = false;

  if(!tree->_children.empty()) {
    vector<string> vln;
    tree->leaf_names(vln);
    
    set<string> sln;
    {
      vector<string>::const_iterator i;
      for(i = vln.begin(); i != vln.end(); i++) {
	sln.insert(*i);
      }
    }

    vector<set<string> >::const_iterator i;
    for(i = clades.begin(); i != clades.end(); i++) {
      if(sln == *i) {
	tree->_highlight = true;
	break;
      }
    }

    list<node_t *>::const_iterator j;
    for(j = tree->_children.begin(); j != tree->_children.end(); j++) {
      compare_tree(leaves, clades, *j);
    }
  } else {
    vector<string>::const_iterator i;
    for(i = leaves.begin(); i != leaves.end(); i++) {
      if(tree->_label == *i) {
	tree->_highlight = true;
	break;
      }
    }
  }
}  

//----------------------------------------------------------------------

void
Tree::compare_clades(const Tree &tree2) const {
  vector<set<string> > t2c;
  tree2.build_clades(t2c);
  compare_clades(t2c);
}

//----------------------------------------------------------------------

void
Tree::compare_clades(const vector<set<string> > &t2c) const {
  vector<set<string> > t1c;
  build_clades(t1c);
  
  int n = 0;

  vector<set<string> >::const_iterator i;
  vector<set<string> >::const_iterator j;

  for(i = t1c.begin(); i != t1c.end(); i++) {
    for(j = t2c.begin(); j != t2c.end(); j++) {
      if(*i == *j) {
	n++;
	break;
	
#if 0
	set<string>::const_iterator k;
	for(k=(*i).begin(); k!=(*i).end(); k++) {
	  cout << (*k) << ',';
	}
	cout << '\n';
	for(k=(*j).begin(); k!=(*j).end(); k++) {
	  cout << (*k) << ',';
	}
	cout << '\n';
#endif	
      }
    }
    //cout << "====================================\n";
  }

  //cout << t1c.size() << ' ' << t2c.size() << ' ' << n << endl;
}

//----------------------------------------------------------------------

bool
Tree::clade_topology(const set<string> &clade, string &tree) const {
  tree = "";
  set<string> total;
  leaf_names(total);
  return clade_topology(root_, total, clade, tree);
}

//----------------------------------------------------------------------

bool
Tree::clade_topology(node_t *node, const set<string> &total, const set<string> &sought_clade, string &tree) const {
  if(!node) return false;
  set<string> clade;
  node->leaf_names(clade);
  if(clade == sought_clade) {
    Tree tmp(node);
    tree = tmp.gen_str();
    return true;
  } else if(total.size() - clade.size() == sought_clade.size()) {
    set<string> tmp;
    set_difference(total.begin(), total.end(),
		   clade.begin(), clade.end(),
		   inserter(tmp, tmp.begin()));
    if(tmp == sought_clade) {
      if(!gen_str_pruned(node, tree)) abort();
      return true;
    }
  }
  list<node_t *>::const_iterator i;
  for(i = node->_children.begin(); i != node->_children.end(); i++) {
    if(clade_topology(*i, total, sought_clade, tree)) return true;
  }  
  return false;
}

//----------------------------------------------------------------------

bool
Tree::polyroot(Tree &tree) const {
  if(!root_) return false;
  
  const list<node_t*> &children = root_->children();
  
  if(children.size() == 2) {
    node_t n;
    const node_t *n1, *n2, *n3;
    if(children.front()->breadth() > children.back()->breadth()) {
      n1 = children.front()->children().front();
      n2 = children.front()->children().back();
      n3 = children.back();
    } else {
      n1 = children.front();
      n2 = children.back()->children().front();
      n3 = children.back()->children().back();
    }
    int nn1 = n1->breadth();
    int nn2 = n2->breadth();
    int nn3 = n3->breadth();
    if(nn1 >= nn2 && nn2 >= nn3) {
      n.add_node(new node_t(*n1));
      n.add_node(new node_t(*n2));
      n.add_node(new node_t(*n3));
    } else if(nn1 >= nn3 && nn3 >= nn2) {
      n.add_node(new node_t(*n1));
      n.add_node(new node_t(*n3));
      n.add_node(new node_t(*n2));
    } else if(nn2 >= nn1 && nn1 >= nn3) {
      n.add_node(new node_t(*n2));
      n.add_node(new node_t(*n1));
      n.add_node(new node_t(*n3));
    } else if(nn2 >= nn3 && nn3 >= nn1) {
      n.add_node(new node_t(*n2));
      n.add_node(new node_t(*n3));
      n.add_node(new node_t(*n1));
    } else if(nn3 >= nn1 && nn1 >= nn2) {
      n.add_node(new node_t(*n3));
      n.add_node(new node_t(*n1));
      n.add_node(new node_t(*n2));
    } else if(nn3 >= nn2 && nn2 >= nn1) {
      n.add_node(new node_t(*n3));
      n.add_node(new node_t(*n2));
      n.add_node(new node_t(*n1));
    }
    tree = &n;
  } else if(children.size() == 3) {
    tree = *this;
  } else {
    tree = *this;
  }

  return true;
}

//----------------------------------------------------------------------

void
Tree::assign_ids() {
  int internal_id = breadth() + 1;
  int leaf_id = 1;
  assign_ids(root_, internal_id, leaf_id);
}

//----------------------------------------------------------------------

void
Tree::assign_ids(node_t *node, int &internal_id, int &leaf_id) {
  if(node) {
    if(node->_children.empty()) {
      node->_id = leaf_id++;
    } else {
      node->_id = internal_id++;
      list<node_t *>::const_iterator i;
      for(i = node->_children.begin(); i != node->_children.end(); i++) {
	assign_ids(*i, internal_id, leaf_id);
      }
    }
  }
}

//----------------------------------------------------------------------

const node_t *
Tree::find_id(int id) const {
  return find_id(root_, id);
}

//----------------------------------------------------------------------

const node_t *
Tree::find_id(const node_t *node, int id) {
  if(node) {
    if(node->_id == id) return node;
    list<node_t *>::const_iterator i;
    for(i = node->_children.begin(); i != node->_children.end(); i++) {
      const node_t *n = find_id(*i, id);
      if(n) return n;
    }
  }
  return NULL;
}

//----------------------------------------------------------------------

int
Tree::id_lookup(string label) const {
  return id_lookup(root_, label);
}

//----------------------------------------------------------------------

int
Tree::id_lookup(const node_t *node, string label) {
  if(node) {
    if(node->_label == label) return node->_id;
    list<node_t *>::const_iterator i;
    for(i = node->_children.begin(); i != node->_children.end(); i++) {
      int id = id_lookup(*i, label);
      if(id != -1) return id;
    }
  }
  return -1;
}

//----------------------------------------------------------------------

string
Tree::id_lookup(int id) const {
  return id_lookup(root_, id);
}

//----------------------------------------------------------------------

string
Tree::id_lookup(const node_t *node, int id) {
  if(node) {
    if(node->_id == id) return node->_label;
    list<node_t *>::const_iterator i;
    for(i = node->_children.begin(); i != node->_children.end(); i++) {
      string label = id_lookup(*i, id);
      if(!label.empty()) return label;
    }
  }
  return "";
}

//----------------------------------------------------------------------
bool
Tree::reroot(int id) {
  if(!root_ || root_->id() == id) return true;

  //cout << gen_str() << endl;
  
  node_t *p = NULL, *nr = NULL;
  double p_len = 0.0;
  if(!reroot(root_, id, p, p_len, nr)) return false;

  root_->_children.clear();
  delete root_;
  root_ = nr;
  assign_ids();

  //cout << gen_str() << endl;
  
  return true;
}

//----------------------------------------------------------------------


bool
Tree::reroot() {
	int id;
	const node_t * r = this->getSelectedRoot();
	if(r == NULL)
		return false;
	else
		id = r->id();

  if(!root_ || root_->id() == id) return true;

  //cout << gen_str() << endl;
  
  node_t *p = NULL, *nr = NULL;
  double p_len = 0.0;
  if(!reroot(root_, id, p, p_len, nr)) return false;

  root_->_children.clear();
  delete root_;
  root_ = nr;
  assign_ids();

  //cout << gen_str() << endl;
  
  return true;
}

//----------------------------------------------------------------------

// Something may not be quite right about the lengths...
bool
Tree::reroot(node_t *node, int id, node_t *&p, double &p_len, node_t *&nr) {
  if(node) {
    if(node->_id == id) {
      nr = new node_t;
      nr->add_node(node);
      p_len = node->_length = node->_length/2.0;
      p = nr;
      return true;
    }
    list<node_t *>::iterator i;
    for(i = node->_children.begin(); i != node->_children.end(); i++) {
      if(reroot(*i, id, p, p_len, nr)) {
	if(node != root_) {
	  double tmp = node->_length;
	  node->_length = p_len;
	  p_len = tmp;
	  node->remove_node(*i);
	  p->add_node(node);
	} else {
	  list<node_t *>::iterator j;
	  for(j = node->_children.begin(); j != node->_children.end(); j++) {
	    if(j != i) {
	      //(*j)->_length += (*i)->_length / 2.0;
	      (*j)->_length += p_len;
	      p->add_node(*j);
	    }
	  }
	  //(*i)->_length /= 2.0;
	  //(*i)->_length += p_len;
	}
	p = node;
	return true;
      }
    }
  }
  return false;
}

//======================================================================

#ifdef QT

#include <qpainter.h>

//----------------------------------------------------------------------

QtTree::QtTree() {
  font_= QFont("Times New Roman", 10);
  tipWidth = 1;
}

//----------------------------------------------------------------------

QtTree::QtTree(const Tree &t)
  : Tree(t)
{}

//----------------------------------------------------------------------

QtTree::~QtTree() {}

//----------------------------------------------------------------------

QtTree &
QtTree::operator = (const Tree &t) {
  Tree::operator = (t);
  return *this;
}

//----------------------------------------------------------------------

void
QtTree::draw(QPainter *p, int sx, int sy, int ex, int ey) const {
  if(!root_) return;

  p->save();
  
  if(draw_root_) {
    if(root_->_selected) {
      p->setPen(Qt::red);
      p->setBrush(Qt::red);
    } else {
      if(root_->_highlight) {
	p->setPen(Qt::blue);
	p->setBrush(Qt::blue);
      } else {
	p->setPen(Qt::black);
	p->setBrush(Qt::black);
      }
    }

    if(flip_vert_axis_) {
      p->drawLine(int(root_->x()+10+.5), int(root_->y()+.5),
		  int(root_->x()+.5),    int(root_->y()+.5));
    } else {
      p->drawLine(int(root_->x()-10+.5), int(root_->y()+.5),
		  int(root_->x()+5),     int( root_->y()+.5));
    }
  }

  draw(p, sx, sy, ex, ey, root_);

  p->restore();
}

//----------------------------------------------------------------------

void
QtTree::draw(QPainter *p, int sx, int sy, int ex, int ey, const node_t *node) const {
  if(!node) return;

#if 1
  if(!flip_hor_axis_ && !flip_vert_axis_) {
    if(node->_x > ex ||
       node->_min_y > ey ||
       node->_max_x < sx ||
       node->_max_y < sy) {
      //cout << "e: " << '[' << sx << ',' << sy << "] [" << ex << ',' << ey << "]   " << '[' << node->_min_x << ',' << node->_max_x << "] [" << node->_min_y << ',' << node->_max_y << "]\n";
      return;
    }
  }
#else
  if((node->_x > ex || flip_vert_axis_ && node->_x < ex) ||
     (node->_min_y > ey || flip_hor_axis_ && node->_min_y < ey) ||
     (node->_max_x < sx || flip_vert_axis_ && node->_max_x > sx) ||
     (node->_max_y < sy || flip_hor_axis_ && node->_max_y > sy)) {
    return;
  }
#endif

  p->save();

  if(node->_selected) {
    p->setPen(Qt::red);
    p->setBrush(Qt::red);
  } else {
    if(node->_highlight) {
      p->setPen(Qt::blue);
      p->setBrush(Qt::blue);
    } else {
      p->setPen(Qt::black);
      p->setBrush(Qt::black);
    }
  }

   QPen pen = p->pen();
   pen.setWidth(tipWidth);
   p->setPen(pen);

  p->fillRect(int(node->x()-tipWidth/2-.5), int(node->y()-tipWidth/2-.5), tipWidth + 2, tipWidth+2, p->brush());
  //p->drawPie(int(node->x()-3+.5), int(node->y()-3+.5), 6, 6, 0, 5760);
  
  p->setFont(font_);

  if(node->_children.empty()) {
    //cout << "writing: " << node->_label << '\t' << node->_selected << '\n';
    QRect r = p->fontMetrics().boundingRect(node->_label.c_str());
    if(flip_vert_axis_) {
      p->drawText(int(node->x()-5-p->fontMetrics().width(node->_label.c_str())+.5), int(node->y()+r.height()/2.0+.5), node->_label.c_str());
    } else {
      p->drawText(int(node->x()+5+.5), int(node->y()+r.height()/2.0+.5),
		  node->_label.c_str());
    }
  } else {
    if(draw_style == square) {
      p->drawLine(int(node->x()+.5), int(node->_children.front()->y()-tipWidth/2+.5),
		  int(node->x()+.5), int(node->_children.back()->y()+tipWidth/2+.5));
    }
    list<node_t *>::const_iterator i;
    for(i = node->_children.begin(); i != node->_children.end(); i++) {
      if(draw_style == square) {
	p->drawLine(int(node->x()+.5), int((*i)->y()+.5),
		    int((*i)->x()+.5), int((*i)->y()+.5));
      } else if(draw_style == curve) {
	int x = int(node->x()+.5);
	int y = int(node->y()+.5);
	int w = int((*i)->x()-node->x()+.5);
	int h = int((*i)->y()-node->y()+.5);
	if(h > 0) {
	  p->drawArc(x, y-h, w*2, h*2, 16*180, 16*90);
	} else {
	  p->drawArc(x, y-h, w*2, h*2, 16*180, -16*90);
	}
      }
      draw(p, sx, sy, ex, ey, *i);
    }
    if(draw_node_id) {
      QString str;
      str.setNum(node->id());
      p->drawText(int(node->x()+.5)-15, int(node->y()+.5), str);
    }
  }
  p->restore();
}

//----------------------------------------------------------------------

void
QtTree::max_label_dim(const QFontMetrics &fm, int &width, int &height) const {
  width = 0;
  height = 0;
  max_label_dim(root_, fm, width, height);
}

//----------------------------------------------------------------------

void
QtTree::max_label_dim(const node_t *node, const QFontMetrics &fm, int &width, int &height) const {
  if(!node) return;
  
  if(node->_children.empty()) {
    QRect r = fm.boundingRect(node->_label.c_str());
    width = r.width() - r.x();
    height = r.height();
  } else {
    list<node_t *>::const_iterator i;
    for(i = node->_children.begin(); i != node->_children.end(); i++) {
      int w = 0, h = 0;
      max_label_dim(*i, fm, w, h);
      width = max(width, w);
      height = max(height, h);
    }
  }
}

const QFont & QtTree::getFont() const {
  return font_;

}

void QtTree::setFont(QFont f) {
  font_ = f;
}
#endif

//======================================================================
