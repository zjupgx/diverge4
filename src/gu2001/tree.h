#ifndef _TREE_H_
#define _TREE_H_

#include <string>
#include <list>
#include <map>
#include <set>
#include <vector>
#include <iosfwd>

#include <stdio.h>

class TreeView;

class node_t {
public:
  node_t();
  node_t(const node_t &n);
  node_t(const std::string &label, double length=0.0);
  ~node_t();

  node_t &operator = (const node_t &n);

  void add_node(node_t *n);
  bool remove_node(node_t *n); // Does not delete the node
  int height() const;
  int n_nodes() const;
  int n_children() const;
  int breadth() const;
  double max_branch_length() const;
  double total_branch_length() const;
  int id() const; 
  void leaf_names(std::vector<std::string> &names) const;
  void leaf_names(std::set<std::string> &names) const;
  
  //void print(int indent = 0) const;

  double x() const { return _x; }
  double y() const { return _y; }

  double min_x() const { return _min_x; }
  double max_x() const { return _max_x; }
  
  double min_y() const { return _min_y; }
  double max_y() const { return _max_y; }

  bool selected() const { return _selected; }
  bool highlight() const { return _highlight; }

  std::string label() const { return _label; }
  double length() const { return _length; }
  const std::list<node_t *> &children() const { return _children; }

  void generate_profile(int * profile, double * branch, int total_leaves);
  void reorder_id(int & id);

  bool parent_id(int child_id, int &parent_id);

  bool parents(int child_id, std::vector<int> & parents);

  bool contain_id(int id);

  int id_depth(int id);

  
  bool flipNode();

private:
  std::list<node_t *> _children;
  std::string _label;
  double _length;  // branch length for this node
  double _x, _y;
  double _min_x, _max_x;
  double _min_y, _max_y;
  bool _selected;
  bool _highlight;
  int _id;
  
  friend class Tree;
  friend class QtTree;
  friend class PSTree;
};

class Tree {
public:
  Tree();
  Tree(const Tree &t);
  Tree(const node_t *n);
  ~Tree();

  Tree &operator = (const Tree &t);

  bool import(std::string tree_str);
  
  bool load(const std::string str,const std::string filename);
  bool save(const std::string filename);

  bool load_augment(const std::string filename);
  bool load_highlights(const std::string filename);
  void highlight(const std::set<std::string> &highlights);
  bool has_branch_length() {return has_lengths;}
  void set_has_branch_length(bool has_brch) {has_lengths = has_brch;}
  int height() const;  // height of tree including leaves
  int breadth() const; // number of leaves
  double max_branch_length() const;
  double total_branch_length() const;

  bool draw_root() const;
  void set_draw_root(bool);
  
  void req_size(int text_width, int text_height, int &req_width, int &req_height) const;
  void calc_position(int va_width, int va_height, int text_width, int text_height, bool use_branch_lengths = true);
  void shift(int x_shift, int y_shift);
  void flip(bool hor_axis, bool vert_axis);

  const node_t *select(int x, int y, bool exclusive=false, int radius=3);
  void select_all();
  void unselect_all();

  typedef enum { square, curve } draw_style_t;
  void set_draw_style(draw_style_t);

  void set_draw_node_id(bool);

  void set_filename(std::string filename)  { filename_ = filename; }
  std::string filename() const { return filename_; }

  void leaf_names(std::vector<std::string> &names) const;
  void leaf_names(std::set<std::string> &names) const;
  std::string gen_str(bool include_lengths = true) const;
  bool gen_str_pruned(const node_t *node, std::string &tree, bool include_lengths = true) const;
  bool gen_str_wrt_seq(const std::vector<std::string> &names, std::string &str) const;

  void compare_tree(const Tree &tree2);
  void compare_clades(const Tree &tree2) const;
  void compare_clades(const std::vector<std::set<std::string> > &clades) const;
  void build_clades(std::vector<std::set<std::string> > &clades) const;
  bool clade_topology(const std::set<std::string> &clade, std::string &tree) const;

  bool polyroot(Tree &tree) const;
  
  const node_t *find_id(int id) const;
  int id_lookup(std::string label) const;
  std::string id_lookup(int id) const;
  
  const node_t *root() const { return root_; }
  bool reroot(int id);
  bool reroot();

  void generate_tree_profile(int * profile, double * branch, int total_leaves);
  void reorder_id();

  bool flipNode(int id);
  bool flipNode();

  bool parent_id(int child_id, int & parent_id);

  bool parents(int id, std::vector<int> &parents);

  int id_depth(int id);




  const node_t *getSelectedRoot() const;
  const node_t *getSelectedRoot(node_t * node) const;

  void assign_ids();

protected:
  node_t *root_;
  int ch;
  bool has_lengths;
  std::string filename_;

  bool draw_root_;
  
  draw_style_t draw_style;
  bool draw_node_id;

  bool flip_hor_axis_, flip_vert_axis_;
  
  bool import(std::string &tree_str, node_t *parent);
  double parse_length(std::string &tree_str);
  void calc_position(node_t *node, bool use_branch_lengths,
		     double dx, double dy, double margin,
		     double max_tree_width, int va_width_m, int va_width,
		     int va_height, int text_height,
		     double c_branch_length, double max_branch_length,
		     double cx, double &cy);
  void shift(int x_shift, int y_shift, node_t *tree);
  void flip(bool hor_axis, bool vert_axis, node_t *tree);
  void select(node_t *node, bool state);
  node_t *find(int x, int y, node_t *node, int radius=3);
  void augment(node_t *node, const std::map<std::string, std::string> &aug_map);
  void highlight(node_t *node, const std::set<std::string> &highlights);
  void gen_str(std::string &str, bool include_lengths, const node_t *tree) const;
  bool gen_str_pruned(std::string &str, bool include_lengths, const node_t *node, const node_t *prune_node) const;
  bool gen_str_wrt_seq(std::string &str, const std::vector<std::string> &names, const node_t *tree) const;
  static void build_clades(std::vector<std::set<std::string> > &clades, const node_t *tree);
  static void compare_tree(const std::vector<std::string> &leaves, const std::vector<std::set<std::string> > &clades, node_t *tree);
 
  static void assign_ids(node_t*, int&, int&);
  static const node_t *find_id(const node_t *node, int id);
  static int id_lookup(const node_t *node, std::string label);
  static std::string id_lookup(const node_t *node, int id);
  bool reroot(node_t *node, int id, node_t *&p, double &p_len, node_t *&nr);
  bool clade_topology(node_t *node, const std::set<std::string> &total, const std::set<std::string> &sought_clade, std::string &tree) const;

};

#ifdef QT

#include <qfontmetrics.h>

class QPainter;

class QtTree : public Tree {
public:
  QtTree();
  QtTree(const Tree &t);
  ~QtTree();

  QtTree &operator = (const Tree &t);

  void draw(QPainter *p, int sx, int sy, int ex, int ey) const;
  void max_label_dim(const QFontMetrics &fm, int &width, int &height) const;
  void max_label_dim(const node_t *node, const QFontMetrics &fm, int &width, int &height) const;

  const QFont & getFont() const;
  void setFont(QFont f);
  int getTipWidth() {return tipWidth;}
  void setTipWidth(int w) {tipWidth = w;}

protected:
  QFont font_;
  int tipWidth;
  void draw(QPainter *p, int sx, int sy, int ex, int ey, const node_t *node) const;

  friend class TreeView;
};

#endif

#endif
