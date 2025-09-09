#define _CRT_SECURE_NO_WARNINGS

#include <string>
#include <vector>

#include <math.h>
#include <ctype.h>

#include <stdio.h>

#ifndef WIN32
#include <values.h>
#else
#include <float.h>
#endif

#include "cluster.h"
#include "sequence.h"

using namespace std;

#ifndef WIN32
#  define isfinite(x) \
     (sizeof (x) == sizeof (float)                                            \
      ? __finitef (x)                                                         \
      : sizeof (x) == sizeof (double)                                         \
      ? __finite (x) : __finitel (x))
#else
#define isfinite(x) _finite(x)
#endif

#if 0
static void
print_matrix(const double *d, size_t dim, size_t n, const vector<string> &labels) {
  for(size_t i=0; i<n; i++) {
    cout << '\t' << labels[i];
  }
  cout << '\n';
  for(size_t i=0; i<n; i++) {
    cout << labels[i];
    for(size_t j=0; j<=i; j++) {
      cout << "\t--";
    }
    for(size_t j=i+1; j<n; j++) {
      cout << '\t' << d[i * dim + j];
    }
    cout << '\n';
  }
}
#endif

static void
construct_distance_matrix(const vector<sequence_t> &seqs,
			  clustering_method_t method,
			  bool ignore_gaps, double *dist, double &max_dist) {
  max_dist = 0.0;
  size_t n = seqs.size();
  if(n == 0) return;
  size_t len = seqs.front().sequence.size();
  for(size_t s1=0; s1<n-1; s1++) {
    for(size_t s2=s1+1; s2<n; s2++) {
      int d=0;
      for(size_t i=0; i<len; i++) {
	if(seqs[s1].sequence[i] != seqs[s2].sequence[i]) d++;
      }
      double p = d / double(len);
      double q;
      double dij = 0.0;
      switch(method) {
      case pDistance:
	dij = p;
	break;
      case Poisson:
	dij = -log(1.0 - p);
	break;
      case Kimura:
	q = 1.0 - p - 0.2 * p*p;
	if(q < 0.0) q = 0.0;
	dij = -log(q);
	//printf("%f %f\n", p, dist[s1 * n + s2]);
	break;
      }
      if(dij <= 0.0) dij = 0.0;
      dist[s1 * n + s2] = dij;
      if(isfinite(dij) && max_dist < dij) max_dist = dij;
    }
  }
}

static void
nj_cluster(vector<string> &labels, double *d, double max_dist, double *r, size_t dim, size_t n,
	   string &tree, bool branch_lengths) {
  if(n == 0) {
    return;
  } else if(n == 1) {
    tree = labels[0];
    return;
  } else if(n == 2) {
    if(branch_lengths) {
      double dh = d[0 * dim + 1]/2.0;
      if(!isfinite(dh)) dh = max_dist;
      if(dh <= 0.0) dh = 0.0;
      char di[16], dj[16];
      sprintf(di, "%f", dh);
      sprintf(dj, "%f", dh);
      tree = '(' + labels[0] + ':' + di + ',' + labels[1] + ':' + dj + ')';
    } else {
      tree = '(' + labels[0] + ',' + labels[1] + ')';
    }
    return;
  }
  
  //print_matrix(d, dim, n, labels);

  size_t i, j, i2, j2, a, b;
  int min_i = -1, min_j = -1;

  for(i=0; i<n; i++) {
    r[i] = 0.0;
    for(j=0; j<n; j++) {
      if(i < j) {
	r[i] += d[i * dim + j];
      } else if(i > j) {
	r[i] += d[j * dim + i];
      }
    }
  }

#if 0
  cout << "R:";
  for(i=0; i<n; i++) {
    cout << '\t' << r[i];
  }
  cout << '\n';
#endif
  
  double min_m = DBL_MAX;
  for(i=0; i<n-1; i++) {
    for(j=i+1; j<n; j++) {
      double m = d[i * dim + j] - (r[i] + r[j]) / (n - 2);
      //cout << i << '\t' << j << '\t' << m << '\n';
      if(isfinite(m) && m < min_m) {
	min_m = m;
	min_i = i;
	min_j = j;
      }
    }
  }

  if(min_i == -1 && min_j == -1) {
    //cout << "min = " << min_m << " at " << min_i << ' ' << min_j << endl;
    min_i = 0;
    min_j = 1;
    min_m = max_dist;
  }

  for(a=b=min_i; a<n; a++) {
    if(a == min_i) {
      if(branch_lengths) {
	double vi = d[min_i * dim + min_j] / 2.0 + (r[min_i] - r[min_j]) / (2 * (n - 2));
	double vj = d[min_i * dim + min_j] - vi;
	if(!isfinite(vi)) vi = max_dist;
	if(!isfinite(vj)) vj = max_dist;
	if(vi <= 0.0) vi = 0.0;
	if(vj <= 0.0) vj = 0.0;
	char di[16], dj[16];
	sprintf(di, "%f", vi);
	sprintf(dj, "%f", vj);
	//printf("%f\t%d\t%s\n", vi, vi <= 0.0, di);
	//printf("%f\t%d\t%s\n", vj, vj <= 0.0, dj);
	labels[b++] = '(' + labels[min_i] + ':' + di + ',' + labels[min_j] + ':' + dj + ')';
      } else {
	labels[b++] = '(' + labels[min_i] + ',' + labels[min_j] + ')';
      }
    } else if(a != min_j) {
      labels[b++] = labels[a];
    }
  }

  double dij = d[min_i * dim + min_j];
  if(!isfinite(dij)) dij = max_dist;
  for(i2=i=0; i<n-1; i++) {
    if(i == min_i) {
      for(j2=i2+1, j=i+1; j<n; j++) {
	if(j < min_j) {
	  d[i2 * dim + j2] = (d[i * dim + j] + d[j * dim + min_j] - dij) / 2.0;
	  j2++;
	} else if(j > min_j) {
	  d[i2 * dim + j2] = (d[i * dim + j] + d[min_j * dim + j] - dij) / 2.0;
	  j2++;
	}
      }
      i2++;
    } else if(i != min_j) { // i != min_i && i != min_j
      for(j2=i2+1, j=i+1; j<n; j++) {
	if(j == min_i) {
	  d[i2 * dim + j2] = (d[i * dim + j] + d[i * dim + min_j] - dij) / 2.0;
	  j2++;
	} else if(j != min_j) { // j != min_i && j != min_j
	  d[i2 * dim + j2] = d[i * dim + j];
	  j2++;
	}
      }
      i2++;
    }
  }
  
  nj_cluster(labels, d, max_dist, r, dim, n-1, tree, branch_lengths);
}

static void
preprocess(vector<sequence_t> &seqs, double *dist, double &max_dist, vector<string> &labels,
	   clustering_method_t method, bool ignore_gaps, bool branch_lengths) {
  if(ignore_gaps) clean_gaps(seqs);

  int s, i, j;

  for(s=0; s<(int)seqs.size(); s++) {
    int len = seqs[s].sequence.size();
    for(i=0; i<len; i++) {
      seqs[s].sequence[i] = toupper(seqs[s].sequence[i]);
    }
  }

  construct_distance_matrix(seqs, method, ignore_gaps, dist, max_dist);

  labels.resize(seqs.size());
  for(i=0; i<(int)seqs.size(); i++) {
    for(j=0; j<(int)seqs[i].label.size(); j++) {
      if(isspace(seqs[i].label[j])) seqs[i].label[j] = '_';
    }
    labels[i] = seqs[i].label;
  }
}
	   
void
nj_cluster(vector<sequence_t> seqs, string &tree,
	   clustering_method_t method, bool ignore_gaps, bool branch_lengths) {
  tree = "";
  size_t dim = seqs.size();
  if(dim == 0) return;
  double *dist = new double[dim * dim];
  double *r = new double[dim];
  double max_dist;
  vector<string> labels;
  preprocess(seqs, dist, max_dist, labels, method, ignore_gaps, branch_lengths);
  nj_cluster(labels, dist, max_dist, r, dim, dim, tree, branch_lengths);
  delete[] dist;
  delete[] r;
}

void
sample(const vector<sequence_t> &seqs, vector<sequence_t> &ss) {
  if(seqs.empty()) return;
  int len = seqs.front().sequence.size();
  int num = seqs.size();
  ss.resize(num);
  int i, j;
  for(j=0; j<num; j++) {
    ss[j].label = seqs[j].label;
    ss[j].sequence.resize(len);
  }
  for(i=0; i<len; i++) {
#ifndef WIN32
    int pos = random() % len;
#else
    int pos = rand() % len;
#endif
    for(j=0; j<num; j++) {
      ss[j].sequence[i] = seqs[j].sequence[pos];
    }
  }
}

void
nj_bootstrap(vector<sequence_t> seqs, vector<string> &trees,
	     clustering_method_t method, int nsamples, bool ignore_gaps) {
  size_t dim = seqs.size();
  if(dim == 0) return;
  double *dist = new double[dim * dim];
  double *r = new double[dim];
  double max_dist;
  vector<string> labels;
  vector<sequence_t> ss;
  for(int i=0; i<nsamples; i++) {
    //cout << "I: " << i << endl;
    sample(seqs, ss);
    //print_seqs(ss);
    preprocess(ss, dist, max_dist, labels, method, ignore_gaps, false);
    string ss_tree;
    nj_cluster(labels, dist, max_dist, r, dim, dim, ss_tree, false);
    //cout << ss_tree << '\n';
    trees.push_back(ss_tree);
  }
  delete[] dist;
  delete[] r;
}
