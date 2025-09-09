#ifndef _CLUSTER_H_
#define _CLUSTER_H_

#include <string>
#include <vector>

#include "sequence.h"

typedef enum { pDistance, Poisson, Kimura } clustering_method_t;

void nj_cluster(std::vector<sequence_t> seqs,
		std::string &tree,
		clustering_method_t matrix_id,
		bool ignore_gaps = true,
		bool branch_lengths = true);

void nj_bootstrap(std::vector<sequence_t> seqs,
		  std::vector<std::string> &trees,
		  clustering_method_t matrix_id,
		  int nsamples,
		  bool ignore_gaps = true);

void sample(const std::vector<sequence_t> &seqs, std::vector<sequence_t> &ss);

#endif
