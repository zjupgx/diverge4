#ifndef _TYPE_TWO_H_
#define _TYPE_TWO_H_

#include <vector>
#include <string>

#include "tree.h"
#include "sequence.h"

bool
type_two_compute(const std::vector<Tree> &trees, const std::vector<sequence_t> &sequences,
		 std::vector<std::vector<double> > &summary, std::vector<std::vector<double> > &rets2);

// New function for gap-weighted posterior probability calculation
bool
type_two_compute_with_gaps(const std::vector<Tree> &trees, 
                          const std::vector<sequence_t> &clean_sequences,
                          const std::vector<sequence_t> &original_sequences,  
                          const std::vector<int> &kept_positions,
                          std::vector<std::vector<double> > &summary, 
                          std::vector<std::vector<double> > &rets2);

#endif
