#ifndef _ANCESTRAL_COMPUTE_H_
#define _ANCESTRAL_COMPUTE_H_

#include <vector>
#include <string>
#include "tree.h"
#include "sequence.h"


bool ancestral_compute(Tree &tree,
		       const std::vector<sequence_t> &sequences,
			   std::vector<sequence_t> &sequencesOut, std::string & warning);


#endif
 