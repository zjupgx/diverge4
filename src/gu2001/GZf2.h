#ifndef _GZF2_H_
#define _GZF2_H_

#include <string>
#include <vector>

#include "matrix.h"

bool GZf2_compute(const std::vector<DVector> &numSub,
		  std::vector<std::vector<double> > &summary,
		  std::vector<std::vector<double> > &rets);

#endif
