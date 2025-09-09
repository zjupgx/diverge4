#ifndef _GZ97_H_
#define _GZ97_H_

#include "matrix.h"

bool gz97_compute(bool jttF,
	double alpha,
	const CMatrix2D& alignment,
	int treeE,
	const IMatrix2D& tree,
	const DVector& freq,
	const DMatrix2D& prob,
	DVector& ret);

bool gz97_compute_adapter_for_rvs(bool jttF,
	double alpha,
	const CMatrix2D& alignment,
	int treeE,
	const IMatrix2D& tree,
	const DVector& freq,
	const DMatrix2D& prob,
	DVector& ret,
	double& inferred_alpha);

bool ancestral_inference(bool jttF,
	double alpha,
	const CMatrix2D& alignment,
	int treeE,
	const IMatrix2D& tree,
	const DVector& freq,
	const DMatrix2D& prob,
	CMatrix2D& alignmentOut);

bool parse_tree(const char* tree_str,
	int seqNum,
	int& treeE,
	IMatrix2D& tree);

#endif
