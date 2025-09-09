#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gz97.h"
#include "matrix.h"


/*----------------------------------------------------------------------*/
static double gz97_inferred_alpha;
static bool common_compute(bool jttF, double alpha, const CMatrix2D &alignment,
			   int treeE, const IMatrix2D &tree,
			   const DVector &freq_, const DMatrix2D &prob_, BMatrix2D &seqN,
			   DVector &length, IMatrix2D &branches, IVector &parent,
			   IMatrix2D &best);
static void estimate_branch_lengths(int seqNum, int treeE,
				    const DMatrix2D &distance,
				    const IMatrix2D &tree,
				    IMatrix2D &branches,
				    IVector &parent, DVector &length);
static void compute_trans_matrix(int seqNum, const DMatrix2D &prob,
				 DMatrix3D &probability, DVector &length);
static bool reconstruct(int treeE,
			const BMatrix2D &seqN,
			const DMatrix3D &probability,
			const IMatrix2D &branches,
			const IMatrix2D &tree,
			const DVector &freqObs, const DVector &freq,
			IMatrix2D &best);
static void binaryGo(int branchI, int groupI, int j, int seqNum,
		     IMatrix3D &group, IMatrix2D &noOfGroup,
		     const IMatrix2D &tree);
static void matrix(const DMatrix2D &prob, DMatrix2D &transition);
static int change(char amnio);
static char inv_change(int amnio);
static double totalProb(int site, int treeE, const BMatrix2D &seqN,
			const DMatrix3D &probability,
			const IMatrix2D &branches,
			const IMatrix2D &tree, 
			const DVector &freq);
static void jttf(const BMatrix2D &seqN,
		 DVector &freqObs, DVector &freq,
		 DMatrix2D &prob);
static void stepDistance(int treeE, const IMatrix2D &best,
			 const BMatrix2D &seqN,
			 const DVector &length,
			 const IMatrix2D &branches,
			 DVector &distStep1);
static void gammaDistanceNew1(const BMatrix2D &seqN,
			      const DVector &length, const DVector &distStep1);
static double gammln(double xx);
static double sullivan(double mean, double z, const DVector &distStep);
static void rateF(double alpha, const BMatrix2D &seqN,
		  const DVector &length, const DVector &distStep1);
static double nodeByNode(int site, int node, int xxx, int treeE,
			 const BMatrix2D &seqN,
			 const DMatrix3D &probability,
			 const IMatrix2D &branches,
			 const IMatrix2D &tree,
			 const DVector &freq);

static const char aminoAcid[20] = { 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
				    'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V' };

/*----------------------------------------------------------------------*/
// number of sites N = rets.size()
bool
gz97_compute_adapter_for_rvs(bool jttF, double alpha, const CMatrix2D &alignment,
			     int treeE, const IMatrix2D &tree, const DVector &freq, const DMatrix2D &prob, 
			     DVector &ret, double &inferred_alpha) {
  bool result = gz97_compute(jttF, alpha, alignment, treeE, tree, freq, prob, ret);
  inferred_alpha = gz97_inferred_alpha;
  return result;
}
/*----------------------------------------------------------------------*/

bool
gz97_compute(bool jttF, double alpha, const CMatrix2D &alignment,
	     int treeE, const IMatrix2D &tree,
	     const DVector &freq, const DMatrix2D &prob, DVector &ret) {
  int seqNum = alignment.rows();
  int seqLength = alignment.cols();
  BMatrix2D seqN(seqNum, seqLength);

  DVector length;
  IMatrix2D branches;   // (x,0) = Left, (x,1) = Right
  IVector parent;
  IMatrix2D best;
  
  if(!common_compute(jttF, alpha, alignment, treeE, tree, freq, prob, seqN, length, branches, parent, best)) {
    return false;
  }

  parent.clear();
  
  DVector distStep1(seqLength);

  stepDistance(treeE, best, seqN, length, branches, distStep1);
  gammaDistanceNew1(seqN, length, distStep1);

  ret = distStep1;
  
  return true;
}

/*----------------------------------------------------------------------*/

// `common_compute` 一系列的预处理和计算工作,主要包括以下几个步骤:

// 1. 计算序列之间的 Gamma 距离:
//    - 遍历所有序列对,计算每对序列之间的距离。
//    - 使用 Gamma 距离公式计算距离值,其中 `alpha` 是 Gamma 分布的参数。

// 2. 估计进化树的分支长度:
//    - 调用 `estimate_branch_lengths` 函数,根据序列之间的距离和进化树拓扑结构,计算出每个分支的长度。
//    - 同时计算出每个节点的父节点信息 `parent`。

// 3. 构建序列状态矩阵 `seqN`:
//    - 根据输入的多序列比对 `alignment`,将每个序列中的氨基酸或核苷酸编码为数字。
//    - 将编码后的序列状态存储在 `seqN` 矩阵中。

// 4. 计算频率和概率矩阵:
//    - 如果 `jttF` 为 true,则调用 `jttf` 函数计算观察频率 `freqObs` 和概率矩阵 `prob`。
//    - 否则,使用输入的频率 `freq_` 和概率矩阵 `prob_`。
//    - 然后,调用 `compute_trans_matrix` 函数计算转移概率矩阵 `probability`。

// 5. 重构最优分支分配:
//    - 调用 `reconstruct` 函数,根据进化树拓扑结构、序列状态和概率矩阵,计算出最优的分支分配情况 `best`。


bool
ancestral_inference(bool jttF,
		    double alpha,
		    const CMatrix2D &alignment,
		    int treeE,
		    const IMatrix2D &tree,
		    const DVector &freq,
		    const DMatrix2D &prob,
		    CMatrix2D &alignmentOut) {
  int i, j;
  int seqNum = alignment.rows();
  int seqLength = alignment.cols();
  BMatrix2D seqN(seqNum, seqLength);

  DVector length;
  IMatrix2D branches;   // (x,0) = Left, (x,1) = Right
  IVector parent;
  IMatrix2D best;
  
  if(!common_compute(jttF, alpha, alignment, treeE, tree, freq, prob, seqN, length, branches, parent, best)) {
    return false;
  }

  int n = tree.rows();

#ifdef DEBUG_DIVERGE
  for(i=0; i<n; i++) {
    printf("%d\t%d\t%d\n", tree(i, 0), i, tree(i, 1));
  }
#endif

  alignmentOut.resize(seqNum*2, seqLength);
  for(i = 0; i < seqNum; i++) {
    for(j = 0; j < seqLength; j++) {
      alignmentOut(i, j) = inv_change(seqN(i, j));
    }
  }
  for(i = 0; i < seqNum; i++) {
    for(j = 0; j < seqLength; j++) {
      alignmentOut(i+seqNum, j) = inv_change(best(j, i));
    }
  }

#ifdef DEBUG_DIVERGE
  for(i = 0; i < alignmentOut.rows(); i++) {
    printf("%d\t", i);
    for(j = 0; j < alignmentOut.cols(); j++) {
      printf("%c", alignmentOut(i, j));
    }
    putchar('\n');
  }
  
  puts("done");
#endif
  
  return true;
}

/*----------------------------------------------------------------------*/

bool
parse_tree(const char *tree_str, int seqNum, int &treeE, IMatrix2D &tree) {
  const char *p = tree_str;
  int i;
  
  tree.resize(seqNum * 2, 2);
  tree = 0;

  i = seqNum + 1;
  int max = i;
  int left = 1;
  int right = 0;
  char aChar;
  int tmp;

  while(*p++ != '(') {}
  while((aChar = *p++) != '\r' && aChar != '\n' && aChar != '\0') {
    if(aChar != ' ' && aChar != ',') {
      if(aChar == '(') {
	left++;
	if(tree(i, 0) == 0) {
	  tree(i, 0) = i + 1;
	  i++;
	  max++;
	} else {
	  if(tree(i, 1) == 0) {
	    i = tree(i, 1) = ++max;
	  } else {
	    i = treeE = ++max;
	  }
	}
      } else {
	if(aChar == ')') {
	  int j;
	  right++;
	  for(j = seqNum + 1; j <= seqNum * 2 - 2; j++) {
	    if(tree(j, 0) == i || tree(j, 1) == i) {
	      break;
	    }
	  }
	  i = j;
	} else {
	  p--;
	  int n = 0;
	  if(tree(i, 0) == 0) {
	    tmp = sscanf(p, "%d%n", &tree(i, 0), &n);
	  } else {
	    if(tree(i, 1) == 0) {
	      tmp = sscanf(p, "%d%n", &tree(i, 1), &n);
	    } else {
	      tmp = sscanf(p, "%d%n", &treeE, &n);
	    }
	  }
	  p += n;
	}
      }
    }
  }

  if(left != right) {
    return false;
  }

  return true;
}

/*----------------------------------------------------------------------*/

static bool
common_compute(bool jttF, double alpha, const CMatrix2D &alignment,
	     int treeE, const IMatrix2D &tree,
	     const DVector &freq_, const DMatrix2D &prob_, BMatrix2D &seqN,
	     DVector &length, IMatrix2D &branches, IVector &parent,
	     IMatrix2D &best) {
  int i, j, k;
  int seqNum = alignment.rows();
  int seqLength = alignment.cols();

  DVector freq = freq_;
  DMatrix2D prob = prob_;
    
  /* COMPUTING GAMMA DISTANCES */

  DMatrix2D distance(seqNum - 1, seqNum);
  
  for(i = 0; i < seqNum - 1; i++) {
    for(j = i + 1; j < seqNum; j++) {
      distance(i, j) = 0.0;
      for(k = 0; k < seqLength; k++) {
	if(alignment(i, k) != alignment(j, k))
	  distance(i, j) += 1.0 / seqLength;
      }
      distance(i, j) = alpha * (pow(1.0 - distance(i, j), -1.0 / alpha) - 1.0);
    }
  }

  estimate_branch_lengths(seqNum, treeE, distance, tree, branches, parent, length);
  
#ifdef DEBUG_DIVERGE
  FILE *fp;
  fp=stdout;

  fprintf(fp, "tree\n");
  for(i = 0; i < 2 * seqNum; i++) {
	  fprintf(fp, "n = %d, left = %d, right = %d\n", i, tree(i, 0), tree(i, 1));
  }

  for(i = 0; i < 2 * seqNum - 3; i++) {
	  fprintf(fp, "site = %d, length = %f\n", i, length(i));
  }
#endif

  seqN.resize(seqNum, seqLength);
  for(i = 0; i < seqNum; i++) {
    for(j = 0; j < seqLength; j++) {
      seqN(i, j) = change(alignment(i, j));
    }
  }

  DVector freqObs(20);

  if(jttF) jttf(seqN, freqObs, freq, prob);

  DMatrix3D probability;
  compute_trans_matrix(seqNum, prob, probability, length);
  prob.clear();

  if(!reconstruct(treeE, seqN, probability,
		  branches, tree, freqObs, freq, best)) {
    return false;
  }

  //probability.clear();
  //freqObs.clear();
  //freq.clear();
  
  return true;
}

/*----------------------------------------------------------------------*/

static void
estimate_branch_lengths(int seqNum, int treeE, const DMatrix2D &distance,
			const IMatrix2D &tree, IMatrix2D &branches,
			IVector &parent, DVector &length) {
  int i, j, k, k1, k2;
  
  /* LEAST-SQUARE BRANCH LENGTH ESTIMATION */
  /* determine the four groups */

  branches.resize(seqNum * 2 - 3, 2);

  for(i = 0, j = seqNum + 1; i <= seqNum - 3; i++, j++) {
    branches(i, 0) = j;
    branches(i, 1) = tree(j, 0);
  }

  for(i = seqNum - 2, j = seqNum + 1; i <= seqNum * 2 - 5; i++, j++) {
    branches(i, 0) = j;
    branches(i, 1) = tree(j, 1);
  }

  branches(seqNum * 2 - 4, 0) = seqNum + 1;
  branches(seqNum * 2 - 4, 1) = treeE;

  parent.resize(seqNum * 2 - 1);
  
  for(i = 1; i <= seqNum * 2 - 2; i++) {
    for(j = seqNum + 1; j <= seqNum * 2 - 2; j++) {
      if((tree(j, 0) - i) * (tree(j, 1) - i) == 0)
	parent(i) = j;
    }
  }
  parent(treeE) = seqNum + 1;
  parent(seqNum + 1) = treeE;

  IMatrix3D group(seqNum * 2 - 3, 5, seqNum - 2);
  IMatrix2D noOfGroup(seqNum * 2 - 3, 5);
  DMatrix2D w(seqNum * 2 - 3, 6);

  group = 0;
  noOfGroup = 0;

  int min = 0, minimum;
  IVector x(5);
  for(i = 0; i < seqNum * 2 - 3; i++) {
    IVector box(seqNum + 1);
    for(j = 0; j < seqNum + 1; j++) {
      box(j) = j;
    }
    minimum = seqNum + 1;
    if(branches(i, 1) > seqNum) {
      x(3) = tree(branches(i, 1), 0);
      x(4) = tree(branches(i, 1), 1);
      x(1) = parent(branches(i, 0));
      if(branches(i, 0) == seqNum + 1)
		  x(1) = tree(seqNum + 1, 0);
      x(2) = tree(branches(i, 0), 1);
      if(x(2) == branches(i, 1))
		  x(2) = tree(branches(i, 0), 0);

      for(k = 1; k <= 4; k++) {
	if(x(k) > seqNum && minimum > abs(seqNum + 1 - x(k))) {
	  minimum = abs(seqNum + 1 - x(k));
	  min = k;
	}
      }
      for(k = 1; k <= 4; k++) {
	if(k != min) {
	  if(x(k) > seqNum) {
	    binaryGo(i, k, x(k), seqNum, group, noOfGroup, tree);
	  } else {
	    noOfGroup(i, k) = 1;
	    group(i, k, 0) = x(k);
	  }
	}
      }
      noOfGroup(i, min) = seqNum - noOfGroup(i, 1) - noOfGroup(i, 2) - noOfGroup(i, 3) - noOfGroup(i, 4);


      for(j = 1; j <= 4; j++) {
	if(j != min) {
	  for(k = 0; k < noOfGroup(i, j); k++) {
	    box(group(i, j, k)) = -1;
	  }
	}
      }
      j = 0;
      for(k = 1; k <= seqNum; k++) {
	if(box(k) != -1) {
	  group(i, min, j++) = box(k);
	}
      }

      w(i, 0) = -.5 / noOfGroup(i, 1) / noOfGroup(i, 2);
      w(i, 1) = 1.0 * (noOfGroup(i, 2) * noOfGroup(i, 3) + noOfGroup(i, 1) * noOfGroup(i, 4) + .0) / ((noOfGroup(i, 1) + noOfGroup(i, 2)) * (noOfGroup(i, 3) + noOfGroup(i, 4)) * (2.0 * noOfGroup(i, 1) * noOfGroup(i, 3)) + .0);
      w(i, 2) = 1.0 * (noOfGroup(i, 1) * noOfGroup(i, 3) + noOfGroup(i, 2) * noOfGroup(i, 4) + .0) / ((noOfGroup(i, 1) + noOfGroup(i, 2)) * (noOfGroup(i, 3) + noOfGroup(i, 4)) * (2.0 * noOfGroup(i, 1) * noOfGroup(i, 4)) + .0);
      w(i, 3) = 1.0 * (noOfGroup(i, 1) * noOfGroup(i, 3) + noOfGroup(i, 2) * noOfGroup(i, 4) + .0) / ((noOfGroup(i, 1) + noOfGroup(i, 2)) * (noOfGroup(i, 3) + noOfGroup(i, 4)) * (2.0 * noOfGroup(i, 2) * noOfGroup(i, 3)) + .0);
      w(i, 4) = 1.0 * (noOfGroup(i, 2) * noOfGroup(i, 3) + noOfGroup(i, 1) * noOfGroup(i, 4) + .0) / ((noOfGroup(i, 1) + noOfGroup(i, 2)) * (noOfGroup(i, 3) + noOfGroup(i, 4)) * (2.0 * noOfGroup(i, 2) * noOfGroup(i, 4)) + .0);
      w(i, 5) = -.5 / (noOfGroup(i, 3) * noOfGroup(i, 4) + .0);
    } else {
      noOfGroup(i, 4) = -1;
      noOfGroup(i, 3) = 1;
      group(i, 3, 0) = branches(i, 1);
      x(1) = parent(branches(i, 0));
      if(branches(i, 0) == seqNum + 1)
	x(1) = tree(seqNum + 1, 0);
      x(2) = tree(branches(i, 0), 1);
      if(x(2) == branches(i, 1))
	x(2) = tree(branches(i, 0), 0);

      for(k = 1; k <= 4; k++) {
	if(x(k) > seqNum && minimum > abs(seqNum + 1 - x(k))) {
	  minimum = abs(seqNum + 1 - x(k));
	  min = k;
	  break;
	}
      }

      for(k = 1; k <= 2; k++) {
	if(k != min) {
	  if(x(k) > seqNum) {
	    binaryGo(i, k, x(k), seqNum, group, noOfGroup, tree);
	  } else {
	    noOfGroup(i, k) = 1;
	    group(i, k, 0) = x(k);
	  }
	}
      }
      noOfGroup(i, min) = seqNum - noOfGroup(i, 1) - noOfGroup(i, 2) - noOfGroup(i, 3);
      for(j = 1; j <= 3; j++) {
	if(j != min) {
	  for(k = 0; k < noOfGroup(i, j); k++)
	    box(group(i, j, k)) = -1;
	}
      }
      for(k = 1, j = 0; k <= seqNum; k++) {
	if(box(k) != -1) {
	  group(i, min, j++) = box(k);
	}
      }

      w(i, 0) = .5 / noOfGroup(i, 1);
      w(i, 1) = .5 / noOfGroup(i, 2);
      w(i, 2) = -.5 / (noOfGroup(i, 1) * noOfGroup(i, 2));
    }
  }

  /* branch lengths estimation */

  length.resize(seqNum * 2 - 3);
  length = 0;
  
  for(i = 0; i < seqNum * 2 - 3; i++) {
    int left = 0, right = 0;
	
    if(branches(i, 1) > seqNum) {
      for(j = 1; j <= seqNum - 1; j++) {
	for(k1 = 1; k1 <= 4; k1++) {
	  for(k2 = 0; k2 < noOfGroup(i, k1); k2++) {
	    if(j == group(i, k1, k2))
	      left = k1;
	  }
	}

	for(k = j + 1; k <= seqNum; k++) {
	  for(k1 = 1; k1 <= 4; k1++) {
	    for(k2 = 0; k2 < noOfGroup(i, k1); k2++) {
	      if(k == group(i, k1, k2))
		right = k1;
	    }
	  }
	  if((left == 1 && right == 2) || (left == 2 && right == 1))
	    length(i) += distance(j - 1, k - 1) * w(i, 0);
	  else if((left == 1 && right == 3) || (left == 3 && right == 1))
	    length(i) += distance(j - 1, k - 1) * w(i, 1);
	  else if((left == 1 && right == 4) || (left == 4 && right == 1))
	    length(i) += distance(j - 1, k - 1) * w(i, 2);
	  else if((left == 2 && right == 3) || (left == 3 && right == 2))
	    length(i) += distance(j - 1, k - 1) * w(i, 3);
	  else if((left == 2 && right == 4) || (left == 4 && right == 2))
	    length(i) += distance(j - 1, k - 1) * w(i, 4);
	  else if((left == 3 && right == 4) || (left == 4 && right == 3))
	    length(i) += distance(j - 1, k - 1) * w(i, 5);
	}
      }
    } else {
      for(j = 1; j <= seqNum - 1; j++) {
	for(k1 = 1; k1 <= 3; k1++) {
	  for(k2 = 0; k2 < noOfGroup(i, k1); k2++) {
	    if(j == group(i, k1, k2))
	      left = k1;
	  }
	}

	for(k = j + 1; k <= seqNum; k++) {
	  for(k1 = 1; k1 <= 3; k1++) {
	    for(k2 = 0; k2 < noOfGroup(i, k1); k2++) {
	      if(k == group(i, k1, k2))
		right = k1;
	    }
	  }
	  if((left == 1 && right == 3) || (left == 3 && right == 1))
	    length(i) += distance(j - 1,k - 1) * w(i, 0);
	  else if((left == 2 && right == 3) || (left == 3 && right == 2))
	    length(i) += distance(j - 1, k - 1) * w(i, 1);
	  else if((left == 1 && right == 2) || (left == 2 && right == 1))
	    length(i) += distance(j - 1, k - 1) * w(i, 2);
	}
      }
    }
  }
}

/*----------------------------------------------------------------------*/

static void
compute_trans_matrix(int seqNum, const DMatrix2D &prob, DMatrix3D &probability, DVector &length) {
  int i, j, k;

  /* COMPUTE TRANSITION MATRIX FOR EACH BRANCH */

  for(i = 0; i < seqNum * 2 - 3; i++) {
    if(length(i) < 0.0)
      length(i) = 0.0;
  }

  IVector order(seqNum * 2 - 3);
  order(0) = 0;
  for(i = 1; i < seqNum * 2 - 3; i++) {
    if(length(i) < length(order(i - 1))) {
      if(length(i) < length(order(0))) {
	for(k = i; k > 0; k--) {
	  order(k) = order(k - 1);
	}
	order(0) = i;
      } else {
	for(j = i - 1; j >= 0; j--) {
	  if(length(i) >= length(order(j))) {
	    for(k = i; k > j + 1; k--) {
	      order(k) = order(k - 1);
	    }
	    order(j + 1) = i;
	    break;
	  }
	}
      }
    } else {
      order(i) = i;
    }
  }

  probability.resize(2 * seqNum - 3, 20, 20);

  DMatrix2D transition(20,20);
  for(i = 0; i < 20; i++) {
    for(j = 0; j < 20; j++) {
      if(i == j) {
	transition(i,j) = 1.0;
      } else {
	transition(i,j) = 0.0;
      }
    }
  }


#ifdef DEBUG_DIVERGE
  for(i = 0; i < (2 * seqNum - 3); i++) {
    printf("order(%d) = %d\n", i, order(i));
  }
#endif

  for(i = 0; length(order(i)) * 100.0 < 1.0 && (i < 2 * (seqNum - 3)); i++) {
#ifdef DEBUG_DIVERGE
    printf("Enter loop\n");
#endif

    for(j = 0; j < 20; j++) {
      for(k = 0; k < 20; k++) {
	if(j != k) {
	  probability(order(i), j, k) = prob(j,k) * length(order(i)) * 100;
	} else {
	  double v = 1.0 - length(order(i)) * 100 * (1.0 - prob(j,k));
#ifdef DEBUG_DIVERGE
	  printf("i = %d, j = %d, k = %d, order(i) = %d, length(order(i)) = %f, prob(j, k) = %f, probability(order(i), j, k) = %f, v = %f\n", 
	       i, j, k, order(i), length(order(i)), prob(j, k), probability(order(i), j, k), v);
#endif
	  
	  probability(order(i), j, k) = v;
	}
      }
    }
  }

  // Do not change i... the next block depends on it
  for(j = 1; j <= 100 * length(order(2 * seqNum - 4)); j++) {
    matrix(prob, transition);

    do {
      double l = 100 * length(order(i));
      if((int)l != j)
	break;

      for(int k1 = 0; k1 < 20; k1++) {
	for(int k2 = 0; k2 < 20; k2++) {
	  double p = 0.0;
	  if(k1 != k2) {
	    p = transition(k1,k2) + prob(k1,k2) * (l - j);
	  } else {
	    p = transition(k1,k2) - (1.0 - prob(k1,k2)) * (l - j);
	  }
	  probability(order(i), k1, k2) = p;
	}
      }
      i++;
    } while(i <= seqNum * 2 - 4);
  }
}

/*----------------------------------------------------------------------*/

static bool
reconstruct(int treeE,
	    const BMatrix2D &seqN,
	    const DMatrix3D &probability,
	    const IMatrix2D &branches,
	    const IMatrix2D &tree,
	    const DVector &freqObs, const DVector &freq,
	    IMatrix2D &best) {
  /* RECONSTRUTING ANCESTRAL STATES */

  int i, j, k;
  
  int seqNum = seqN.rows();
  int seqLength = seqN.cols();

  best.resize(seqLength, seqNum);

  for(i = 0; i < seqLength; i++) {
#ifdef DEBUG_DIVERGE
    printf("%d / %d\n", i+1, seqLength);
#endif
    double probPath = totalProb(i, treeE, seqN, probability, branches, tree, freq);
    for(j = seqNum + 1; j <= seqNum * 2 - 2; j++) {
      double max1 = .0;
      int max = -1;
      for(k = 0; k < 20; k++) {
	double iiii = nodeByNode(i, j, k, treeE, seqN, probability, branches, tree, freqObs);
	double p =  iiii / probPath;
	if(p > max1) {
	  max1 = p;
	  max = k;
	}
      }
      best(i, j - seqNum) = max;
    }
  }
  
  return true;
}

/*----------------------------------------------------------------------*/

/* DETERMINE FOUR GROUPS */

static void
binaryGo(int branchI, int groupI, int j, int seqNum,
	 IMatrix3D &group, IMatrix2D &noOfGroup,
	 const IMatrix2D &tree) {
  if(tree(j, 0) > seqNum) {
    binaryGo(branchI, groupI, tree(j, 0), seqNum, group, noOfGroup, tree);
  }
  if(tree(j, 0) < seqNum + 1) {
    group(branchI, groupI, noOfGroup(branchI, groupI)) = tree(j, 0);
    noOfGroup(branchI, groupI)++;
  }
  
  if(tree(j, 1) > seqNum) {
    binaryGo(branchI, groupI, tree(j, 1), seqNum, group, noOfGroup, tree);
  }
  if(tree(j, 1) < seqNum + 1) {
    group(branchI, groupI, noOfGroup(branchI, groupI)) = tree(j, 1);
    noOfGroup(branchI, groupI)++;
  }
}

/*----------------------------------------------------------------------*/

/* MATRIX MULTIPLICATION */

static void
matrix(const DMatrix2D &prob, DMatrix2D &transition) {
  int i, j, k;
  DMatrix2D temp(20,20);
  for(i = 0; i < 20; i++) {
    for(j = 0; j < 20; j++) {
      temp(i,j) = .0;
      for(k = 0; k < 20; k++) {
	temp(i,j) += transition(i,k) * prob(k,j);
      }
    }
  }
  transition = temp;
}

/*----------------------------------------------------------------------*/

/* CHANGE from amino acid to number */

static int
change(char amnio) {
  int rv = -1;
  
  if(amnio == 'A') rv = 0;
  else if(amnio == 'R') rv = 1;
  else if(amnio == 'N') rv = 2;
  else if(amnio == 'D') rv = 3;
  else if(amnio == 'C') rv = 4;
  else if(amnio == 'Q') rv = 5;
  else if(amnio == 'E') rv = 6;
  else if(amnio == 'G') rv = 7;
  else if(amnio == 'H') rv = 8;
  else if(amnio == 'I') rv = 9;
  else if(amnio == 'L') rv = 10;
  else if(amnio == 'K') rv = 11;
  else if(amnio == 'M') rv = 12;
  else if(amnio == 'F') rv = 13;
  else if(amnio == 'P') rv = 14;
  else if(amnio == 'S') rv = 15;
  else if(amnio == 'T') rv = 16;
  else if(amnio == 'W') rv = 17;
  else if(amnio == 'Y') rv = 18;
  else if(amnio == 'V') rv = 19;

  return rv;
}

/*----------------------------------------------------------------------*/

/* CHANGE from number to amino acid */

static char
inv_change(int amnio) {
  int rv = -1;
  
  if(     amnio ==  0) rv = 'A';
  else if(amnio ==  1) rv = 'R';
  else if(amnio ==  2) rv = 'N';
  else if(amnio ==  3) rv = 'D';
  else if(amnio ==  4) rv = 'C';
  else if(amnio ==  5) rv = 'Q';
  else if(amnio ==  6) rv = 'E';
  else if(amnio ==  7) rv = 'G';
  else if(amnio ==  8) rv = 'H';
  else if(amnio ==  9) rv = 'I';
  else if(amnio == 10) rv = 'L';
  else if(amnio == 11) rv = 'K';
  else if(amnio == 12) rv = 'M';
  else if(amnio == 13) rv = 'F';
  else if(amnio == 14) rv = 'P';
  else if(amnio == 15) rv = 'S';
  else if(amnio == 16) rv = 'T';
  else if(amnio == 17) rv = 'W';
  else if(amnio == 18) rv = 'Y';
  else if(amnio == 19) rv = 'V';

  return rv;
}

/*----------------------------------------------------------------------*/

static double
totalProb(int site, int treeE, const BMatrix2D &seqN, const DMatrix3D &probability,
	  const IMatrix2D &branches, const IMatrix2D &tree, const DVector &freq) {
  int i, j, k, k1 = -1, k2 = -1;
  double sum1, sum2, sum3;

  int seqNum = seqN.rows();
  
  DMatrix2D likelihood(seqNum * 2 - 1, 21);
  IVector nfailures(seqNum - 3);
  nfailures = 0;
  
  for(i = 1; i < seqNum + 1; i++) {
    for(j = 0; j < 20; j++) {
      likelihood(i, j) = 0.0;
      if(j == seqN(i - 1, site))
	likelihood(i, j) = 1.0;
    }
    likelihood(i, 20) = 1.0;
  }

  for(i = seqNum + 1; i < seqNum * 2 - 1; i++) {
    likelihood(i, 20) = -1.0;
  }

out600:
  for(i = seqNum + 2; i < seqNum * 2 - 1; i++) {
    if((likelihood(tree(i, 0), 20) > .0) && (likelihood(tree(i, 1), 20) > .0) && (likelihood(i, 20) < .0)) {
      likelihood(i, 20) = 1.0;
      for(j = 0; j < seqNum * 2 - 3; j++) {
	if(branches(j, 0) == i && branches(j, 1) == tree(i, 0))
	  k1 = j;
	if(branches(j, 0) == i && branches(j, 1) == tree(i, 1))
	  k2 = j;
      }
      for(j = 0; j < 20; j++) {
	sum1 = sum2 = 0.0;
	for(k = 0; k < 20; k++) {
	  sum1 += likelihood(tree(i, 0), k) * probability(k1, j, k);
	  sum2 += likelihood(tree(i, 1), k) * probability(k2, j, k);
	}
	likelihood(i, j) = sum1 * sum2;
      }
    }
  }
  for(i = seqNum + 2; i < seqNum * 2 - 1; i++) {
    if(likelihood(i, 20) < 0.0) {
      //if(nfailures(i - (seqNum + 2)) < 10) {
	//nfailures(i - (seqNum + 2))++;
	goto out600;
      //}
    }
  }

  for(j = 0; j < seqNum * 2 - 3; j++) {
    if(branches(j, 0) == seqNum + 1 && branches(j, 1) == tree(seqNum + 1, 1))
      k1 = j;
  }

  double sum = 0.0;
  for(j = 0; j < 20; j++) {
    sum1 = sum2 = sum3 = 0.0;
    for(k = 0; k < 20; k++) {
      sum1 += likelihood(tree(seqNum + 1, 0), k) * probability(0, j, k);
      sum2 += likelihood(tree(seqNum + 1, 1), k) * probability(k1, j, k);
      sum3 += likelihood(treeE, k) * probability(seqNum * 2 - 4, j, k);
    }
    likelihood(seqNum + 1, j) = sum1 * sum2 * sum3 * freq(j);
    sum += likelihood(seqNum + 1, j);
  }

  return sum;
}

/*----------------------------------------------------------------------*/

static void
jttf(const BMatrix2D &seqN, DVector &freqObs, DVector &freq, DMatrix2D &prob) {
  int i, j;
  double aDouble;

  int seqNum = seqN.rows();
  int seqLength = seqN.cols();
  
  aDouble = 1.0 / (seqLength * seqNum);
  freqObs = 0;
  for(i = 0; i < seqNum; i++) {
    for(j = 0; j < seqLength; j++) {
      freqObs(seqN(i, j)) += aDouble;
    }
  }
  for(i = 0; i < 20; i++) {
    prob(i,i) = 1.0;
    for(j = 0; j < 20; j++) {
      if(i != j) {
	prob(i,j) *= freqObs(j) / freq(j);
	prob(i,i) -= prob(i,j);
      }
    }
  }
  freq = freqObs;
}

/*----------------------------------------------------------------------*/

static void
stepDistance(int treeE, const IMatrix2D &best, const BMatrix2D &seqN,
	     const DVector &length, const IMatrix2D &branches,
	     DVector &distStep1) {
  int i, site, x, y;

  double B = 0.0;
  int seqNum = seqN.rows();
  int seqLength = seqN.cols();
  int M = 2 * seqNum - 3;
  
  for(i = 0; i < 2 * seqNum - 3; i++) {
    B += length(i);
  }
  
  IVector set(2 * seqNum - 3);
  for(site = 0; site < seqLength; site++) {
    int sum = 0;
    for(i = 0; i < 2 * seqNum - 3; i++)
      set(i) = 2 * seqNum;
    for(i = 0; i < seqNum * 2 - 4; i++) {
      if(branches(i, 0) > seqNum) {
	x = best(site, branches(i, 0) - seqNum);
      } else {
	x = seqN(branches(i, 0) - 1, site);
      }
      if(branches(i, 1) > seqNum) {
	y = best(site, branches(i, 1) - seqNum);
      } else {
	y = seqN(branches(i, 1) - 1, site);
      }
      if(x != y) {
	set(sum++) = i;
      }
    }
    x = best(site, 1);
    if(treeE > seqNum)
      y = best(site, treeE - seqNum);
    else {
      y = seqN(treeE - 1, site);
    }
    if(x != y) {
      set(sum++) = 2 * seqNum - 4;
    }
    if(sum == 0) {
      distStep1(site) = 0.0;
    } else {
      int m = sum;
      double a = (m + .0) / B;
      double b = (M - 1.0) / B;
      double c = (a + b) / 2.0;
      double c0 = c;
      double oka = B;
      double okb = B;
      for(i = 0; i < m; i++) {
	if(length(set(i)) != .0) {
	  oka -= length(set(i)) / (1.0 - exp(0.0 - length(set(i)) * a));
	  okb -= length(set(i)) / (1.0 - exp(0.0 - length(set(i)) * b));
	} else {
	  oka -= 1.0 / a;
	  okb -= 1.0 / b;
	}
      }

      do {
	double okc = B;
	for(i = 0; i < m; i++) {
	  if(length(set(i)) != .0) {
	    okc -= length(set(i)) / (1.0 - exp(0.0 - (length(set(i)) + .0000001) * c));
	  } else {
	    okc -= 1.0 / c;
	  }
	}
	if(okc > .0) {
	  c0 = c;
	  b = c;
	  c = (c + a) / 2.0;
	}
	if(okc < .0) {
	  c0 = c;
	  a = c;
	  c = (c + b) / 2.0;
	}
      } while((c0 - c) > 0.001 || (c - c0) > 0.001);
      distStep1(site) = c * B;
    }
  }
}

/*----------------------------------------------------------------------*/

static void
gammaDistanceNew1(const BMatrix2D &seqN, const DVector &length,
		  const DVector &distStep1) {
  int i;
  double mean, variance, moment;

  double alpha = 0.0;
  double sum = 0.0;
  double sumSquare = 0.0;

  int seqLength = seqN.cols();
  
  for(i = 0; i < seqLength; i++) {
    sum += distStep1(i);
    sumSquare += distStep1(i) * distStep1(i);
  }
  mean = sum / (seqLength + .0);
  variance = sumSquare / (seqLength - 1.0) - sum * sum / (seqLength - 1.0) / (seqLength + .0);
  moment = mean * mean / (variance - mean);
  if(moment > .0) {
    alpha = sullivan(mean, moment, distStep1);
#ifdef DEBUG_DIVERGE
    printf("=== gammaDistanceNew, alpha = %f\n", alpha);
#endif
  } else {
    // Unable to estimate alpha;
  }
  gz97_inferred_alpha = alpha;
  rateF(alpha, seqN, length, distStep1);
}


/*----------------------------------------------------------------------*/

static double
gammln(double xx) {
  double x, tmp, ser;
  static double cof[6] = { 76.18009173, -86.50532033,     24.01409822,
			   -1.231739516,  0.120868003e-2, -0.536382e-5 };
  int j;

  x = xx - 1.0;
  tmp = x + 5.5;
  tmp -= (x + 0.5) * log(tmp);
  ser = 1.0;
  for(j = 0; j <= 5; j++) {
    x += 1.0;
    ser += cof[j] / x;
  }

  return -tmp + log(2.50662827465 * ser);
}

/*----------------------------------------------------------------------*/

static double
sullivan(double mean, double z, const DVector &distStep) {
  int k;
  double a, x, y, fx, fy, fz, fa;

  int seqLength = distStep.size();
  
  z *= 2;
  x = .001;
  y = z - .01;
  fx = .0;
  fy = .0;
  fz = .0;
  
  for(k = 0; k < seqLength; k++) {
    fx += gammln(distStep(k) + x) - gammln(x) - gammln(distStep(k) + 1) - x * log(1.0 + mean / x) + distStep(k) * log(mean / (mean + x));
    fy += gammln(distStep(k) + y) - gammln(y) - gammln(distStep(k) + 1) - y * log(1.0 + mean / y) + distStep(k) * log(mean / (mean + y));
    fz += gammln(distStep(k) + z) - gammln(z) - gammln(distStep(k) + 1) - z * log(1.0 + mean / z) + distStep(k) * log(mean / (mean + z));
  }

  do {
    a = (x + y) / 2.0;
    fa = .0;
    for(k = 0; k < seqLength; k++) {
      fa += gammln(distStep(k) + a) - gammln(a) - gammln(distStep(k) + 1) - a * log(1.0 + mean / a) + distStep(k) * log(mean / (mean + a));
    }
    if(fa < fy) {
      x = a;
      fx = fa;
    } else {
      z = y;
      fz = fy;
      y = a;
      fy = fa;
    }

    a = (y + z) / 2.0;
    fa = .0;
    for(k = 0; k < seqLength; k++) {
      fa += gammln(distStep(k) + a) - gammln(a) - gammln(distStep(k) + 1) - a * log(1.0 + mean / a) + distStep(k) * log(mean / (mean + a));
    }
    if(fa < fy) {
      z = a;
      fz = fa;
    } else {
      x = y;
      fx = fy;
      y = a;
      fy = fa;
    }
  } while((z - x) > .001);

  return y;
}

/*----------------------------------------------------------------------*/

static void
rateF(double alpha, const BMatrix2D &seqN, const DVector &length,
      const DVector &distStep1) {
  int i, j;
  double mutability[20] = { 1.0, .83, 1.04, .86, .44, .84, .77, .50, .91, 1.03,
			    .54, .72, .93, .51, .58, 1.17, 1.07, .25, .5, .98 };
  double B, sum, sum1;

  int seqNum = seqN.rows();
  int seqLength = seqN.cols();
  
  DVector xx(seqLength), xxx(seqLength), muta(seqLength);

  B = .0;
  sum = .0;
  for(i = 0; i < 2 * seqNum - 3; i++)
    B += length(i);
  for(i = 0; i < seqLength; i++) {
    xx(i) = (alpha + distStep1(i) - 1) / (alpha + B);
    if(xx(i) < .0)
      xx(i) = .0;
    sum += xx(i);
  }
  for(i = 0; i < seqLength; i++) {
    xx(i) *= (seqLength + .0) / sum;
  }

  sum = 0.0;
  sum1 = 0.0;
  for(i = 0; i < seqLength; i++) {
    muta(i) = 0.0;
    for(j = 0; j < seqNum; j++) {
      muta(i) += (mutability[seqN(j, i)] + .0) / (seqNum + .0);
    }
    xxx(i) = xx(i) / muta(i);
    sum += muta(i) / (seqLength + .0);
    sum1 += xxx(i) / (seqLength + .0);
  }
}

/*----------------------------------------------------------------------*/

static double
nodeByNode(int site, int node, int xxx, int treeE, const BMatrix2D &seqN,
	   const DMatrix3D &probability, const IMatrix2D &branches,
	   const IMatrix2D &tree, const DVector &freq) {
  int i, j, k, k1 = -1, k2 = -1;
  double sum1, sum2, sum3;

  int seqNum = seqN.rows();
  
  DMatrix2D likelihood(seqNum * 2 - 1, 21);
  IVector nfailures(seqNum - 3);
  nfailures = 0;

  for(i = 1; i < seqNum + 1; i++) {
    for(j = 0; j < 20; j++) {
      likelihood(i, j) = 0.0;
      if(j == seqN(i - 1, site))
	likelihood(i, j) = 1.0;
    }
    likelihood(i, 20) = 1.0;
  }

  //for(i = seqNum + 1; i < seqNum * 2 - 1; i++) {
  for(i = seqNum + 1; i < seqNum * 2 - 1; i++) {
    likelihood(i, 20) = -1.0;
  }

out600:
  for(i = seqNum + 2; i < seqNum * 2 - 1; i++) {
	  int myTemp1 = tree(i, 0); 
	  int myTemp2 = tree(i, 1);
    if((likelihood(myTemp1, 20) > .0) && (likelihood(myTemp2, 20) > .0) && (likelihood(i, 20) < .0)) {
      likelihood(i, 20) = 1.0;
      for(j = 0; j < seqNum * 2 - 3; j++) {
	if(branches(j, 0) == i && branches(j, 1) == tree(i, 0))
	  k1 = j;
	if(branches(j, 0) == i && branches(j, 1) == tree(i, 1))
	  k2 = j;
      }
      for(j = 0; j < 20; j++) {
	sum1 = sum2 = 0.0;
	for(k = 0; k < 20; k++) {
	  sum1 += likelihood(tree(i, 0), k) * probability(k1, j, k);
	  sum2 += likelihood(tree(i, 1), k) * probability(k2, j, k);
	}
	likelihood(i, j) = sum1 * sum2;

      }
      if(node == i) {
	for(j = 0; j < 20; j++) {
	  if(j != xxx)
	    likelihood(i, j) = 0.0;
	}
      }

    }
  }

  for(i = seqNum + 2; i < seqNum * 2 - 1; i++) {
    if(likelihood(i, 20) < 0.0) {
    //if(nfailures(i - (seqNum + 2)) < 10) {
	//nfailures(i - (seqNum + 2))++;
	  goto out600;
    //}
    }
  }

  for(j = 0; j < seqNum * 2 - 3; j++) {
    if(branches(j, 0) == seqNum + 1 && branches(j, 1) == tree(seqNum + 1, 1))
      k1 = j;
  }

	/////////////////////From Here//////////////
  if((seqNum+1)!=node) {
   for(j=0;j<20;j++) {
	   sum1=sum2=sum3=.0;
	   for(k = 0; k < 20; k++) {
		  sum1 += likelihood(tree(seqNum + 1, 0), k) * probability(0, j, k);
		  double mytest12 = probability(k1, j, k);
		  double mytest13 = likelihood(tree(seqNum + 1, 1), k);
		  sum2 += likelihood(tree(seqNum + 1, 1), k) * probability(k1, j, k);

		  sum3 += likelihood(treeE, k) * probability(seqNum * 2 - 4, j, k);
      }
	   likelihood(seqNum + 1, j) = sum1 * sum2 * sum3 * freq(j);
   }
  }
  else { 
	  for(j=0;j<20;j++) 
		  likelihood(seqNum+1, j)=0.0;
	  j=xxx;
	  sum1=sum2=sum3=0.0;

	  for(k = 0; k < 20; k++) {
		  sum1 += likelihood(tree(seqNum + 1, 0), k) * probability(0, j, k);
		  sum2 += likelihood(tree(seqNum + 1, 1), k) * probability(k1, j, k);
		  sum3 += likelihood(treeE, k) * probability(seqNum * 2 - 4, j, k);
      }
	  likelihood(seqNum + 1, j) = sum1 * sum2 * sum3 * freq(j);
 }

 double sum=0.0;
 for(j=0;j<20;j++) {
   sum=sum+likelihood(seqNum+1, j);
 }
  
 return(sum);

 //////////////END///////////////

 /*
  double sum = 0.0;
  for(j = 0; j < 20; j++) {
    likelihood(seqNum + 1, j) = 0.0;
    if((seqNum + 1) != node || j == xxx) {
      sum1 = sum2 = sum3 = 0.0;

      for(k = 0; k < 20; k++) {
	sum1 += likelihood(tree(seqNum + 1, 0), k) * probability(0, j, k);
	sum2 += likelihood(tree(seqNum + 1, 1), k) * probability(k1, j, k);
	sum3 += likelihood(treeE, k) * probability(seqNum * 2 - 4, j, k);
      }
      likelihood(seqNum + 1, j) = sum1 * sum2 * sum3 * freq(j);
      sum += likelihood(seqNum + 1, j);
    }
  }
  
  return sum;
  */
}

/*----------------------------------------------------------------------*/
