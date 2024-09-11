#ifndef _GZ97_H_
#define _GZ97_H_


#include <string>
#include <vector>

#include "tree.h"
#include "sequence.h"

bool gu2001_compute(const std::vector<Tree> &trees,
		  const std::vector<sequence_t> &sequences,
		  std::vector<std::vector<double> > &summary,
		  std::vector<std::vector<double> > &rets2);

void Subtree2_cal(double e0[], double e1[], double e2[], int n, double alpha);
void OutputSubtree2(double mle[], double post_profile[], int n);
void PostProfile(double e0[], double e1[], double e2[], double theta, int n, double post_profile[]);
double theta2(double e0[], double e1[], double e2[], int n);
double SE_theta2(double e0[], double e1[], double e2[], double theta, int n);
double loglik_value(double e0[], double e1[],double e2[], double theta, int n);
bool ExpMarkov2(char X1[], int n1, int Left1[], int Right1[],
	 double Blength1[], char X2[], int n2,int Left2[], int Right2[],
	 double Blength2[], double logExp[3], double alpha, double freq[20]);
void gammad(double rate[], double alpha,int K);
double igamma(double z, double alpha);
double gammafun(double x);
bool MarkovProb(char X[], int n, double rate,
		  int Left[], int Right[], double Blength[],
		  double freq[20], double & prob_result);
double tranPr(int i, int j, double v);

#endif
 