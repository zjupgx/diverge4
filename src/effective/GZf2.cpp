/*This is the program for computing the coefficient of functional divergence
under the two-state model*/
#define _CRT_SECURE_NO_WARNINGS
#define WIN32
#include <vector>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "GZf2.h"
#include "matrix.h"


/*------------------------------------------------------------------------*/

using namespace std;

static double covGZ(const DVector &X, const DVector &Y, int n, int index);
static double covGZ(const DVector &X, const DVector &Y, int n, int index,
		    double &theta, double &se, double &r_X, double &r_max,
		    double &z_score);
static void likGZ(const DVector &X, const DVector &Y, int n, double &thetaML, double &alphaML, double &se_theta, double &LRT);
static double lik_cal(const DVector &X, const DVector &Y, int n, double alpha, double theta);
static double K_cal(double X, double Y, double alpha, double D1, double D2);
static double Q_cal(double X, double alpha, double D);
static double average(const DVector &A, int n);
static double var(const DVector &X, int n);
static double covariance(const DVector &X, const DVector &Y, int n);
static void postRatio(const DVector &X, const DVector &Y, int n, DVector &Z, double thetaML, double alphaML);
static double gammafun(double x);

/*------------------------------------------------------------------------*/

bool
GZf2_compute(const vector<DVector> &numSub, vector<vector<double> > &summary, vector<vector<double> > &rets) {
  int groupNum = numSub.size();
  int siteNum = numSub[0].size();
  
  DVector Z(siteNum);
  for(int i = 0; i < groupNum; i++) {
    for(int j = i + 1; j < groupNum; j++) {
      const DVector &X = numSub[i];
      const DVector &Y = numSub[j];

      vector<double> sum(9);
      // sum[0] = theta,    sum[1] = se
      // sum[2] = r_X,      sum[3] = r_max,
      // sum[4] = z_score   (sum[0-4] are model-free
      // sum[5] = thetaML,  sum[6] = alphaML,
      // sum[7] = se_theta, sum[8] = LRT
      
      //double theta, se, r_X, r_max, z_score;
      covGZ(X, Y, siteNum, 0, sum[0], sum[1], sum[2], sum[3], sum[4]);

      likGZ(X, Y, siteNum, sum[5], sum[6], sum[7], sum[8]);
      summary.push_back(sum);
      postRatio(X, Y, siteNum, Z, sum[5], sum[6]);


      int n=Z.size();
      vector<double> Z2(Z.size());
      for(int k=0; k<n; k++) {
	Z2[k] = Z(k);
      }
      rets.push_back(Z2);
      /*Gsite(X, Y, siteNum, Z);*/
    }
  }

  return true;
}

/*------------------------------------------------------------------------*/

static double
covGZ(const DVector &X, const DVector &Y, int n, int index){
  double theta, se, r_X, r_max, z_score;
  return covGZ(X, Y, n, index, theta, se, r_X, r_max, z_score);
}

/*------------------------------------------------------------------------*/

static double
covGZ(const DVector &X, const DVector &Y, int n, int index,
      double &theta, double &se, double &r_X, double &r_max, double &z_score) {
  double alpha;
  double xM, yM, xV, yV, cov;
  double z_X, z_max;

  xM = average(X, n);
  yM = average(Y, n);
  xV = var(X, n);
  yV = var(Y, n);
  cov = covariance(X, Y, n);

  theta = 1.0 - cov / sqrt((xV - xM) * (yV - yM));
  alpha = (xM * xM + yM * yM) / (xV + yV - xM - yM);

  /*statistical testing */
  r_X = cov / sqrt(xV * yV);
  r_max = sqrt((1.0 - xM / xV) * (1.0 - yM / yV));
  z_X = 0.5 * log((1.0 + r_X) / (1.0 - r_X));
  z_max = 0.5 * log((1.0 + r_max) / (1.0 - r_max));
  z_score = (z_X - z_max) * sqrt(n - 3.0);

  se = (1.0 - r_X * r_X) / r_max / sqrt(n - 3.0);

  if(index == 0)
    return 0;
  else if(index == 1)
    return alpha;
  else if(index == 2)
    return theta;
  else
    return -1;
}

/*------------------------------------------------------------------------*/

static void
likGZ(const DVector &X, const DVector &Y, int n, double &thetaML, double &alphaML, double &se_theta, double &LRT) {
  double lik, maxLik, lik0, c;
  double alpha, theta;
  double a_left, a_right, b_left, b_right;
  double a_max, b_max, step, a, b;
  int i, j, ii;

  alpha = covGZ(X, Y, n, 1);
  theta = covGZ(X, Y, n, 2);

  maxLik = lik_cal(X, Y, n, alpha, theta);

  a_left = 0.0;			/* for rho=1/(1+alpha) */
  b_left = 0.0;			/* for theta */
  a_right = 1.0;
  b_right = 1.0;
  a_max = 1.0 / (1.0 + alpha);
  b_max = theta;
  step = 0.1;
  ii = 0;

  do {
    ii++;
    step = (a_right - a_left) / 10.0;
    for(i = 0; i <= 10; i++) {
      for(j = 0; j <= 10; j++) {
	a = a_left + step * i;
	b = b_left + step * j;

	if(a == 0.0)
	  alpha = 20.0;
	else if(a >= 1.0)
	  alpha = 0.01;
	else
	  alpha = 1.0 / a - 1.0;

	if(b >= 1.0)
	  theta = 0.99;
	else if(b <= 0.0)
	  theta = 0.01;
	else
	  theta = b;

	lik = lik_cal(X, Y, n, alpha, theta);

	if(lik > maxLik) {
	  maxLik = lik;
	  a_max = a;
	  b_max = b;
	}
      }
    }
    a_left = a_max - step;
    b_left = b_max - step;
    if(a_left < 0)
      a_left = 0.001;
    if(b_left < 0)
      b_left = 0.001;

    a_right = a_max + step;
    b_right = b_max + step;
  } while(step > 0.001);


  alphaML = 1.0 / a_max - 1.0;
  thetaML = b_max;

  maxLik = lik_cal(X, Y, n, alphaML, thetaML);
  lik0 = lik_cal(X, Y, n, alphaML, 0.001);
  LRT = 2.0 * (maxLik - lik0);
  c = maxLik - lik0;
  if(c < 0.001)
    c = 0.001;
  se_theta = thetaML / sqrt(2 * c);
}

/*-----------------------------------------------------------------------*/

static double
lik_cal(const DVector &X, const DVector &Y, int n, double alpha, double theta) {
  int i;
  double lik, K, Q1, Q2, D1, D2;

  lik = 0.0;
  D1 = average(X, n);
  D2 = average(Y, n);

  for(i = 0; i < n; i++) {
    K = K_cal(X(i), Y(i), alpha, D1, D2);
    Q1 = Q_cal(X(i), alpha, D1);
    Q2 = Q_cal(Y(i), alpha, D2);

    lik += log((1.0 - theta) * K + theta * Q1 * Q2);
  }

  return lik;
}

/*----------------------------------------------------------------------*/

static double
K_cal(double X, double Y, double alpha, double D1, double D2) {
  double gammaXYA, gammaA, gammaX, gammaY;
  double K0, c1, c2, ca;

  gammaXYA = gammafun(X + Y + alpha);
  gammaA = gammafun(alpha);
  gammaX = gammafun(X + 1.0);
  gammaY = gammafun(Y + 1.0);

  K0 = gammaXYA / (gammaA * gammaX * gammaY);
  c1 = D1 / (D1 + D2 + alpha);
  c2 = D2 / (D1 + D2 + alpha);
  ca = alpha / (D1 + D2 + alpha);

  return K0 * pow(c1, X) * pow(c2, Y) * pow(ca, alpha);
}

/*----------------------------------------------------------------------*/

static double
Q_cal(double X, double alpha, double D) {
  double gammaXA, gammaA, gammaX;
  double Q0, c1, ca;

  gammaXA = gammafun(X + alpha);
  gammaA = gammafun(alpha);
  gammaX = gammafun(X + 1.0);

  Q0 = gammaXA / (gammaA * gammaX);
  c1 = D / (D + alpha);
  ca = alpha / (D + alpha);

  return Q0 * pow(c1, X) * pow(ca, alpha);
}

/*------------------------------------------------------------------------*/

static double
average(const DVector &A, int n) {
  double sum;
  int k;

  sum = 0.0;
  for(k = 0; k < n; k++) {
    sum += A(k);
  }

  return sum / n;
}

/*--------------------------------------------------------------------------*/

static double
var(const DVector &X, int n) {
  double ave = 0.0, ave2 = 0.0;
  int k;

  for(k = 0; k < n; k++) {
    ave += X(k) / n;
    ave2 += X(k) * X(k) / n;
  }

  return ave2 - ave * ave;
}

/*--------------------------------------------------------------------------*/

static double
covariance(const DVector &X, const DVector &Y, int n) {
  double aveX = 0.0, aveY = 0.0, aveXY = 0.0;
  int k;

  for(k = 0; k < n; k++) {
    aveX += X(k) / n;
    aveY += Y(k) / n;
    aveXY += X(k) * Y(k) / n;
  }

  return aveXY - aveX * aveY;
}

/*--------------------------------------------------------------------------*/

static void
postRatio(const DVector &X, const DVector &Y, int n, DVector &Z, double thetaML, double alphaML) {
  double aveX, aveY, varX, varY, alpha, theta, cov;
  double c1, c2, gamma, gammaX, gammaY, gammaXY, ratio, r0;
  int k;

  aveX = average(X, n);
  aveY = average(Y, n);

  varX = var(X, n);
  varY = var(Y, n);

  cov = covariance(X, Y, n);

  //alpha = (aveX * aveX + aveY * aveY) / (varX + varY - aveX - aveY);
  alpha = alphaML;
  //theta = 1.0 - cov / sqrt((varX - aveX) * (varY - aveY));
  theta = thetaML;
  if(theta >= 1.0) theta = 0.999999;
  c1 = aveY / (aveX + alpha);
  c2 = aveX / (aveY + alpha);
  r0 = (theta / (1.0 - theta)) * pow(1.0 - c1 * c2, alpha);

  for(k = 0; k < n; k++) {
    gamma = gammafun(alpha);
    gammaX = gammafun(X(k) + alpha);
    gammaY = gammafun(Y(k) + alpha);
    gammaXY = gammafun(X(k) + Y(k) + alpha);
    
    ratio = r0 * pow(1.0 + c1, X(k)) * pow(1.0 + c2, Y(k)) * gammaX * gammaY / (gamma * gammaXY);
    Z(k) = ratio / (1.0 + ratio);
    if(Z(k) <= 0.0) Z(k) = 0.0;
  }
}

/*--------------------------------------------------------------------*/

static double
gammafun(double x) {
  x += 1.0;
  const double y = 1.0 / x;
  const double sqrt_2 = 1.41421356237;
  const double sqrt_pi = 1.77245385091;
  
  return (sqrt_2 * sqrt_pi * pow(x, x - 0.5) * exp(-x) *
	  (1.0 +
	   y / 12.0 +
	   y * y / 288.0 -
	   y * y * y * 139.0 / 51840.0 -
	   y * y * y * y * 5.71 / 24883.2 +
	   y * y * y * y * y * 1.63879 / 2090.18880 +
	   y * y * y * y * y * y * 5.246819 / 75246.796800 -
	   y * y * y * y * y * y * y * 5.34703531 / 9029.61561600)) / (x - 1.0);
}

/*----------------------------------------------------------------------*/
