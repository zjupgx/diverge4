#define _CRT_SECURE_NO_WARNINGS
#define WIN32
#include <string>
#include <vector>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "effective_number.h"
#include "gz97.h"
#include "GZf2.h"

#include "tree.h"
#include "matrices.h"
#include "matrix.h"

#include <iostream>

#include <algorithm>

#ifdef WIN32
#include <minmax.h>
#endif

//----------------------------------------------------------------------

using namespace std;

//----------------------------------------------------------------------

enum SITE_CLASSIFY {IDENTICAL, SAMEGROUP, DIFFGROUP};

//----------------------------------------------------------------------

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

bool type_two_compute(const vector<Tree> &trees, const vector<sequence_t> &sequences,
	vector<vector<double> > &summary, vector<vector<double> > &rets2);
bool
	type_two_divergence(const Tree &t1, const Tree &t2, const vector<sequence_t> &s1,
	const vector<sequence_t> &s2, const double alpha, const vector<double> &tree1_numsub,
	const vector<double> &tree2_numsub, vector<double> &summary, 
	vector<double> &rets);
bool GZ97_compute_adapter(const vector<Tree> &trees, const vector<sequence_t> &sequences,
	vector<double> &alpha, vector< vector<double> > &numSub);
int amino_acid_group(char a);
enum SITE_CLASSIFY amino_acid_compare(char a, char b);
int num_site_with_amino_acid_classify(enum SITE_CLASSIFY cl, const sequence_t & seq1,
	const sequence_t &seq2);
void site_profile();
bool ancestral_inference(const Tree &t, const vector<sequence_t> &seq, 
	sequence_t &ans_seq);
int F(const vector<double> &changes1, const vector<double> &changes2,
	const sequence_t &root1, const sequence_t &root2, enum SITE_CLASSIFY c);
int G(sequence_t &s1, sequence_t &s2, enum SITE_CLASSIFY c);

void construct_seq(const Tree &t, const vector<sequence_t> &seqs, vector<sequence_t> &s);

//----------------------------------------------------------------------

extern double p0;
extern DVector site_muta;

bool effective_number_compute(const vector<Tree> &trees, const vector<sequence_t> &sequences,
	     vector<vector<double> > &rets1, vector<vector<double> > &rets2, 
		 int &EffectiveNumberType1, int &EffectiveNumberType2) 
{
  
//  Type One Computation. 
	int ntrees = trees.size();

	int seqLength = sequences[0].sequence.size();
	vector<DVector> rets(ntrees);

	DVector freq(jtt_freq, 20);
	DMatrix2D prob(jtt_prob, 20, 20);

	vector< DVector > site_muta_groups;
	site_muta_groups.clear();

	for(int i = 0; i < ntrees; i++) 
	{
		int treeE;


		vector<string> taxa;
		trees[i].leaf_names(taxa);

		CMatrix2D alignment(taxa.size(), sequences[0].sequence.size());

		{
			// generate sub_sequences and sub_seqNames
			vector<string>::const_iterator i;
			vector<sequence_t>::const_iterator j;
			int i2;
			for(i2 = 0, i = taxa.begin(); i != taxa.end(); i++) {
				for(j = sequences.begin(); j != sequences.end(); j++) {
					if(j->label == *i) {
						for(int k=0; k<(int)j->sequence.size(); k++) {
							alignment(i2, k) = j->sequence[k];
						}
						i2++;
						break;
					}
				}

				if(j == sequences.end()) {
					abort();
				} 
			}
		}

		string tree_str;
		if(!trees[i].gen_str_wrt_seq(taxa, tree_str)) 
			return false;

		IMatrix2D gu_tree;

		if(!parse_tree(tree_str.c_str(), alignment.rows(), treeE, gu_tree)) {
			return false;
		}
		if(!gz97_compute(true, 2.4, alignment, treeE, gu_tree, freq, prob, rets[i])) {
			return false;
		}

		site_muta_groups.push_back(site_muta);

	}
	
	double p000=0; 

	int iZero=0;
	for(int i=0;i<seqLength;i++)
	{
		if(site_muta_groups[0](i)<0.001 && site_muta_groups[1](i)<0.001) iZero++;
	}


	vector<DVector> &numSub = rets; 


// To computer effective number of sites (Type One) in a iteration process. 
// Some codes are from GZf2_compute() function. 

	vector<double> thetaVector;
	vector<double> seVector;

	int groupNum = numSub.size();
	int siteNum = numSub[0].size();

	DVector Z(siteNum);
	for(int i = 0; i < groupNum; i++) 
	{
		for(int j = i + 1; j < groupNum; j++) 
		{
			
			siteNum = numSub[0].size();
			EffectiveNumberType1 = 0; 
			thetaVector.clear();
			seVector.clear(); 


			DVector &X = numSub[i];
			DVector &Y = numSub[j];

			vector<double> sum(9);
			// sum[0] = theta,    sum[1] = se,
			// sum[2] = r_X,      sum[3] = r_max,
			// sum[4] = z_score   (sum[0-4] are model-free),
			// sum[5] = thetaML,  sum[6] = alphaML,
			// sum[7] = se_theta, sum[8] = LRT

			//double theta, se, r_X, r_max, z_score;
			covGZ(X, Y, siteNum, 0, sum[0], sum[1], sum[2], sum[3], sum[4]);
			likGZ(X, Y, siteNum, sum[5], sum[6], sum[7], sum[8]);
			postRatio(X, Y, siteNum, Z, sum[5], sum[6]);

			thetaVector.push_back(sum[0]);
			seVector.push_back(sum[1]);


			double maxPostProb = - 1.0; 
		    int posPostProb = - 1; 

			for (int k = 0; k < Z.size(); k++)
			{
				if (maxPostProb < Z(k))
				{
					maxPostProb = Z(k); 
					posPostProb = k + 1; 
				}
			}


    		DVector tempX = X; 
        	DVector tempY = Y; 

			while (sum[0] > sum[1] && siteNum > 0)
			{

				DVector tempV = tempX; 

				tempX.resize(siteNum - 1); 

				int ni = 0; 
				for (int k = 0; k < siteNum; k++)
				{
					if (k != posPostProb - 1)
					{
						tempX(ni) = tempV(k); 
						ni++; 
					}
				}

				tempV = tempY; 

				tempY.resize(siteNum - 1); 

				ni = 0; 
				for (int k = 0; k < siteNum; k++)
				{
					if (k != posPostProb - 1)
					{
						tempY(ni) = tempV(k); 
						ni++; 
					}
				}

				tempV.clear(); 


				siteNum--; 
			    EffectiveNumberType1++; 

				//double theta, se, r_X, r_max, z_score;
				covGZ(tempX, tempY, siteNum, 0, sum[0], sum[1], sum[2], sum[3], sum[4]);
				likGZ(tempX, tempY, siteNum, sum[5], sum[6], sum[7], sum[8]);
				postRatio(tempX, tempY, siteNum, Z, sum[5], sum[6]);


				maxPostProb = - 1.0; 
				posPostProb = - 1; 

				for (int k = 0; k < siteNum; k++)
				{
					if (maxPostProb < Z(k))
					{
						maxPostProb = Z(k); 
						posPostProb = k + 1; 
					}
				}
				

				thetaVector.push_back(sum[0]);
				seVector.push_back(sum[1]);


			}

			rets1.push_back(thetaVector); 
			rets1.push_back(seVector); 

		}

	}

	

//  Type Two Computation. 
	ntrees = trees.size();

	int j, totalcount;

	//find out the alpha_ML and number of substitution 
	vector<double> alpha;
	vector< vector<double> > numsubs;
	GZ97_compute_adapter(trees, sequences, alpha, numsubs);
	///////////////////

	totalcount = 0;


	siteNum = numSub[0].size(); 


	thetaVector.clear();
	seVector.clear();


	vector<int> site_used_Array; 

	for (int i = 0; i < siteNum; i++)
	{
	   site_used_Array.push_back(1); 
	}

	
	for(int i = 0; i < ntrees; i++) 
	{
		for(j = i + 1; j < ntrees; j++, totalcount++) 
		{
			
			siteNum = numSub[0].size();
			EffectiveNumberType2 = 0; 
			thetaVector.clear();
			seVector.clear();


        	vector<double> &X = numsubs[i];
    		vector<double> &Y = numsubs[j];


			vector<double> sum;
			vector<double> tempRets;
			vector<sequence_t> s1, s2;
			Tree t1 = trees[i];
			Tree t2 = trees[j];
			construct_seq(t1, sequences, s1);
			construct_seq(t2, sequences, s2);
	        

//          To computer effective number of sites (Type One) in a iteration process. 
//          Some codes are from type_two_divergence() function. 

			int site_num;
			int tree1_seqnum, tree2_seqnum;

			site_num = s1.front().sequence.size();
			tree1_seqnum = s1.size();
			tree2_seqnum = s2.size();

			//---
			double Da, Db;


			Da = 0;
			Db = 0;
			for(int k = 0; k < site_num; k++)
				Da += X[k];
			Da = Da/site_num;

			for(int k = 0; k < site_num; k++)
				Db += Y[k];
			Db = Db/site_num;
			//---
			sequence_t a, b;
			int N, C, R;
			
			ancestral_inference(t1, s1, a);
			ancestral_inference(t2, s2, b);
			
			N = num_site_with_amino_acid_classify(IDENTICAL, a, b);
			C = num_site_with_amino_acid_classify(SAMEGROUP, a, b);
			R = num_site_with_amino_acid_classify(DIFFGROUP, a, b);


			double P = 1 - ((double)N)/site_num;
			double F00_N = ((double)F(X, Y, a, b, IDENTICAL))/site_num;
			double F00_C = ((double)F(X, Y, a, b, SAMEGROUP))/site_num;
			double F00_R = ((double)F(X, Y, a, b, DIFFGROUP))/site_num;

			double d = 0;
			double W = 0; 
			double theta = 0.0; 
			double tmp = 0.0;
			
			while(1) {
				d = -1 * log(1 - (P-theta)/(1-theta));
				W = pow(alpha[totalcount]/(Da+Db+d+alpha[totalcount]), alpha[totalcount]);
				theta = 1 - F00_N/W;
////			if (theta <= 0.0)
			    if (theta <= 0.0 || d <= 0.0 || W <= 0.0)
					break;
				else {
///					if (fabs(tmp - theta) < 0.001)
					if (fabs(tmp - theta) < 0.01)
						break;
					else
						tmp = theta;
				}
			}

			if (d > P)
			{
				d = P; 
				W = pow(alpha[totalcount]/(Da+Db+d+alpha[totalcount]), alpha[totalcount]); 
				theta = 1 - F00_N/W; 
			}
			

			double Z = pow(alpha[totalcount]/(Da+Db+alpha[totalcount]), alpha[totalcount]);

			double theta_se = sqrt( ((double)F00_N * (1-F00_N))/(N+C+R) )/W;

			double Gr = ((double)R)/(R+C);
			double Gc = ((double)C)/(R+C);

			double PIr=0.31213;
			double Ar = (P*Gr - (1-theta) * PIr * (P-theta))/theta;
			//double Ar = (F00_R - (Z-W) * Gr)/(theta * W);
			//double PIr = ( Gr * Z - F00_R )/((1 - theta) * W);


			double Ac = 1 - Ar;
			double PIc = 1 - PIr;

			double h = d/(Da + Db + d + alpha[totalcount]);
			double Q = ((double)R)/(1 + R);

			for(int k = 0; k < site_num; k++) {
				double R;
				char tmp1 = a.sequence[k];
				char tmp2 = b.sequence[k];
				switch(amino_acid_compare(tmp1, tmp2)) {
				case IDENTICAL:
					R = 0.0;
					break;
				case SAMEGROUP:
					R = (theta * Ac)/((1-theta) * PIc * (1-pow((1-h), (X[k] + Y[k] + alpha[totalcount]))));
					break;
				case DIFFGROUP:
					R = (theta * Ar)/((1-theta) * PIr * (1-pow((1-h), (X[k] + Y[k] + alpha[totalcount]))));
					break;
				}

        		tempRets.push_back(R); 

			}

			
   

			thetaVector.push_back(theta); 
			seVector.push_back(theta_se); 


			double maxPostProb = - 1.0; 
			int posPostProb = - 1; 

			for ( int k = 0; k < tempRets.size(); k++)
			{
				if (maxPostProb < tempRets[k])
				{
					maxPostProb = tempRets[k]; 
					posPostProb = k + 1; 
				}
			}
			
			
    		vector<double> tempX = X; 
    		vector<double> tempY = Y; 

			while (theta > theta_se && site_num > 0)
			{

				vector<double> tempVX = tempX; 
				vector<double> tempVY = tempY; 

				tempX.resize(site_num - 1); 
				tempY.resize(site_num - 1); 
				
				int ni = 0; 
				for (int k = 0; k < site_num; k++)
				{
					if (k != posPostProb - 1)
					{
						tempX[ni] = tempVX[k]; 
						tempY[ni] = tempVY[k]; 
						ni++; 
					}
				}
				

				site_num--; 
			    EffectiveNumberType2++; 

				
				Da = 0;
    			Db = 0;
				for( int k = 0; k < site_num; k++)
					Da += tempX[k];
				Da = Da/site_num;

				for( int k = 0; k < site_num; k++)
					Db += tempY[k];
				Db = Db/site_num;
				//---
				

				P = 1 - ((double)N)/site_num;

				
				int resultInt1 = 0;
				int resultInt2 = 0; 
				int resultInt3 = 0; 

				int total_num_site = a.sequence.size(); 
				int calcu_num_site = 0; 
				int sitePos = 0; 

				while (sitePos < total_num_site && calcu_num_site < site_num)
				{
					if (site_used_Array[sitePos] == 1)
					{
						if (tempX[calcu_num_site] == 0 && tempY[calcu_num_site] == 0)
						{
							if (amino_acid_compare(a.sequence.at(sitePos), b.sequence.at(sitePos)) == IDENTICAL)
							{
								resultInt1++; 
							}
							if (amino_acid_compare(a.sequence.at(sitePos), b.sequence.at(sitePos)) == SAMEGROUP)
							{
								resultInt2++; 
							}
							if (amino_acid_compare(a.sequence.at(sitePos), b.sequence.at(sitePos)) == DIFFGROUP)
							{
								resultInt3++; 
							}
						}

						calcu_num_site++; 
					}

					sitePos++; 
				}

				F00_N = ((double) resultInt1) / site_num; 
				F00_C = ((double) resultInt2) / site_num; 
				F00_R = ((double) resultInt3) / site_num; 


				d = 0;
				W = 0; 
				theta = 0.0; 
				tmp = 0.0;

				while(1) {
					d = -1 * log(1 - (P-theta)/(1-theta));
					W = pow(alpha[totalcount]/(Da+Db+d+alpha[totalcount]), alpha[totalcount]);
					theta = 1 - F00_N/W;
////				if (theta <= 0.0)
					if (theta <= 0.0 || d <= 0.0 || W <= 0.0)
						break;
					else {
///						if (fabs(tmp - theta) < 0.001)
						if (fabs(tmp - theta) < 0.01)
							break;
						else
							tmp = theta;
					}
				}

				if (d > P)
				{
				    d = P; 
					W = pow(alpha[totalcount]/(Da+Db+d+alpha[totalcount]), alpha[totalcount]); 
					theta = 1 - F00_N/W; 
				}


				Z = pow(alpha[totalcount]/(Da+Db+alpha[totalcount]), alpha[totalcount]);

				theta_se = sqrt( ((double)F00_N * (1-F00_N))/(N+C+R) )/W;

				Gr = ((double)R)/(R+C);
				Gc = ((double)C)/(R+C);

				PIr=0.31213;
				Ar = (P*Gr - (1-theta) * PIr * (P-theta))/theta;
				//Ar = (F00_R - (Z-W) * Gr)/(theta * W);
				//PIr = ( Gr * Z - F00_R )/((1 - theta) * W);


				Ac = 1 - Ar;
				PIc = 1 - PIr;

				h = d/(Da + Db + d + alpha[totalcount]);
				Q = ((double)R)/(1 + R);

				
				tempRets.clear(); 

				for( int k = 0; k < site_num; k++) {
					double R;
					char tmp1 = a.sequence[k];
					char tmp2 = b.sequence[k];
					switch(amino_acid_compare(tmp1, tmp2)) {
					case IDENTICAL:
						R = 0.0;
						break;
					case SAMEGROUP:
						R = (theta * Ac)/((1-theta) * PIc * (1-pow((1-h), (tempX[k] + tempY[k] + alpha[totalcount]))));
						
						break;
					case DIFFGROUP:
						R = (theta * Ar)/((1-theta) * PIr * (1-pow((1-h), (tempX[k] + tempY[k] + alpha[totalcount]))));
						break;
						
					}

	        		tempRets.push_back(R); 

				}
			
    		
				thetaVector.push_back(theta); 
				seVector.push_back(theta_se); 


				double maxPostProb = - 1.0; 
				int posPostProb = - 1; 
				
				for (int k = 0; k < tempRets.size(); k++)
				{
					if (maxPostProb < tempRets[k])
					{
						maxPostProb = tempRets[k]; 
						posPostProb = k + 1; 
					}
				}
				

				total_num_site = a.sequence.size(); 
				calcu_num_site = 0; 
				sitePos = -1; 

				while (sitePos < total_num_site && calcu_num_site < posPostProb)
				{
					sitePos++; 

					if (site_used_Array[sitePos] == 1)
					{
						calcu_num_site++; 
					}
				}

				if (sitePos < total_num_site)
				{
					site_used_Array[sitePos] = -1; 
				}

			}

			rets2.push_back(thetaVector); 
			rets2.push_back(seVector);

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
		double alpha, alpha1, alpha2;
		double xM, yM, xV, yV, cov;
		double z_X, z_max;

		xM = average(X, n);
		yM = average(Y, n);

		xV = var(X, n);
		yV = var(Y, n);
		cov = covariance(X, Y, n);

		theta = 1.0 - cov / sqrt((xV - xM) * (yV - yM));
///		alpha = (xM * xM + yM * yM) / (xV + yV - xM - yM);
		alpha1 = xM * xM / (xV - xM); 
		alpha2 = yM * yM / (yV - yM); 
		alpha = min(alpha1, alpha2); 


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






