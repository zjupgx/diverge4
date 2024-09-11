#define _CRT_SECURE_NO_WARNINGS
#define WIN32
#include <string>
#include <vector>
#include <iostream>

#include <stdio.h>
#include <stdlib.h>


#include "type_one.h"
#include "gz97.h"
#include "GZf2.h"
#include "tree.h"
#include "matrices.h"
#include "matrix.h"
#include <math.h> 


//#pragma once

using namespace std;

//----------------------------------------------------------------------
extern double p0;
extern DVector site_muta;

//----------------------------------------------------------------------
bool
type_one_compute(const vector<Tree> &trees, const vector<sequence_t> &sequences,
	     vector<vector<double> > &summary, vector<vector<double> > &rets2) 
{
  
  int ntrees = trees.size();
  int seqLength = sequences[0].sequence.size();
  vector<DVector> rets(ntrees);

  DVector freq(jtt_freq, 20);
  DMatrix2D prob(jtt_prob, 20, 20);

  vector< DVector > site_muta_groups;
  site_muta_groups.clear();


  vector< vector<double> > post_probability_array(seqLength, vector<double>(5)); 


  double Alpha1 = 0.0; 
  double Alpha2 = 0.0; 
  double Alpha3 = 0.0; 

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


	double alpha=2.4, se_alpha=0.0;

	double inferred_alpha = 0;
	DVector rets;
	if(!gz97_compute_adapter_for_rvs(true, alpha, alignment,
		treeE, gu_tree, freq, prob, rets, inferred_alpha)) {
			return false;
	}

	if (i == 0)
	{
		Alpha1 = inferred_alpha; 
	}
	else if (i == 1)
	{
		Alpha2 = inferred_alpha; 
	}
	else if (i == 2)
	{
		Alpha3 = inferred_alpha; 
	}

  }

  double p000=0; 

  /*
#ifdef DB
    printf("%d",site_muta_groups[0].size());
#else
	printf("not caught!");
#endif
*/
  int iZero=0;
  for(int i=0;i<seqLength;i++)
  {
	if(site_muta_groups[0](i)<0.001 && site_muta_groups[1](i)<0.001 && site_muta_groups[2](i)<0.001) iZero++;
  }


  double D1 = 0.0; 
  double D2 = 0.0; 
  double D3 = 0.0; 

  vector<double> D(ntrees,0);
  for(int i=0;i<ntrees;i++)
  {
	  for(int j=0;j<seqLength;j++)
	  {
		  D[i] += site_muta_groups[i](j); 
	  }
	  
	  D[i] = D[i] / seqLength; 
  }

  D1 = D[0]; 
  D2 = D[1]; 
  D3 = D[2]; 

 
  vector<DVector> &numSub = rets;
    
  if(!GZf2_compute(numSub, summary, rets2)) 
	  return false; 
  


  double Q0, Q1, Q2, Q3, Qns; 

  double Theta12 = 0.0; 
  double Theta13 = 0.0; 
  double Theta23 = 0.0; 

  double PI0 = 0.0; 

  double Theta1 = 0.0; 
  double Theta2 = 0.0; 
  double Theta3 = 0.0; 


//summary[0] = theta,    summary[1] = se
//summary[2] = r_X,      summary[3] = r_max,
//summary[4] = z_score   (summary[0-4] are model-free)
//summary[5] = thetaML,  summary[6] = alphaML,
//summary[7] = se_theta, summary[8] = LRT
  Theta12 = summary[0][0];
  Theta13 = summary[1][0];
  Theta23 = summary[2][0];

  PI0 = sqrt((1 - Theta12) * (1 - Theta13) * (1 - Theta23)); 


  Theta1 = 1.0 - PI0 / (1 - Theta23); 
  Theta2 = 1.0 - PI0 / (1 - Theta13); 
  Theta3 = 1.0 - PI0 / (1 - Theta12); 

  if (Theta1 < 0)
  {
	  Theta1 = 0; 
  }
  if (Theta2 < 0)
  {
	  Theta2 = 0; 
  }
  if (Theta3 < 0)
  {
	  Theta3 = 0; 
  }


  double Alpha = pow(Alpha1 * Alpha2 * Alpha3, 1.0 / 3.0); 


  double f0 = PI0; 
  double f1 = 1 - Theta23 - PI0; 
  double f2 = 1 - Theta13 - PI0; 
  double f3 = 1 - Theta12 - PI0; 
  double f4 = Theta12 + Theta13 + Theta23 - 2 * (1 - PI0); 


  if (f0 < 0)
  {
	  f0 = 0; 
  }
  if (f1 < 0)
  {
	  f1 = 0; 
  }
  if (f2 < 0)
  {
	  f2 = 0; 
  }
  if (f3 < 0)
  {
	  f3 = 0; 
  }
  if (f4 < 0)
  {
	  f4 = 0; 
  }


  double psi1 = f4 - (1 - Theta1) * Theta2 * Theta3; 
  double psi2 = f4 - (1 - Theta2) * Theta1 * Theta3; 
  double psi3 = f4 - (1 - Theta3) * Theta1 * Theta2; 

  if (psi1 < 0)
  {
	  psi1 = 0; 
  }
  if (psi2 < 0)
  {
	  psi2 = 0; 
  }
  if (psi3 < 0)
  {
	  psi3 = 0; 
  }
  
  for(int i = 0; i < seqLength; i++) 
  {
	  double K123 = GAMMA(site_muta_groups[0](i) + site_muta_groups[1](i) + site_muta_groups[2](i) + Alpha); 
	  K123 /= 0.0 + GAMMA(site_muta_groups[0](i) + 1) * GAMMA(site_muta_groups[1](i) + 1) * GAMMA(site_muta_groups[2](i) + 1) * GAMMA(Alpha); 
	  K123 *= pow(D1 / (D1 + D2 + D3 + Alpha), site_muta_groups[0](i)); 
	  K123 *= pow(D2 / (D1 + D2 + D3 + Alpha), site_muta_groups[1](i)); 
	  K123 *= pow(D3 / (D1 + D2 + D3 + Alpha), site_muta_groups[2](i)); 
	  K123 *= pow(Alpha / (D1 + D2 + D3 + Alpha), Alpha); 
	  
	  double K12 = GAMMA(site_muta_groups[0](i) + site_muta_groups[1](i) + Alpha); 
	  K12 /= GAMMA(site_muta_groups[0](i)+ 1) * GAMMA(site_muta_groups[1](i) + 1) * GAMMA(Alpha); 
	  K12 *= pow(D1 / (D1 + D2 + Alpha), site_muta_groups[0](i)); 
	  K12 *= pow(D2 / (D1 + D2 + Alpha), site_muta_groups[1](i)); 
	  K12 *= pow(Alpha / (D1 + D2 + Alpha), Alpha); 
	  
	  double K13 = GAMMA(site_muta_groups[0](i) + site_muta_groups[2](i) + Alpha); 
	  K13 /= GAMMA(site_muta_groups[0](i) + 1) * GAMMA(site_muta_groups[2](i) + 1) * GAMMA(Alpha); 
	  K13 *= pow(D1 / (D1 + D3 + Alpha), site_muta_groups[0](i)); 
	  K13 *= pow(D3 / (D1 + D3 + Alpha), site_muta_groups[2](i)); 
	  K13 *= pow(Alpha / (D1 + D3 + Alpha), Alpha); 

	  double K23 = GAMMA(site_muta_groups[1](i) + site_muta_groups[2](i) + Alpha); 
	  K23 /= GAMMA(site_muta_groups[1](i) + 1) * GAMMA(site_muta_groups[2](i) + 1) * GAMMA(Alpha); 
	  K23 *= pow(D2 / (D2 + D3 + Alpha), site_muta_groups[1](i)); 
	  K23 *= pow(D3 / (D2 + D3 + Alpha), site_muta_groups[2](i)); 
	  K23 *= pow(Alpha / (D2 + D3 + Alpha), Alpha); 

	  double K1 = GAMMA(site_muta_groups[0](i) + Alpha); 
	  K1 /= GAMMA(site_muta_groups[0](i) + 1) * GAMMA(Alpha); 
	  K1 *= pow(D1 / (D1 + Alpha), site_muta_groups[0](i)); 
	  K1 *= pow(Alpha / (D1 + Alpha), Alpha); 

	  double K2 = GAMMA(site_muta_groups[1](i) + Alpha); 
	  K2 /= GAMMA(site_muta_groups[1](i) + 1) * GAMMA(Alpha); 
	  K2 *= pow(D2 / (D2 + Alpha), site_muta_groups[1](i)); 
	  K2 *= pow(Alpha / (D2 + Alpha), Alpha); 

	  double K3 = GAMMA(site_muta_groups[2](i) + Alpha); 
	  K3 /= GAMMA(site_muta_groups[2](i) + 1) * GAMMA(Alpha); 
	  K3 *= pow(D3/ (D3 + Alpha), site_muta_groups[2](i)); 
	  K3 *= pow(Alpha / (D3 + Alpha), Alpha); 
	  
	  double PXS0 = K123; 
	  double PXS1 = K1 * K23; 
	  double PXS2 = K2 * K13; 
	  double PXS3 = K3 * K12; 
	  double PXS4 = K1 * K2 * K3; 

	  double PX = f0 * PXS0 + f1 * PXS1 + f2 * PXS2 + f3 * PXS3 + f4 * PXS4; 

	  double PS0X = f0 * PXS0 / PX; 
	  double PS1X = f1 * PXS1 / PX; 
	  double PS2X = f2 * PXS2 / PX; 
	  double PS3X = f3 * PXS3 / PX; 
	  double PS4X = f4 * PXS4 / PX; 

	  Q0 = f0 * PXS0 / PX; 
	  Q1 = f1 * PXS1 / PX + psi1 * PXS4 / PX; 
	  Q2 = f2 * PXS2 / PX + psi2 * PXS4 / PX; 
	  Q3 = f3 * PXS3 / PX + psi3 * PXS4 / PX; 
	  Qns = f4 * PXS4 / PX; 
	  
	  post_probability_array[i][0] = PS0X; 
	  post_probability_array[i][1] = PS1X; 
	  post_probability_array[i][2] = PS2X; 
	  post_probability_array[i][3] = PS3X; 
	  post_probability_array[i][4] = PS4X; 
	  
  }

  rets2 = post_probability_array; 
  
  
  return true; 

}


//----------------------------------------------------------------------



// Г 函数
double GAMMA(double x)
{   
	int i; 
	double y, t, s, u; 
	static double a[11] = { 0.0000677106,-0.0003442342, 0.0015397681, 
		                   -0.0024467480, 0.0109736958,-0.0002109075, 
	                        0.0742379071, 0.0815782188, 0.4118402518, 
		                    0.4227843370, 1.0};    
	if (x <= 0.0)    
	{
		printf("err ** x <= 0!\n"); 
		return(-1.0);    
	}    

	y = x;    
	if (y <= 1.0)    
	{
		t = 1.0 / (y * (y + 1.0)); 
		y = y + 2.0; 
	}
	else if (y <= 2.0)
	{
		t = 1.0 / y; 
		y = y + 1.0; 
	}
	else if (y <= 3.0)
	{
		t = 1.0; 
	}
	else
	{
		t = 1.0; 
		while (y > 3.0)
		{
			y = y - 1.0; 
			t = t * y; 
		}
	}

	s = a[0]; 
	u = y - 2.0; 
	for (i = 1; i <= 10; i++)
	{
		s = s * u + a[i]; 
	}
	s = s * t; 

	return(s);

}


//----------------------------------------------------------------------


int trans2(int x)                  // 转换为二进制的数
{
    int m=x, n, b=0, i=0; 

	while (m>0)
	{
		n = m%2;
		m /= 2;
		b += n*pow((double)10, i); 
		i++;
	}
	return b;
}

const char * trans2Str(int x)      // 转换为二进制的字符数
{
	int m=x, n, b=0, i=0; 
	char * w = (char*)"";

	while (m>0)
	{
		n = m%2;
		m /= 2;
		b += n*pow((double)10, i); 
		
		if (n == 0)
		{
			w[i] = '0'; 
		}
		else
		{
			w[i] = '1'; 
		}

		i++;
	}
	w[i] = '\0'; 

	return w;
}

const char * trans16(int x)       // 转换为十六进制的字符数
{
	int h, k, i(0), r;
	char * w = (char*)"";

	h = trans2(x);

	while (h>0)
	{
		r = h % 10000;
		h /= 10000;
		k = (r/1000*8 + r%1000/100*4 + r%100/10*2 + r%10);
		switch (k)
		{
		    case 1:w[i]='1';break;
		    case 2:w[i]='2';break;
    		case 3:w[i]='3';break;
	    	case 4:w[i]='4';break;
		    case 5:w[i]='5';break;
		    case 6:w[i]='6';break;
		    case 7:w[i]='7';break;
		    case 8:w[i]='8';break;
		    case 9:w[i]='9';break;
		    case 10:w[i]='A';break;
		    case 11:w[i]='B';break;
		    case 12:w[i]='C';break;
		    case 13:w[i]='D';break;
		    case 14:w[i]='E';break;
		    case 15:w[i]='F';break;
		    case 16:w[i]='G';break;
		}
		i++;
	}
	w[i] = '\0'; 

	return w; 
}


//----------------------------------------------------------------------


char * strcat(char * s1, char * s2) 
{
	int i,i1,i2;

	i1=strlen(s1);
	i2=strlen(s2);

	char * s3 = new char[i1+i2+1];
	while(*s1!='\0')
		*(s3++)=*(s1++);

	while(*s2!='\0')
		*(s3++)=*(s2++);

	*(s3++)='\0';
	for(i=0;i<i1+i2+1;i++)
		s3--;

	return s3; 
}


//----------------------------------------------------------------------


