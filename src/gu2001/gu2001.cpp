#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "gu2001.h"
#include "gu99.h"

using namespace std;
/*This is the adapter for Gu2001 method
 */
//*******************************************
void Subtree2_cal(double e0[], double e1[], double e2[], int n, double alpha, double * mle, double * post_profile);
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
		  double freq[20], double &prob_result);
double tranPr(int i, int j, double v);
bool alpha_ML(const vector<Tree> &trees, const vector<sequence_t> &sequences, vector<double> &alpha);
//*******************************************


//----------------------------------------------------------------------------
bool gu2001_compute(const std::vector<Tree> &trees,
		    const std::vector<sequence_t> &sequences,
		    std::vector<std::vector<double> > &summary,
		    std::vector<std::vector<double> > &rets2) {
    vector<double> alpha_ml;

    int  Numsite,seqNum,Nseq1,Nseq2,M1,M2,C[3];       // data size
    char **Seq, *X1, *X2,**seqName;      // sequence file

    int * tmpProfile1, * tmpProfile2;
    double * tmpBranchLength1, *tmpBranchLength2;

    double mle[6],  * post_profile;
    int  *Left1,*Right1,*Left2,*Right2;               // tree-profile
    double *Blength1,*Blength2;                       // tree-profile
    int tmp1, tmp2, total_round;
    int i,j,k,l, dump;                                   // computation
    double alpha,logExp[3],*e0,*e1,*e2,freq[20];      // computation

    vector<string> cluster1_taxa, cluster2_taxa;

    if(trees.size() < 2)
	return false;

    //compute alpha_ml first
    if(!alpha_ML(trees, sequences, alpha_ml))
	return false;
    total_round = 0;


    for(tmp1 = 0; tmp1 < trees.size(); tmp1++) {
	for(tmp2 = tmp1 + 1; tmp2 < trees.size(); tmp2++, total_round++) {

	    trees[tmp1].leaf_names(cluster1_taxa);
	    trees[tmp2].leaf_names(cluster2_taxa);

	    Nseq1 = cluster1_taxa.size(); //number of sequences in cluster 1
	    Nseq2 = cluster2_taxa.size(); //number of sequences in cluster 2;
	    Numsite = sequences[0].sequence.size(); // number of sites;

	    // memory allocation
	    seqNum=Nseq1+Nseq2;
	    seqName= (char **)malloc(sizeof(char*)*seqNum);
	    for(i=0;i<seqNum;i++) 
		seqName[i]= (char *)malloc(sizeof(char)*40);
	    Seq= (char **) malloc(sizeof(char*)*seqNum);             
	    for(i=0;i<seqNum;i++)    
		Seq[i]= (char *)malloc(sizeof(char)*(Numsite+1));
	    X1= (char *)malloc(sizeof(char)*Nseq1);
	    X2= (char *)malloc(sizeof(char)*Nseq2);


	    M1=2*Nseq1;            // for a rooted tree, 2n-1 nodes and 2n-2 branches
	    M2=2*Nseq2;
	    Left1= (int *)malloc(sizeof(int)*M1);
	    Left2= (int *)malloc(sizeof(int)*M2);
	    Right1= (int *)malloc(sizeof(int)*M1);
	    Right2= (int *)malloc(sizeof(int)*M2);
   
	    Blength1= (double *)malloc(sizeof(double)*M1);
	    Blength2= (double *)malloc(sizeof(double)*M2);
		   
	    e0= (double *)malloc(sizeof(double)*Numsite);
	    e1= (double *)malloc(sizeof(double)*Numsite);
	    e2= (double *)malloc(sizeof(double)*Numsite);

	    //tree profile input
	    M1=2*Nseq1-1;                  // number of nodes (unrooted)
	    M2=2*Nseq2-1;                  // number of nodes (unrooted)

	    tmpProfile1 = (int *) malloc(sizeof(int) * M1 * 3);
	    tmpBranchLength1 = (double *) malloc(sizeof(double) * M1);

	    tmpProfile2 = (int *) malloc(sizeof(int) * M2 * 3);
	    tmpBranchLength2 = (double *) malloc(sizeof(double) * M2);


	    //initilize the profile
	    ((Tree)trees[tmp1]).generate_tree_profile(tmpProfile1, tmpBranchLength1, Nseq1);
	    //debug




	    ((Tree)trees[tmp2]).generate_tree_profile(tmpProfile2, tmpBranchLength2, Nseq2);



	    for (i=1;i<=M1;i++){
		Left1[i] = tmpProfile1[3 *(i-1)];
		Right1[i] = tmpProfile1[3 *(i-1) + 2];
		Blength1[i] = tmpBranchLength1[i-1];
	    }
  
	    for (i=1;i<=M2;i++){
		Left2[i] = tmpProfile2[3 *(i-1)];
		Right2[i] = tmpProfile2[3 *(i-1) + 2];
		Blength2[i] = tmpBranchLength2[i-1];
	    }
   

	    // sequence file input (Gu-Zhang format)

	    i = 0;
		   
	    //process tree 1
	    for(j = 0, i = 0; j < Nseq1; j++, i++) {
		std::string str = trees[tmp1].id_lookup(j+1);
		for(k = 0; k < str.size(); k++)
		    seqName[i][k] = str[k];
		seqName[i][k] = 0;
			   
		//search for the match sequence based on the name
		for(k = 0; k < sequences.size(); k++) {
		    if(sequences[k].label == str)
			break;
		}
			   
		if(k >= sequences.size())
		    return false; //not found, something wrong

		for(l = 0; l < sequences[k].sequence.size(); l++)
		    Seq[i][l] = sequences[k].sequence[l];
			   
		Seq[i][l]=0;
	    }

	    for(j = 0; j < Nseq2; j++, i++) {
		std::string str = trees[tmp2].id_lookup(j+1);
		for(k = 0; k < str.size(); k++)
		    seqName[i][k] = str[k];
		seqName[i][k] = 0;
		//search for the match sequence based on the name
		for(k = 0; k < sequences.size(); k++) {
		    if(sequences[k].label == str)
			break;
		}
			   
		if(k >= sequences.size())
		    return false; //not found, something wrong

		for(l = 0; l < sequences[k].sequence.size(); l++)
		    Seq[i][l] = sequences[k].sequence[l];
			   
		Seq[i][l]=0;
	    }


	    //***********************************
	    //    Computation
	    //***********************************

   
	    alpha=alpha_ml[total_round];
   
	    for (i=0; i<20;i++)
		freq[i]=0.05;

	    for (k=0; k<Numsite;k++){
//		if(Progress::wasCancelled())
//		    return false;
				
		for (j=0; j<Nseq1; j++) X1[j]=Seq[j][k];
		for (j=0; j<Nseq2; j++) X2[j]=Seq[j+Nseq1-1][k];
		if(!ExpMarkov2(X1,Nseq1,Left1,Right1,Blength1,X2,Nseq2,Left2,Right2,Blength2,logExp,alpha,freq))
		  return false;
		e0[k]=logExp[0];
		e1[k]=logExp[1];
		e2[k]=logExp[2];

//		Progress::increment(100);
	    }
			
	    post_profile=(double *)malloc(sizeof(double)*Numsite);
	    Subtree2_cal(e0,e1,e2,Numsite,alpha, mle, post_profile);
			
	    vector<double> sum(4);
	    sum[0] = mle[0];
	    sum[1] = mle[1];
	    sum[2] = mle[2];
	    sum[3] = mle[5];
			
	    summary.push_back(sum);
			
	    vector<double> Z2(Numsite);
			
	    for(i = 0; i < Numsite; i++)
		Z2[i] = post_profile[i];
			
	    rets2.push_back(Z2);
			
	    free(post_profile);

	    ///deallocate the memory
	    free(tmpBranchLength2);
	    free(tmpProfile2);
	    free(tmpBranchLength1);
	    free(tmpProfile1);
		    
	    free(e2);
	    free(e1);
	    free(e0);
	    free(Blength2);
	    free(Blength1);
	    free(Right2);
	    free(Right1);
	    free(Left2);
	    free(Left1);
	    free(X2);
	    free(X1);

	    for(i=0;i<seqNum;i++)    
		free(Seq[i]);
	    free(Seq);
	    for(i=0;i<seqNum;i++) 
		free(seqName[i]);
	    free(seqName);
	}
    }
    
    return true;
}

//*******************************************
void Subtree2_cal(double e0[], double e1[], double e2[], int n, double alpha, double * mle, double * post_profile)   
{
    /*for estimating theta*/
    mle[0]=alpha;                           // gamma parameter
    mle[1]=theta2(e0,e1,e2,n);              // theta
#ifdef DEBUG_DIVERGE
    printf("\n");
    printf("mle-5 is OK \n");
#endif

    mle[2]=SE_theta2(e0,e1,e2,mle[1],n);    // standard error of theta
    mle[3]=loglik_value(e0,e1,e2,0.0,n);      // log-lik unde H0 (theta=0)
    mle[4]=loglik_value(e0,e1,e2,mle[1],n); // log-lik
    mle[5]=mle[4]-mle[3];                   // log-lik difference

    PostProfile(e0,e1,e2,mle[1],n,post_profile);   // posterior profile
    // Rate difference profile
    OutputSubtree2(mle, post_profile,n);    // Output
}

//--------------------------------------------------------------------------
void OutputSubtree2(double mle[], double post_profile[], int n)
{
#ifdef DEBUG_DIVERGE
    FILE *fp;
    int  i;

    fp = stdout;
    fprintf(fp,"Results from the program Subtree:\n");
	
    fprintf(fp,"THIS IS FOR TESTING:\n");
    fprintf(fp,"ML Estimates:\n");

    fprintf(fp,"alpha=  ");
    fprintf(fp,"%f", mle[0]);
    fprintf(fp,"\n");

    fprintf(fp,"theta=  ");
    fprintf(fp,"%f", mle[1]);
    fprintf(fp,"\n");

    fprintf(fp,"s.e. theta=  ");
    fprintf(fp,"%f", mle[2]);
    fprintf(fp,"\n");

    fprintf(fp,"LRT=  ");
    fprintf(fp,"%f", mle[5]);
    fprintf(fp,"\n");
	
    fprintf(fp,"Posterior profile:\n");
    for (i=0; i<n; i++) {
	fprintf(fp,"%f", post_profile[i]);
	fprintf(fp,"\n");
    }

#endif
}

//---------------------------------------------------------------------------
void PostProfile(double e0[], double e1[], double e2[],
		 double theta, int n, double post_profile[])
{
    double p,q;
    int k;
    // printf("\n");
    for (k=0; k<n; k++) {
	q=theta*exp(e1[k]+e2[k]);
	p=q+(1.0-theta)*exp(e0[k]);
	post_profile[k]=q/p;
    }
}

double loglik_value(double e0[], double e1[],double e2[], double theta, int n)  // log-lik
{
    int k;
    double p,lnlik;

    lnlik=0.0;
	
    for (k=0;k<n;k++) {
	p=theta*exp(e1[k]+e2[k])+(1.0-theta)*exp(e0[k]);
	//if (p<0.000000001) lnlik=lnlik-20.723;	
	lnlik=lnlik+log(p);
    }
	
    return(lnlik);
}


//*********************************************************
bool ExpMarkov2(char X1[], int n1, int Left1[], int Right1[],
		double Blength1[], char X2[], int n2,int Left2[], int Right2[],
		double Blength2[], double logExp[3], double alpha, double freq[20])
//*********************************************************
{
    int i,k,K;
    double rate[16],ln_p1,ln_p2,Exp_s[3],v;
    
    K=8;
    gammad(rate,alpha,K); //K-discrete model (rates)
    for(i=0;i<3;i++) 
	Exp_s[i]=0.0;

    for (k=0;k<K;k++) {
	v=rate[k];
	if (!MarkovProb(X1,n1,v,Left1,Right1,Blength1,freq, ln_p1))
	  return false;
	if ( !MarkovProb(X2,n2,v,Left2,Right2,Blength2,freq, ln_p2) )
	  return false;
	 
	Exp_s[0]=Exp_s[0]+exp(ln_p1+ln_p2)/K;
	Exp_s[1]=Exp_s[1]+exp(ln_p1)/K;
	Exp_s[2]=Exp_s[2]+exp(ln_p2)/K;
    }

    for(i=0;i<3;i++) 
	logExp[i]=log(Exp_s[i]);

    return true;
}



/*----------------------------------------------------------------*/

void gammad(double rate[], double alpha,int K)
{
    double z,c,step,value;
    int i,j;
	
    step=0.01;
    z=0.0;
    
    for (i=1;i<K;i++) {
	j=0;
    loop:   z=z+step;
	j=j+1;
	value=igamma(z,alpha);
	c=1.0*i/K;

	if (value<c) goto loop;
		
	rate[i]=z/alpha;
    }

    rate[0]=0.1*rate[1];
    rate[K-1]=2.0*rate[K-1];
}

/*------------------------------------------------------------------*/
double igamma(double z, double alpha)
{
    double I,a, b,x_a, x_b;
    int k,N;
	
    N=100;
    I=0.0;
	
    for (k=0;k<N;k++) {
	x_a=z*k/N;
	x_b=z*(1.0+k)/N;
	if (k==0) x_a=0.001;
	a=pow(x_a,alpha-1.0)*exp(-x_a);
	b=pow(x_b,alpha-1.0)*exp(-x_b);
	I=I+0.5*(a+b)*z/N;
    }

    I=I/gammafun(alpha);
    return (I);
}

//----------------------------------------------------------------
double gammafun(double x)
{

    double value,y;

    x=x+1.0;
    value=sqrt(2.00)*sqrt(3.1416)*pow(x,x-0.5)*exp(-x);
    y=1.0/x;
    value=value*(1.0+y/12.0+y*y/288.0-y*y*y*139.0/51840.0
                 -y*y*y*y*5.71/24883.2+y*y*y*y*y*1.63879/2090.18880
                 + y*y*y*y*y*y*5.246819/75246.796800
                 - y*y*y*y*y*y*y*5.34703531/9029.61561600);

    value=value/(x-1.0);
    return (value);
}

//---------------------------------------------------------------------------
bool MarkovProb(char X[], int n, double rate,
		  int Left[], int Right[], double Blength[],
		  double freq[20], double & prob_result)
// return the log-prob value
{
    double **Prob;
    double bL, bR, ProbL, ProbR,p_markov,lnP;
    int i,k,ia,iL,iR;
    int leftNode,rightNode;
    char aminoAcid[20]={'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'};

#ifdef DEBUG_DIVERGE
    FILE *fp;
    fp = stdout;
    for(i = 0; i < n; i++) {
	fprintf(fp, "X[%d]  = %c\n", i, X[i]);
    }
#endif

    // post-order traversal algorithm in a rooted tree
    Prob=(double **) malloc(sizeof(double*)*(2*n));
    for (i=0;i<2*n;i++)
	Prob[i]= (double *) malloc(sizeof(double)*20);

    // prob[k][ia] is the prob at node k conditional of amino acid ia
    // Initial set: Prob=0
    for (k=0;k<=2*n-1;k++) {
	for (ia=0;ia<20;ia++) {
	    Prob[k][ia]=0.0;
	}
    }

    // Set Prob=1 for external nodes: k=1,...,n
    for (k=1;k<=n;k++){
	for (ia=0;ia<20;ia++){
	    if (X[k-1]==aminoAcid[ia]) 
		Prob[k][ia]=1.0; //X starting at k=0 to k=n-1
	}
    }

#ifdef DEBUG_DIVERGE
    fprintf(fp, "for external node\n");
    for(k = 0; k <= (2*n-1); k++) {
	for(ia = 0; ia < 20; ia++)
	    fprintf(fp, "%f ", Prob[k][ia]);
	fprintf(fp, "\n");
    }
#endif

    // Iteration for internal nodes
    for (k=n+1;k<=2*n-1;k++){
	leftNode=Left[k];
	rightNode=Right[k];
	bL=Blength[leftNode]*rate;
	bR=Blength[rightNode]*rate;
	// Node k has two offsprings, leftNode and rightNode, with branch
	// lengths bL and bR, respectively.

#ifdef DEBUG_DIVERGE
	fprintf(fp, "k = %d, lefNode = %d, rightNode = %d, bL = %f, bR= %f\n", k, leftNode, rightNode, bL, bR);
#endif

	for (ia=0;ia<20;ia++) {
	    for (iL=0;iL<20;iL++) {
		for (iR=0;iR<20;iR++){
		    ProbL=tranPr(iL,ia,bL);
		    ProbR=tranPr(iR,ia,bR);
		    Prob[k][ia]=Prob[k][ia] + Prob[leftNode][iL]*Prob[rightNode][iR]*ProbL*ProbR;
		}
	    }
	}
    }

#ifdef DEBUG_DIVERGE
    fprintf(fp, "for internal node\n");
    for(k = 0; k <= (2*n-1); k++) {
	for(ia = 0; ia < 20; ia++)
	    fprintf(fp, "%f ", Prob[k][ia]);
	fprintf(fp, "\n");
    }


#endif

    p_markov=0.0;
    for (ia=0;ia<20;ia++)
	p_markov=p_markov+freq[ia]*Prob[2*n-1][ia];
    if(p_markov == 0.0) {
#ifdef DEBUG_DIVERGE
	fprintf(fp, "site with p_markov 0\n");
#endif
//	QMessageBox::warning(NULL, "Gu2001 Method", "Gu2001 method failed.  Please re-check input sequence data and tree information.");
	return false;
    }
    lnP=log(p_markov);

    prob_result = lnP;
    return true;
}

//---------------------------------------------------------------------------


double theta2(double e0[], double e1[], double e2[], int n)
{
    double theta,error,H,a,b,c,f,g,grad;
    int k;
    /* iteration for theta*/
    theta=0.5;
    error=0.0001;

 loop: 
    grad=0.0;
    H=0.0;
	
    for ( k=0; k<n; k++ ) {
	a=exp(e1[k]);
	b=exp(e2[k]);
	c=exp(e0[k]);
	f=a*b-c;
	g=f/(theta*f+c);
	grad=grad+g;
	H=H-g*g;
#ifdef DEBUG_DIVERGE
	printf("k = %d, a = %f, b = %f, c = %f, f = %f, g = %f, grad = %f, H = %f\n", k, a, b, c, f, g, grad, H);
#endif	    
    }
	
    theta=theta-grad/H;
#ifdef DEBUG_DIVERGE
    printf(" theta = %f \n",theta);
    printf(" grad = %f \n", grad);
    printf(" H = %f \n", H);
#endif
	
    if(grad/H > 0) {
	if (grad/H>error) goto loop;
	return(theta);
    }
    else {
	if ((-1 * grad/H)>error) goto loop;
	return(theta);
    }
}



//----------------------------------------------------------------------------
double SE_theta2(double e0[], double e1[], double e2[], double theta, int n)           // standard error of theta
{
    double se,H,a,b,c,f,g,grad;
    int k;

    grad=0.0;
    H=0.0;

    for (k=0;k<n; k++){
	a=exp(e1[k]);
	b=exp(e2[k]);
	c=exp(e0[k]);
	f=a*b-c;
	g=f/(theta*f+c);
	grad=grad+g;
	H=H-g*g;
    }

    se=sqrt(-1.0/H);
    return(se);
}


/*--------------------------------------------------------------------------*/
//  Transition probability from i to j, given branch length v.
double tranPr(int i, int j, double v)
{
    double p;

    /* One-parameter model for amino acid changes */
    if (i==j) 
	p=0.05+0.95*exp(-v/0.95);
    else 
	p=0.05-0.05*exp(-v/0.95);

    return(p);
}


/*---------------------------------------------------------------------------*/
bool
alpha_ML(const vector<Tree> &trees, const vector<sequence_t> &sequences, vector<double> &alpha) {
    vector<vector<double> > summary, rets2;
    int i;

    if(!gu99_compute(trees, sequences, summary, rets2))
	return false;
	
    for(i = 0; i < summary.size(); i++)
	alpha.push_back(summary[i][6]);

    return true;
}



