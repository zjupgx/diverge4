#include <string>
#include <vector>

#include <stdio.h>
#include <assert.h>

#include <iostream>
#include <fstream>


#include "common.h"
#include "gu2001.h"
#include "tree.h"
#include "cluster.h"
#include "tree.h"
#include "sequence.h"

//----------------------------------------------------------------------

using namespace std;

vector<sequence_t> sequences;
vector<int> removed;
vector<Tree> trees;
vector<int> kept;


//----------------------------------------------------------------------
bool
preprocess(vector<string> &names);


//----------------------------------------------------------------------

void
calculate() {
	int i;
	bool trees_branch_zero;
	vector<string> t_names;
	
	trees_branch_zero = false;
	for( i = 0; i < trees.size(); i++ ) {
		if(!trees[i].has_branch_length()) {
			trees_branch_zero = true;
			break;
		}
	}
	if (trees_branch_zero) {
		//QMessageBox::warning(NULL, "Gu2001 Method", "Gu2001 method failed. The input trees do not contain branch length!");
		return;
	}
	
	if(!preprocess(t_names)) return;
	
	
	vector<string> names;
	vector<summary_t> summary;
	vector<result_t> results;
}


//----------------------------------------------------------------------

bool
preprocess(vector<string> &names) {
	if(trees.size() < 2) {
		//    QMessageBox::information(NULL, "Gu2001 Method", "Need to create at least two clusters first.");
		return false;
	}
	if(sequences.empty()) {
		//    QMessageBox::information(NULL, "Gu2001 Method", "Need to load aligned sequences first.");
		return false;
	}
	
	names.resize(trees.size());
	for(size_t i=0; i<trees.size(); i++) {
		string str = trees[i].filename();
		
		size_t j;
		
		j = str.find_last_of('/');
		if(j != string::npos) str.erase(0, j+1);
		
		j = str.find_last_of('.');
		if(j != string::npos) str.erase(j, string::npos);
		
		names[i] = str;
	}
	
	return true;
}

//----------------------------------------------------------------------

bool
process(const vector<Tree> &trees,
		const vector<sequence_t> &sequences,
		const vector<string> &names,
		vector<string> &names2,
		vector<summary_t> &summary2,
		vector<result_t> &results2) {
	vector<vector<double> > summary, rets;
	
	if(!gu2001_compute(trees, sequences, summary, rets)) {
		return false;
	}
	
	if(!summary.empty()) {
		char *names[4] = { "AlphaML", "ThetaML", "SE Theta", "LRT Theta" };
		
		size_t nsum = summary[0].size();
		assert(nsum == 4);
		
		summary_t s;
		s.values.resize(summary.size());
		
		for(size_t i=0; i<nsum; i++) {
			s.name = names[i];
			for(size_t j=0; j<summary.size(); j++) {
				s.values[j] = summary[j][i];
			}
			summary2.push_back(s);
		}
	}
	
	if(!rets.empty()) {
		size_t i, j, ngroups = names.size();
		for(i = 0; i < ngroups; i++) {
			for(j = i + 1; j < ngroups; j++) {
				string str = names[i] + '/' + names[j];
				names2.push_back(str);
			}
		}
		
		size_t nrets = rets[0].size();
		
		result_t r;
		r.values.resize(rets.size());
		
		for(i = 0; i < nrets; i++) {
			r.pos = i;
			for(j = 0; j < rets.size(); j++) {
				r.values[j] = rets[j][i];
			}
			results2.push_back(r);
		}
	}
	
	return true;
}

//----------------------------------------------------------------------


//----------------------------------------------------------------------
int main(int argc, char *argv []){
	if(argc<3){
		printf("Usage: %s alignment.aln output.txt cluseter1.tree cluseter2.tree ... clusetern.tree",argv[0]);
		exit(0);
	}

	vector<string> treeFiles;
	
	string alignmentFile=argv[1]; 
	char *output=argv[2];
	for(int ii=3;ii<argc;ii++){
		treeFiles.push_back(argv[ii]);
	}
	
	
	FILE *fp;
	fp=fopen(output, "w");
	
	for(int i=0;i<argc-3;i++){
		Tree tree;
		if(!tree.load(treeFiles[i])) {
			printf("load tree error:%s\n", treeFiles[i]);
			exit(0);
		}else{
			trees.push_back(tree);	
		}
	}
	
	if(!load_sequences(alignmentFile, sequences)) {
		printf("unable to load sequences alignment: %s\n", alignmentFile);
		exit(0);
	}
	clean_gaps(sequences, removed, kept);
	
	//calculate gu99		
	vector<string> t_names;
	if(!preprocess(t_names)) return;
	
	vector<string> names;
	vector<summary_t> summary;
	vector<result_t> results;
	
	if(!process(trees, sequences, t_names, names, summary, results)) return;
	
	int n = trees.size(), j, k;
	vector<string> cluster_names;
	//DMatrix2D theta_values(n, n);
	
	//parameters
	fprintf(fp,"Parameter");
	for(k=0, i=0; i<n; i++) {
		for(j=i+1; j<n; j++) {
			//theta_values(j, i) = theta_values(i, j) = summary[k++].values[5];
			//		theta_values(j, i) = theta_values(i, j) = summary[0].values[k++];
			fprintf(fp,"\t%s\/%s",trees[i].filename().c_str(), trees[j].filename().c_str());
		}
		//cluster_names.push_back(trees[i].filename().c_str());
	}
	
	fprintf(fp,"\n");
	
	int iflds=summary.size();
	for(i=0;i<iflds;i++){
		fprintf(fp,"%s",summary[i].name.c_str());
		for(j=0;j<n;j++){
			fprintf(fp,"\t%f",summary[i].values[j]);
		}
		fprintf(fp,"\n");
	}
	
	fprintf(fp,"\n");
	fprintf(fp,"\n");
	
	
	//result
	fprintf(fp,"Position");
	for(k=0, i=0; i<n; i++) {
		for(j=i+1; j<n; j++) {
			//theta_values(j, i) = theta_values(i, j) = summary[k++].values[5];
			//		theta_values(j, i) = theta_values(i, j) = summary[0].values[k++];
			fprintf(fp,"\t%s\/%s",trees[i].filename().c_str(), trees[j].filename().c_str());
		}
		//cluster_names.push_back(trees[i].filename().c_str());
	}
	fprintf(fp,"\n");
	
	int ipos=results.size();
	for(i=0;i<ipos;i++){
		fprintf(fp,"%d",kept[i]+1);
		for(j=0;j<n;j++){
			fprintf(fp,"\t%f",results[i].values[j]);
		}
		fprintf(fp,"\n");
	}
	
	fprintf(fp,"\n");
	fprintf(fp,"\n");
	fclose(fp);
	
	exit(0);
}

