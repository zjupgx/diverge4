#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdexcept>

#include "ancestral_compute.h"
#include "asymmetric_test_one.h"
#include "common.h"
#include "gz97.h"
#include "GZf2.h"
#include "tree.h"
#include "cluster.h"
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h> 
namespace py = pybind11;
//----------------------------------------------------------------------

using namespace std;

bool
preprocess(vector<string> &names,vector<Tree> &trees,vector<sequence_t> &sequences) {
	try{
		if (sequences.empty()){
			throw std::invalid_argument("Need to load aligned sequences first.");
		}
		if (trees.size() !=3 ){
			throw std::invalid_argument("Need to create three clusters first.");
		}
	}
	catch(const char* &e){
		cout<<e<<endl;
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
	try{	
		if(!asymmetric_test_one_compute(trees, sequences,rets)) {
			throw std::runtime_error("Error in asymmetric_test_one_compute module!");
		}
		
		if(rets.empty()) 
		{
			throw std::runtime_error("Result is empty");
		}
		else{
			size_t i, j; 
			names2.push_back("Coefficient of Correlation of Theta"); 
			names2.push_back("Delta Variation"); 

			size_t nrets = rets.size();

			result_t r;
			r.values.resize(rets[0].size());

			for(i = 0; i < nrets; i++) 
			{
				r.pos = i;
				for(j = 0; j < rets[0].size(); j++) 
				{
					r.values[j] = rets[i][j];
				}
				results2.push_back(r);
			}

		}
	}	
	catch(const std::exception& e){
		cout << e.what() << endl;
	}
	return true;
}

class AsymCalculator {
	public:
		std::vector<sequence_t> sequences;
    std::vector<int> kept;
    std::vector<int> removed;
    std::vector<std::string> r_names;
		std::vector<summary_t> summary;
    std::vector<result_t> results;
		std::vector<std::string> argv;
		std::vector<std::string> filenames;
		std::vector<std::string> treeStrs;
		std::string alignmentFile;
		std::vector<Tree> trees;
		std::vector<Tree> trees2;
		AsymCalculator(const std::vector<std::string>& args , const std::vector<std::string>& file_names) : argv(args), filenames(file_names){
			alignmentFile = argv[0];
			for (int ii = 1; ii < argv.size(); ii++) {
					treeStrs.push_back(argv[ii]);
			}
			for (int i = 0; i < argv.size() - 1; i++) {
				Tree tree;
				Tree polytree;
				if (!tree.load(treeStrs[i],filenames[i])){return;} else{
					if(!tree.polyroot(polytree)) {return;}
					else{trees.push_back(polytree);trees2.push_back(tree);}
					}
			}
		}
		// AsymCalculator(){}
		void calculate() {
			if (!load_sequences(alignmentFile, sequences)){return;}
			clean_gaps(sequences, removed, kept);	
			vector<string> t_names; 
			if (!preprocess(t_names,trees2,sequences)) {return;}
			if (!process(trees, sequences, t_names, r_names, summary, results)) {return;}
			return;
		}
		py::array _results() const {
        // convert results to numpy array
        py::array_t<double, py::array::c_style> arr({results.size(), results[0].values.size()});
        auto buffer = arr.mutable_unchecked<2>();
        for (size_t i = 0; i < results.size(); ++i) {
            for (size_t j = 0; j < results[0].values.size(); ++j) {
                buffer(i, j) = results[i].values[j];
            }
        }
        return arr;
    }
		//get argv
		std::vector<std::string>& _argv() {
				return argv;
		}
		//get alignmentFile
		const std::string& _alignmentFile() const {
				return alignmentFile;
		}
		
    const std::vector<std::string>& _r_names() const {
        return r_names;
    }

    const std::vector<int>& _kept() const {
        return kept;
    }
};

PYBIND11_MODULE(_asymcpp, m) {
		m.doc() = "\
		pybind11 Asymmetric Test: Delta Variation class plugin:\n\
		- create_calculator: Create new AsymCalculator\n\
		- calculate: Complete the calculation process\n\
		- _r_names: Obtain row_name\n\
		- _results: Obtain site-specific profile posterior probability results\n\
		";
    py::class_<AsymCalculator>(m, "AsymCalculator")
        // .def(py::init<>())
				.def(py::init<const std::vector<std::string> &,const std::vector<std::string> &>())
        .def("calculate", &AsymCalculator::calculate)
				.def("_argv", &AsymCalculator::_argv)
				.def("_alignmentFile", &AsymCalculator::_alignmentFile)
				.def("_results", &AsymCalculator::_results)
				.def("_r_names", &AsymCalculator::_r_names);
		m.def("create_calculator", [](const py::list& args,const py::list& file_names){
			std::vector<std::string> argv;
			std::vector<std::string> filenames;
			for (auto& arg : args) {
					argv.push_back(arg.cast<std::string>());
			}
			for (auto& filename : file_names) {
					filenames.push_back(filename.cast<std::string>());
			}
        return new AsymCalculator(argv,filenames);
    });
}
