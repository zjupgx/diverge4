#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <assert.h>
#include "common.h"
#include "tree.h"
#include "cluster.h"
#include "matrix.h"
#include "false_discovery_rate.h"
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <stdexcept>
namespace py = pybind11;

//----------------------------------------------------------------------
using namespace std;
//----------------------------------------------------------------------
bool
preprocess(vector<string> &names,vector<Tree> &trees,vector<sequence_t> &sequences) {
  try{
		if (trees.size()<2){
			throw std::invalid_argument("Need to create at least two clusters first.");
		}
	}
	catch(const char* &e){
		cout<<e<<endl;
	}
	try{
		if (sequences.empty()){
			throw std::invalid_argument("Need to load aligned sequences first.");
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
		       vector<result_t> &results1, 
			   vector<result_t> &results2, 
			   vector<result_t> &results3, 
			   vector<result_t> &results4) 
{
  vector<vector<double> > summary, rets1, rets2, rets3, rets4;
  try{
		if (!false_discovery_rate_compute(trees, sequences, rets1, rets2, rets3, rets4)) 
		{
			throw std::runtime_error("Error in false_discovery_rate_compute module!");
		}


		if(rets1.empty()) 
		{
			throw std::runtime_error("Result_1 is empty!");
		}
		else{

			size_t i, j, ngroups = names.size();

			names2.push_back("Probability Cutoff"); 
				for(i = 0; i < ngroups; i++) 
			{
				for(j = i + 1; j < ngroups; j++) 
				{
							string str = "FDR of " + names[i] + '/' + names[j];
								names2.push_back(str);
						}
				}

				size_t nrets = rets1[0].size();
			
				result_t r;
				r.values.resize(rets1.size());

				for(i = 0; i < nrets; i++) {
						r.pos = i;
						for(j = 0; j < rets1.size(); j++) {
							r.values[j] = rets1[j][i];
						}
						results1.push_back(r);
				}

		}

		if(rets2.empty()) 
		{
			throw std::runtime_error("Result_2 is empty!");
		}
		else{

			size_t i, j; 

			if (rets1.empty())
			{
				names2.push_back("Probability Cutoff"); 
				names2.push_back("False Discovery Rate"); 
			}


			size_t nrets = rets2[0].size();

			result_t r;

			r.values.resize(rets2.size());

			for(i = 0; i < nrets; i++) {
				r.pos = i;
				for(j = 0; j < rets2.size(); j++) {
					r.values[j] = rets2[j][i];
				}
				results2.push_back(r);
			}

		}


		if(rets3.empty()) 
		{
			throw std::runtime_error("Result_3 is empty!");
		}
		else{
			size_t i, j;

				size_t nrets = rets3[0].size();
			
				result_t r;
				r.values.resize(rets3.size());

				for(i = 0; i < nrets; i++) {
						r.pos = i;
						for(j = 0; j < rets3.size(); j++) {
							r.values[j] = rets3[j][i];
						}
						results3.push_back(r);
				}
		}

		if(rets4.empty()) 
		{
			throw std::runtime_error("Result_4 is empty!");
		}
		else{
			size_t i, j; 

			size_t nrets = rets4[0].size();

			result_t r;
			r.values.resize(rets4.size());

			for(i = 0; i < nrets; i++) {
				r.pos = i;
				for(j = 0; j < rets4.size(); j++) {
					r.values[j] = rets4[j][i];
				}
				results4.push_back(r);
			}
		}
		}
	catch (const std::exception& e)
	{
		cout << e.what() << endl;
	}
  return true;

}


class FdrCalculator {
	public:
		std::vector<sequence_t> sequences;
    std::vector<int> kept;
    std::vector<int> removed;
    std::vector<std::string> r_names;
    std::vector<result_t> results1, results2, results3, results4;
		std::vector<std::string> argv;
		std::vector<std::string> filenames;
		std::vector<std::string> treeFiles;
		std::string alignmentFile;
		std::vector<Tree> trees;
		std::vector<Tree> trees2;
		FdrCalculator(const std::vector<std::string>& args , const std::vector<std::string>& file_names) : argv(args), filenames(file_names) {
			alignmentFile = argv[0];
			for (int ii = 1; ii < argv.size(); ii++) {
					treeFiles.push_back(argv[ii]);
			}
			for (int i = 0; i < argv.size() - 1; i++) {
				Tree tree;
				Tree polytree;
				if (!tree.load(treeFiles[i],filenames[i])){return;} else{
					if(!tree.polyroot(polytree)) {return;}
					else{trees.push_back(polytree);trees2.push_back(tree);}
					}
			}
		}
		// FdrCalculator(){}
		void calculate() {
			
			if (!load_sequences(alignmentFile, sequences)){return;}
			clean_gaps(sequences, removed, kept);	
			vector<string> t_names; 
			if (!preprocess(t_names,trees2,sequences)) {return;}
			if (!process(trees, sequences, t_names, r_names, results1, results2, results3, results4)) {return;}
			return;
			}

		py::array _results1() const {
        // convert results1 to numpy array
        py::array_t<double, py::array::c_style> arr({results1.size(), results1[0].values.size()});
        auto buffer = arr.mutable_unchecked<2>();
        for (size_t i = 0; i < results1.size(); ++i) {
            for (size_t j = 0; j < results1[0].values.size(); ++j) {
                buffer(i, j) = results1[i].values[j];
            }
        }
        return arr;
    }
		py::array _results2() const {
        // convert results2 to numpy array
        py::array_t<double, py::array::c_style> arr({results2.size(), results2[0].values.size()});
        auto buffer = arr.mutable_unchecked<2>();
        for (size_t i = 0; i < results2.size(); ++i) {
            for (size_t j = 0; j < results2[0].values.size(); ++j) {
                buffer(i, j) = results2[i].values[j];
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

PYBIND11_MODULE(_fdrcpp, m) {
		m.doc() = "\
		pybind11 fdr class plugin:\n\
		- create_calculator: Create new FdrCalculator\n\
		- calculate: Complete the calculation process\n\
		- _results1: Obtain type one fdr\n\
		- _results2: Obtain type two fdr\n\
		- _r_names: Obtain row_name\n\
		- _kept: Obtain site position\n\
		";
    py::class_<FdrCalculator>(m, "FdrCalculator")
        // .def(py::init<>())
				.def(py::init<const std::vector<std::string> &,const std::vector<std::string> &>())
        .def("calculate", &FdrCalculator::calculate)
				.def("_argv", &FdrCalculator::_argv)
				.def("_results1", &FdrCalculator::_results1)
				.def("_results2", &FdrCalculator::_results2)
				.def("_kept", &FdrCalculator::_kept)
				.def("_r_names", &FdrCalculator::_r_names);
		m.def("create_calculator", [](const py::list& args,const py::list& file_names){
			std::vector<std::string> argv;
			std::vector<std::string> filenames;
			for (auto& arg : args) {
					argv.push_back(arg.cast<std::string>());
			}
			for (auto& filename : file_names) {
					filenames.push_back(filename.cast<std::string>());
			}
			return new FdrCalculator(argv, filenames);
		});
}
