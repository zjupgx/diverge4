#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "common.h"
#include "tree.h"
#include "cluster.h"
#include "matrix.h"
#include "effective_number.h"
#include <pybind11/stl.h>
#include <pybind11/pybind11.h> 
#include <pybind11/numpy.h>     
#include <stdexcept>
namespace py = pybind11;
//----------------------------------------------------------------------
using namespace std;

bool
preprocess(vector<string> &names,vector<Tree> &trees,vector<sequence_t> &sequences)
{
	try{
		if (sequences.empty()){
			throw std::invalid_argument("Need to load aligned sequences first.");
		}
		if (trees.size() !=2 ){
			throw std::invalid_argument("Need to create two clusters first.");
		}
	}
	catch(const std::exception& e){
		cout << e.what() << endl;
	}

	names.resize(trees.size());
	for (size_t i = 0; i < trees.size(); i++) {
		string str = trees[i].filename();

		size_t j;

		j = str.find_last_of('/');
		if (j != string::npos) str.erase(0, j + 1);

		j = str.find_last_of('.');
		if (j != string::npos) str.erase(j, string::npos);

		names[i] = str;
	}

	return true;
}


bool
process(const vector<Tree>& trees,
	const vector<sequence_t>& sequences,
	const vector<string>& names,
	vector<string>& names2,
	vector<result_t>& results1,
	vector<result_t>& results2,
	int& EffectiveNumberType1, int& EffectiveNumberType2)
{
	vector<vector<double> > rets1, rets2;
	try{
		if (!effective_number_compute(trees, sequences, rets1, rets2, EffectiveNumberType1, EffectiveNumberType2))
			{
				throw std::runtime_error("Error in effective_number_compute module!");
			}

		if (rets1.empty()){throw std::runtime_error("Results_1 is empty!");}
		else{

			size_t i, j;

			names2.push_back("Theta");
			names2.push_back("Theta Standard Deviation");

			size_t nrets = rets1[0].size();

			result_t r;
			r.values.resize(rets1.size());

			for (i = 0; i < nrets; i++) {
				r.pos = i;
				for (j = 0; j < rets1.size(); j++) {
					r.values[j] = rets1[j][i];
				}
				results1.push_back(r);
			}
		}

		if (rets2.empty()){throw std::runtime_error("Results_2 is empty!");}
		else{

			size_t i, j;

			if (rets1.empty())
			{
				names2.push_back("Theta");
				names2.push_back("Theta Standard Deviation");
			}

			size_t nrets = rets2[0].size();

			result_t r;
			r.values.resize(rets2.size());

			for (i = 0; i < nrets; i++) {
				r.pos = i;
				for (j = 0; j < rets2.size(); j++) {
					r.values[j] = rets2[j][i];
				}
				results2.push_back(r);
			}

		}
	}	
	catch (const std::exception& e)
	{
		cout << e.what() << endl;
	}
	return true;

}

class EffCalculator {
	public:
		std::vector<sequence_t> sequences;
    std::vector<int> kept;
    std::vector<int> removed;
    std::vector<std::string> r_names;
		std::vector<summary_t> summary;
    std::vector<result_t> results1,results2;
		std::vector<std::string> argv;
		std::vector<std::string> filenames;
		std::vector<std::string> treeStrs;
		std::string alignmentFile;
		std::vector<Tree> trees;
		std::vector<Tree> trees2;
		EffCalculator(const std::vector<std::string>& args , const std::vector<std::string>& file_names) : argv(args), filenames(file_names)  {
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
		// EffCalculator(){}
		void calculate() {
			
			if (!load_sequences(alignmentFile, sequences)){return;}
			clean_gaps(sequences, removed, kept);	
			vector<string> t_names; 
			if (!preprocess(t_names,trees2,sequences)) {return;}
			int EffectiveNumberType1 = 0;
			int EffectiveNumberType2 = 0;
			if (!process(trees, sequences, t_names, r_names, results1, results2, EffectiveNumberType1, EffectiveNumberType2)) {
				return;
			}
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

PYBIND11_MODULE(_effectivecpp, m) {
		m.doc() = "\
		pybind11 Effective Number of Sites Related to Functional class plugin:\n\
		- create_calculator: Create new EffCalculator\n\
		- calculate: Complete the calculation process\n\
		- _r_names: Obtain row_name\n\
		- _results1: Obtain Type-one Site-Specific Profile results\n\
		- _results2: Obtain Type-two Site-Specific Profile results\n\
		";
    py::class_<EffCalculator>(m, "EffCalculator")
        // .def(py::init<>())
				.def(py::init<const std::vector<std::string> &,const std::vector<std::string> &>())
        .def("calculate", &EffCalculator::calculate)
				.def("_argv", &EffCalculator::_argv)
				.def("_alignmentFile", &EffCalculator::_alignmentFile)
				.def("_results1", &EffCalculator::_results1)
				.def("_results2", &EffCalculator::_results2)
				.def("_r_names", &EffCalculator::_r_names);
		m.def("create_calculator", [](const py::list& args,const py::list& file_names){
			std::vector<std::string> argv;
			std::vector<std::string> filenames;
			for (auto& arg : args) {
					argv.push_back(arg.cast<std::string>());
			}
			for (auto& filename : file_names) {
					filenames.push_back(filename.cast<std::string>());
			}
			return new EffCalculator(argv, filenames);
		});
}

