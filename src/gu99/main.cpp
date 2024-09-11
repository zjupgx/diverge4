#define _CRT_SECURE_NO_WARNINGS
#define WIN32
#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <assert.h>
#include "common.h"
#include "gu99.h"
#include "tree.h"
#include "cluster.h"
#include "sequence.h"
#include "matrix.h"
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;
//----------------------------------------------------------------------
using namespace std;


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
		const std::vector<std::string> &names,
		std::vector<std::string> &names2,
		std::vector<summary_t> &summary2,
		std::vector<result_t> &results2) {
	vector<vector<double> > summary, rets;

	try
	{
		if (!gu99_compute(trees, sequences, summary, rets))
		{
			throw std::exception("Error in gu99_compute module!");
		}
	}
	catch (const char *&e)
	{
		cout << e << endl;
	}

	try
	{
		if (summary.empty())
		{
			throw std::exception("Summary is empty.");
		}
		else
		{
			const char *names[9] = {"MFE Theta", "MFE se", "MFE r X", "MFE r max",
															"MFE z score", "ThetaML", "AlphaML", "SE Theta",
															"LRT Theta"};

			size_t nsum = summary[0].size();
			assert(nsum == 9);

			summary_t s;
			s.values.resize(summary.size());

			for (size_t i = 0; i < nsum; i++)
			{
				s.name = names[i];
				for (size_t j = 0; j < summary.size(); j++)
				{
					s.values[j] = summary[j][i];
				}
				summary2.push_back(s);
			}
		}
	}
	catch (const char *&e)
	{
		cout << e << endl;
	}
	try
	{
		if (rets.empty())
		{
			throw std::exception("Result is empty.");
		}
		else
		{
			size_t i, j, ngroups = names.size();
			for (i = 0; i < ngroups; i++)
			{
				for (j = i + 1; j < ngroups; j++)
				{
					string str = names[i] + '/' + names[j];
					names2.push_back(str);
				}
			}

			size_t nrets = rets[0].size();

			result_t r;
			r.values.resize(rets.size());

			for (i = 0; i < nrets; i++)
			{
				r.pos = i;
				for (j = 0; j < rets.size(); j++)
				{
					r.values[j] = rets[j][i];
				}
				results2.push_back(r);
			}
		}
	}
	catch (const char *&e)
	{
		cout << e << endl;
	}
	return true;
}

//vector<int> kept
//vector<result_t> results;
//vector<string> r_names;  row name
//vector<summary_t> summary;
//t_name: tree file name without extension

class Gu99Calculator {
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
		// Gu99Calculator {(const std::vector<std::string>& args):argv(args);(const std::vector<std::string>& file_names):filenames(file_names)} {
		Gu99Calculator (const std::vector<std::string>& args , const std::vector<std::string>& file_names) : argv(args), filenames(file_names) {
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
		// Gu99Calculator(){}
		void calculate() {
			
			if (!load_sequences(alignmentFile, sequences)){printf( "load sequence error" );return;}
			clean_gaps(sequences, removed, kept);
			vector<string> t_names;
			if (!preprocess(t_names,trees2,sequences)) {return;}
			if (!process(trees, sequences, t_names, r_names, summary, results)) {return;}
			return;
		}
		py::list _summary() const {
        py::list params;
        for (const auto& summary_ : summary) {
            py::dict param;
            param["name"] = summary_.name;
            // convert values to numpy array
            py::array_t<double> values(summary_.values.size());
            auto buffer = values.mutable_unchecked<1>();
            for (size_t i = 0; i < summary_.values.size(); ++i) {
                buffer(i) = summary_.values[i];
            }
            param["values"] = values;
            params.append(param);
        }
        return params;
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

PYBIND11_MODULE(_gu99cpp, m) {
		m.doc() = "\
		pybind11 gu99cpp class plugin:\n\
		- create_calculator: Create new Gu99Calculator\n\
		- calculate: Complete the calculation process\n\
		- _summary: Obtain related parameters summary\n\
		- _r_names: Obtain row_name\n\
		- _kept: Obtain site position\n\
		- _results: Obtain site-specific profile posterior probability results\n\
		";
    py::class_<Gu99Calculator>(m, "Gu99Calculator")
        // .def(py::init<>())
				// .def(py::init<const std::vector<std::string>&>())
				.def(py::init<const std::vector<std::string> &,const std::vector<std::string> &>())
        .def("calculate", &Gu99Calculator::calculate)
				.def("_summary", &Gu99Calculator::_summary)
				.def("_argv", &Gu99Calculator::_argv)
				.def("_alignmentFile", &Gu99Calculator::_alignmentFile)
				.def("_results", &Gu99Calculator::_results)
				.def("_kept", &Gu99Calculator::_kept)
				.def("_r_names", &Gu99Calculator::_r_names);
		// m.def("create_calculator", [](const py::list& args,const py:list& file_names) {
    //     std::vector<std::string> argv;
    //     std::vector<std::string> filenames;
		// 		for (auto& arg : args) {
    //         argv.push_back(arg.cast<std::string>());
    //     }
		// 		for (auto& filename : file_names) {
		// 				filenames.push_back(filename.cast<std::string>());
		// 		}
    //     return new Gu99Calculator(argv,filenames);
    // });
		m.def("create_calculator", [](const py::list& args,const py::list& file_names){
			std::vector<std::string> argv;
			std::vector<std::string> filenames;
			for (auto& arg : args) {
					argv.push_back(arg.cast<std::string>());
			}
			for (auto& filename : file_names) {
					filenames.push_back(filename.cast<std::string>());
			}
			return new Gu99Calculator(argv, filenames);
		});
}  