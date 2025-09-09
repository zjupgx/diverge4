#define _CRT_SECURE_NO_WARNINGS
#include <string>
#include <vector>
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include "common.h"
#include "type_two.h"
#include "tree.h"
#include "cluster.h"
#include "tree.h"
#include "sequence.h"
#include <pybind11/stl.h>
#include <pybind11/pybind11.h> 
#include <pybind11/numpy.h>
#include <stdexcept>
namespace py = pybind11;
//----------------------------------------------------------------------
using namespace std;
int preprocess(vector<string> &names,vector<Tree> &trees,vector<sequence_t> &sequences) {
	try{
		if (trees.size()<2){
			throw std::invalid_argument("Need to create at least two clusters first.");
		}
		else if (sequences.empty()){
			throw std::invalid_argument("Need to load aligned sequences first.");
		}
	}
	catch(const std::exception& e){
		cout<<e.what()<<endl;
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

bool
process(const vector<Tree> &trees,
		const vector<sequence_t> &sequences,
		const vector<string> &names,
		vector<string> &names2,
		vector<summary_t> &summary2,
		vector<result_t> &results2) {
	vector<vector<double> > summary, rets;
	try
	{
		if(!type_two_compute(trees, sequences, summary, rets)) {
			throw std::runtime_error("Error in type_two_compute module!");
		}
		
		if(summary.empty()) {
			throw std::runtime_error("Summary is empty!");
		}
		else{
			/*char *names[21] = { "Da", "Db", "N", "C",
		"R", "p", "d", "W", "Z", "Alpha ML", "Theta-II", "Theta SE",
			"Gr", "Gc", "h", "Q", "Ar", "PIr", "F00,N", "F00,R", "F00,C"};*/
			
			const char *names[21] = { "Da", "Db", "N", "C",
				"R", "Alpha ML", "Theta-II", "Theta SE", "Ar", "PIr",
				"p", "d", "W", "Z", "Gr", "Gc", "h", "Q", "F00,N", "F00,R", "F00,C" };
			
			size_t nsum = summary[0].size();
			//assert(nsum == 9);
			
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
		
		if(rets.empty()) {
			throw std::runtime_error("Result is empty.");
		}
		else{
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
	}
	catch (const std::exception& e)
	{
		cout << e.what() << endl;
	}
	return true;
}
//----------------------------------------------------------------------




class T2Calculator {
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
		T2Calculator (const std::vector<std::string>& args , const std::vector<std::string>& file_names) : argv(args), filenames(file_names) {
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
		// T2Calculator(){}
		void calculate() {
			
			if (!load_sequences(alignmentFile, sequences)){return;}
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

PYBIND11_MODULE(_type2cpp, m) {
		m.doc() = "\
		pybind11 type2cpp class plugin:\n\
		- create_calculator: Create new T2Calculator\n\
		- calculate: Complete the calculation process\n\
		- _summary: Obtain related parameters summary\n\
		- _r_names: Obtain row_name\n\
		- _kept: Obtain site position\n\
		- _results: Obtain site-specific profile posterior probability results\n\
		";
    py::class_<T2Calculator>(m, "T2Calculator")
        // .def(py::init<>())
				.def(py::init<const std::vector<std::string> &,const std::vector<std::string> &>())
        .def("calculate", &T2Calculator::calculate)
				.def("_summary", &T2Calculator::_summary)
				.def("_argv", &T2Calculator::_argv)
				.def("_alignmentFile", &T2Calculator::_alignmentFile)
				.def("_results", &T2Calculator::_results)
				.def("_kept", &T2Calculator::_kept)
				.def("_r_names", &T2Calculator::_r_names);
		m.def("create_calculator", [](const py::list& args,const py::list& file_names){
			std::vector<std::string> argv;
			std::vector<std::string> filenames;
			for (auto& arg : args) {
					argv.push_back(arg.cast<std::string>());
			}
			for (auto& filename : file_names) {
					filenames.push_back(filename.cast<std::string>());
			}
			return new T2Calculator(argv, filenames);
		});
}