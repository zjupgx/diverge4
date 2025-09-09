#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <assert.h>
#include <thread>
#include <future>
#include <mutex>
#include <exception>
#include "common.h"
#include "gu99.h"
#include "tree.h"
#include "cluster.h"
#include "sequence.h"
#include "matrix.h"
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <stdexcept>
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
		const vector<string> &names,
		vector<string> &names2,
		vector<summary_t> &summary2,
		vector<result_t> &results2) {
	vector<vector<double> > summary, rets;

	try
	{
		if (!gu99_compute(trees, sequences, summary, rets))
		{
			throw std::runtime_error("Error in gu99_compute module!");
		}
	}
	catch (const std::exception& e)
	{
		cout << e.what() << endl;
	}

	try
	{
		if (summary.empty())
		{
			throw std::runtime_error("Summary is empty!");
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
	catch (const std::exception& e)
	{
		cout << e.what() << endl;
	}
	try
	{
		if (rets.empty())
		{
			throw std::runtime_error("Result is empty.");
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
	catch (const std::exception& e)
	{
		cout << e.what() << endl;
	}
	return true;
}

//----------------------------------------------------------------------


//----------------------------------------------------------------------

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

//----------------------------------------------------------------------
// Batch Calculator for Parallel Processing
//----------------------------------------------------------------------

struct TaskResult {
    bool success;
    std::string error_message;
    std::vector<summary_t> summary;
    std::vector<result_t> results;
    std::vector<std::string> r_names;
    std::vector<int> kept;
    int task_id;
    
    TaskResult() : success(false), task_id(-1) {}
};

struct TaskInput {
    std::vector<std::string> argv;
    std::vector<std::string> filenames;
    int task_id;
    
    TaskInput(const std::vector<std::string>& args, 
              const std::vector<std::string>& files, 
              int id) 
        : argv(args), filenames(files), task_id(id) {}
};

class Gu99BatchCalculator {
private:
    std::vector<TaskInput> tasks;
    std::vector<TaskResult> task_results;
    size_t max_threads;
    std::mutex results_mutex;
    
    // Worker function for each thread
    TaskResult calculate_single_task(const TaskInput& task) {
        TaskResult result;
        result.task_id = task.task_id;
        
        try {
            // Create a separate calculator for this task
            Gu99Calculator calculator(task.argv, task.filenames);
            
            // Perform calculation
            calculator.calculate();
            
            // Extract results safely
            result.summary = calculator.summary;
            result.results = calculator.results;
            result.r_names = calculator.r_names;
            result.kept = calculator.kept;
            result.success = true;
            
        } catch (const std::exception& e) {
            result.success = false;
            result.error_message = std::string("Task ") + std::to_string(task.task_id) + 
                                 " failed: " + e.what();
        } catch (...) {
            result.success = false;
            result.error_message = std::string("Task ") + std::to_string(task.task_id) + 
                                 " failed with unknown error";
        }
        
        return result;
    }
    
public:
    Gu99BatchCalculator(size_t num_threads = 0) {
        if (num_threads == 0) {
            max_threads = std::thread::hardware_concurrency();
            if (max_threads == 0) max_threads = 1;
        } else {
            max_threads = num_threads;
        }
    }
    
    void add_task(const std::vector<std::string>& args, 
                  const std::vector<std::string>& filenames) {
        int task_id = static_cast<int>(tasks.size());
        tasks.emplace_back(args, filenames, task_id);
    }
    
    void calculate_batch() {
        if (tasks.empty()) {
            return;
        }
        
        task_results.clear();
        task_results.resize(tasks.size());
        
        // Determine actual number of threads to use
        size_t num_threads = std::min(max_threads, tasks.size());
        
        // Create thread pool
        std::vector<std::future<TaskResult>> futures;
        futures.reserve(tasks.size());
        
        // Submit all tasks
        for (const auto& task : tasks) {
            futures.push_back(
                std::async(std::launch::async, 
                          &Gu99BatchCalculator::calculate_single_task, 
                          this, task)
            );
        }
        
        // Collect results
        for (size_t i = 0; i < futures.size(); ++i) {
            try {
                TaskResult result = futures[i].get();
                if (result.task_id >= 0 && 
                    result.task_id < static_cast<int>(task_results.size())) {
                    task_results[result.task_id] = std::move(result);
                }
            } catch (const std::exception& e) {
                // Handle individual task failure
                TaskResult failed_result;
                failed_result.task_id = static_cast<int>(i);
                failed_result.success = false;
                failed_result.error_message = std::string("Future exception: ") + e.what();
                task_results[i] = std::move(failed_result);
            }
        }
    }
    
    // Python interface methods
    py::list get_all_results() const {
        py::list batch_results;
        
        for (const auto& task_result : task_results) {
            py::dict task_dict;
            task_dict["task_id"] = task_result.task_id;
            task_dict["success"] = task_result.success;
            task_dict["error_message"] = task_result.error_message;
            
            if (task_result.success) {
                // Summary
                py::list summary_list;
                for (const auto& summary_item : task_result.summary) {
                    py::dict summary_dict;
                    summary_dict["name"] = summary_item.name;
                    py::array_t<double> values(summary_item.values.size());
                    auto buffer = values.mutable_unchecked<1>();
                    for (size_t i = 0; i < summary_item.values.size(); ++i) {
                        buffer(i) = summary_item.values[i];
                    }
                    summary_dict["values"] = values;
                    summary_list.append(summary_dict);
                }
                task_dict["summary"] = summary_list;
                
                // Results
                if (!task_result.results.empty()) {
                    py::array_t<double, py::array::c_style> results_array(
                        {task_result.results.size(), task_result.results[0].values.size()});
                    auto results_buffer = results_array.mutable_unchecked<2>();
                    for (size_t i = 0; i < task_result.results.size(); ++i) {
                        for (size_t j = 0; j < task_result.results[0].values.size(); ++j) {
                            results_buffer(i, j) = task_result.results[i].values[j];
                        }
                    }
                    task_dict["results"] = results_array;
                }
                
                // Row names and kept positions
                task_dict["r_names"] = task_result.r_names;
                task_dict["kept"] = task_result.kept;
            }
            
            batch_results.append(task_dict);
        }
        
        return batch_results;
    }
    
    size_t get_num_tasks() const {
        return tasks.size();
    }
    
    size_t get_max_threads() const {
        return max_threads;
    }
    
    void clear_tasks() {
        tasks.clear();
        task_results.clear();
    }
    
    std::vector<bool> get_success_status() const {
        std::vector<bool> status;
        status.reserve(task_results.size());
        for (const auto& result : task_results) {
            status.push_back(result.success);
        }
        return status;
    }
    
    std::vector<std::string> get_error_messages() const {
        std::vector<std::string> messages;
        messages.reserve(task_results.size());
        for (const auto& result : task_results) {
            messages.push_back(result.error_message);
        }
        return messages;
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
		- Gu99BatchCalculator: Parallel batch processing of multiple tasks\n\
		- create_batch_calculator: Create new Gu99BatchCalculator\n\
		";
    
    // Original single calculator
    py::class_<Gu99Calculator>(m, "Gu99Calculator")
        .def(py::init<const std::vector<std::string> &,const std::vector<std::string> &>())
        .def("calculate", &Gu99Calculator::calculate)
				.def("_summary", &Gu99Calculator::_summary)
				.def("_argv", &Gu99Calculator::_argv)
				.def("_alignmentFile", &Gu99Calculator::_alignmentFile)
				.def("_results", &Gu99Calculator::_results)
				.def("_kept", &Gu99Calculator::_kept)
				.def("_r_names", &Gu99Calculator::_r_names);
    
    // New batch calculator
    py::class_<Gu99BatchCalculator>(m, "Gu99BatchCalculator")
        .def(py::init<>(), "Create batch calculator with default thread count")
        .def(py::init<size_t>(), "Create batch calculator with specified thread count")
        .def("add_task", &Gu99BatchCalculator::add_task, 
             "Add a task to the batch queue",
             py::arg("args"), py::arg("filenames"))
        .def("calculate_batch", &Gu99BatchCalculator::calculate_batch,
             "Execute all tasks in parallel")
        .def("get_all_results", &Gu99BatchCalculator::get_all_results,
             "Get results from all completed tasks")
        .def("get_num_tasks", &Gu99BatchCalculator::get_num_tasks,
             "Get total number of tasks")
        .def("get_max_threads", &Gu99BatchCalculator::get_max_threads,
             "Get maximum number of threads")
        .def("clear_tasks", &Gu99BatchCalculator::clear_tasks,
             "Clear all tasks and results")
        .def("get_success_status", &Gu99BatchCalculator::get_success_status,
             "Get success status for all tasks")
        .def("get_error_messages", &Gu99BatchCalculator::get_error_messages,
             "Get error messages for all tasks");
    
    // Factory functions
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
    
    m.def("create_batch_calculator", [](size_t num_threads = 0){
        return new Gu99BatchCalculator(num_threads);
    }, py::arg("num_threads") = 0, "Create a new batch calculator");
}  