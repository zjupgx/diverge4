from . import _gu2001cpp,_gu99cpp,_type2cpp,_asymcpp,_fdrcpp,_effectivecpp,_rvscpp,_typeOneAnalysiscpp
from Bio import Phylo
import pandas as pd

import matplotlib.pyplot as plt
import numpy as np

def pre_check_tree(*tree_files):
    for tree_file in tree_files:
        tree = Phylo.read(tree_file, 'newick')
        for clade in tree.find_clades():
            clade.branch_length = 1
        # get tree depth
        tree_depth = max(tree.depths().values())
        if tree_depth < 3:
            raise Exception(f'tree depth is less than 3, please check your tree file{tree_file}')
            return False    
    return True

class Gu99():
    def __init__(self,aln_file,*tree_files) -> None:
        self.input = [aln_file]
        self.input.extend(iter(tree_files))
        if pre_check_tree(*tree_files):
            self._calculate()
            self._summary = self.summary()
    def _calculate(self):
        """
        Create a new Gu99 Calculator and Complete the calculation steps.
        """
        self.calculator = _gu99cpp.create_calculator(self.input)
        self.calculator.calculate()
    def summary(self):
        calculator = self.calculator
        summary_results = pd.DataFrame(columns=calculator._r_names())
        dicts = calculator._summary()
        for dict in dicts:
          summary_results.loc[dict["name"]] = dict["values"]
        return summary_results
    def results(self):
        calculator = self.calculator
        results = pd.DataFrame(calculator._results(),columns=calculator._r_names(),index=[i+1 for i in calculator._kept()])
        return results
    def fundist(self):
        if len(self.input)<=3:
            return None
        else:
            import numpy as np
            n = len(self.input)-1  
            theta = np.zeros((n, n))
            

            k = 0
            for i in range(n):
                for j in range(i + 1, n):
                    theta[j, i] = theta[i, j] = self._summary.iloc[0,k]
                    k += 1

            # 功能距离
            B = 0.0
            d = np.zeros((n, n))
            for i in range(n - 1):
                for j in range(i + 1, n):
                    d[j, i] = d[i, j] = -np.log(1.0 - theta[i, j])
                    B += d[i, j] / (n - 1)

            dist_results = []
            for k in range(n):
                bk = 0.0
                for i in range(n - 1):
                    if i == k:
                        continue
                    for j in range(i + 1, n):
                        if j == k:
                            continue
                        bk += d[i, j] / (n - 2)
                dist_results.append(B - bk)
            
            return pd.DataFrame(dist_results,index=self.calculator._r_names())
class Gu2001():
    def __init__(self,aln_file,*tree_files) -> None:
        self.input = [aln_file]
        self.input.extend(iter(tree_files))
        if pre_check_tree(*tree_files):
            self._calculate()
    def _calculate(self):
        """
        Create a new Gu2001 Calculator and Complete the calculation steps.
        """
        self.calculator = _gu2001cpp.create_calculator(self.input)
        self.calculator.calculate()
    def _help(self):
        print(help(_gu2001cpp))
    def summary(self):
        calculator = self.calculator
        summary_results = pd.DataFrame(columns=calculator._r_names())
        dicts = calculator._summary()
        for dict in dicts:
          summary_results.loc[dict["name"]] = dict["values"]
        summary_results.index.name = "Parameters"
        return summary_results
    def results(self):
        calculator = self.calculator
        results = pd.DataFrame(calculator._results(),columns=calculator._r_names(),index=[i+1 for i in calculator._kept()])
        results.index.name = "Position"
        return results 
class Type2():
    def __init__(self,aln_file,*tree_files) -> None:
        self.input = [aln_file]
        self.input.extend(iter(tree_files))
        if pre_check_tree(*tree_files):
            self._calculate()
    def _calculate(self):
        """
        Create a new type2 Calculator and Complete the calculation steps.
        """
        self.calculator = _type2cpp.create_calculator(self.input)
        self.calculator.calculate()
    def summary(self):
        calculator = self.calculator
        summary_results = pd.DataFrame(columns=calculator._r_names())
        dicts = calculator._summary()
        for dict in dicts:
          summary_results.loc[dict["name"]] = dict["values"]
        summary_results.index.name = "Parameters"
        return summary_results
    def results(self):
        calculator = self.calculator
        results = pd.DataFrame(calculator._results(),columns=calculator._r_names(),index=[i+1 for i in calculator._kept()])
        results.index.name = "Position"
        return results
    def _help(self):
        print(help(_type2cpp))
class Asym():
    def __init__(self,aln_file,*tree_files) -> None:
        self.input = [aln_file]
        self.input.extend(iter(tree_files))
        self._calculate()
    def _calculate(self):
        """
        Create a new asym Calculator and Complete the calculation steps.
        """
        self.calculator = _asymcpp.create_calculator(self.input)
        self.calculator.calculate()
    def results(self):
        calculator = self.calculator
        results = pd.DataFrame(calculator._results(),columns=calculator._r_names(),index = list(range(1,np.size(calculator._results(),0)+1)))
        results.index.name = "Cluster Number of Outgroup" 
        return results
    def _help(self):
        print(help(_asymcpp))
class Effective():
    def __init__(self,aln_file,*tree_files) -> None:
        self.input = [aln_file]
        self.input.extend(iter(tree_files))
        self._calculate()
        self.type1_effective_number = self.calculator._results1().shape[0]
        self.type2_effective_number = self.calculator._results2().shape[0]
        print(f"Type1 Effective Number of Sites is {self.type1_effective_number}, Type2 Effective Number of Sites is {self.type2_effective_number}")
    def _calculate(self):
        """
        Create a new Effective number of site Calculator and Complete the calculation steps.
        """
        self.calculator = _effectivecpp.create_calculator(self.input)
        self.calculator.calculate()
    def type1_results(self):
        calculator = self.calculator
        results = pd.DataFrame(calculator._results1(),columns=calculator._r_names())
        results.index.name = "Number"
        return results
    def type2_results(self):
        calculator = self.calculator
        results = pd.DataFrame(calculator._results2(),columns=calculator._r_names())
        results.index.name = "Number"
        return results
    def _help(self):
        print(help(_effectivecpp))
class Fdr():
    def __init__(self,aln_file,*tree_files) -> None:
        self.input = [aln_file]
        self.input.extend(iter(tree_files))
        self._calculate()
    def _calculate(self):
        """
        Create a new fdr Calculator and Complete the calculation steps.
        """
        self.calculator = _fdrcpp.create_calculator(self.input)
        self.calculator.calculate()
    def type1_results(self):
        calculator = self.calculator
        results = pd.DataFrame(calculator._results1(),columns=calculator._r_names())
        results = results.set_index(results.columns[0],drop=True)
        return results
    def type2_results(self):
        calculator = self.calculator
        results = pd.DataFrame(calculator._results2(),columns=calculator._r_names())
        results = results.set_index(results.columns[0],drop=True)
        return results
    def _help(self):
        print(help(_fdrcpp))     
class Rvs():

    def __init__(self,aln_file,*tree_files) -> None:
        self.input = [aln_file]
        self.input.extend(iter(tree_files))
        self._calculate()
    def _calculate(self):
        """
        Create a new rvs Calculator and Complete the calculation steps.
        """
        self.calculator = _rvscpp.create_calculator(self.input)
        self.calculator.calculate()
    def summary(self):
        calculator = self.calculator
        summary_results = pd.DataFrame(columns=calculator._s_names())
        dicts = calculator._summary()
        for dict in dicts:
          summary_results.loc[dict["name"]] = dict["values"]
        summary_results.index.name = "Parameters"
        return summary_results
    def results(self):
        """
        (Parameter:Xk: Number of Changes; Rk: Posterior Mean of Evolutionary Rate)
        """
        calculator = self.calculator
        results = pd.DataFrame(calculator._results(),columns=calculator._r_names(),index=[i+1 for i in calculator._kept()])
        return results
    def _help(self):
        print(help(_rvscpp))
class TypeOneAnalysis():
    """
    A class for type one analysis
    Note:Under the two-state model(functional divergence unrelated F0 or related F1),\
		there are eight possible combined states for three duplicate clusters, \
		which can be reduced to five nondegenerate patterns.\n\
		S0=(F0, F0, F0) means no type-one divergence occurred in any clusters.\n\
		S1=(F1, F0, F0) means type-one functional divergence occurred only in cluster 1, \n\
		similarly: S2 =(F0, F1,F0) and S3 =(F0, F0, F1). \n\
		The final pattern S4 is for the rest of four states, each of which has two or three clusters that have experienced type-one functional divergence.
    """
    def __init__(self,aln_file,*tree_files) -> None:
        self.input = [aln_file]
        self.input.extend(iter(tree_files))
        self._calculate()
    def _calculate(self):
        """
        Create a new TypeOneAnalysis Calculator and Complete the calculation steps.
        """
        self.calculator = _typeOneAnalysiscpp.create_calculator(self.input)
        self.calculator.calculate()
    def summary(self):
        calculator = self.calculator
        summary_results = pd.DataFrame(columns=calculator._s_names())
        dicts = calculator._summary()
        for dict in dicts:
          summary_results.loc[dict["name"]] = dict["values"]
        summary_results.index.name = "Parameters"
        return summary_results
    def results(self):
        calculator = self.calculator
        results = pd.DataFrame(calculator._results(),columns=calculator._r_names(),index=[i+1 for i in calculator._kept()])
        results.index.name = "Position"
        return results
    def _help(self):
        print(help(_typeOneAnalysiscpp))
if __name__ == "__main__":
    # test gu99
    T = Gu99("../test_data/CASP.aln", "../test_data/cl1.tree", "../test_data/cl2.tree","../test_data/cl3.tree")
    summary = T.summary()
    results = T.results()
    print(summary)
    print("------------")
    print(results)
    # test gu2001
    T = Gu2001("../test_data/CASP.aln", "../test_data/treelength/cl1.tree", "../test_data/treelength/cl2.tree","../test_data/treelength/cl3.tree")
    summary = T.summary()
    results = T.results()
    print(summary)
    print("------------")
    print(results)
    # test type2 
    T = Type2("../test_data/CASP.aln", "../test_data/cl1.tree", "../test_data/cl2.tree","../test_data/cl3.tree")
    summary = T.summary()
    results = T.results()
    print(summary)
    print("------------")
    print(results)
    # test asym
    T = Asym("../test_data/CASP.aln", "../test_data/cl1.tree", "../test_data/cl2.tree","../test_data/cl3.tree")
    results = T.results()
    print(results)
    # test effective
    T = Effective("../test_data/CASP.aln", "../test_data/cl1.tree", "../test_data/cl2.tree")
    # T._help()
    results1 = T.type1_results()
    results2 = T.type2_results()
    print(results1)
    print("------------")
    print(results2)
    # test fdr
    T = Fdr("../test_data/CASP.aln", "../test_data/cl1.tree", "../test_data/cl2.tree","../test_data/cl3.tree")
    results1 = T.type1_results()
    results2 = T.type2_results()
    print(results1)
    print("------------")
    print(results2)
    # test rvs
    T = Rvs("../test_data/CASP.aln", "../test_data/cl1.tree")
    summary = T.summary()
    results = T.results()
    print(summary)
    print("------------")
    print(results)
    # test type one analysis
    T = TypeOneAnalysis("../test_data/CASP.aln", "../test_data/cl1.tree", "../test_data/cl2.tree","../test_data/cl3.tree")
    summary = T.summary()
    results = T.results()
    print(summary)
    print("------------")
    print(results)
    pass