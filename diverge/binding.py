from . import _gu2001cpp,_gu99cpp,_type2cpp,_asymcpp,_fdrcpp,_effectivecpp,_rvscpp,_typeOneAnalysiscpp
from Bio import Phylo
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

def read_tree(tree_file):
    try:
        tree = Phylo.read(tree_file, 'newick')
    except ValueError:
        try:
            tree = Phylo.read(tree_file, 'nexus')
        except ValueError:
            raise ValueError(f"Tree file '{tree_file}' is not in newick or nexus format.")
    # set brach length to 1 for all clade
    # for clade in tree.get_nonterminals():
    #     clade.branch_length = 1
    return tree

def get_colnames(r_names):
    columns_names = []
    for name in r_names:
        # if os.path.isfile(name):
        #     col_name = os.path.basename(name)
        #     # col_name = os.path.splitext(col_name)[1]
        # else:
        #     col_name = name
        # columns_names.extend([col_name])
        columns_names.extend([name])
    return columns_names

def load_tree_file(tree_file,check=True):
    tree = read_tree(tree_file)
    if check:
        check_tree(tree)
    return tree.format("newick").strip("\n")


def check_tree_file(*tree_files):
    for tree_file in tree_files:
        tree = Phylo.read(tree_file, 'newick')
        tree_depth = max([len(tree.trace(tree.root,clade)) for clade in tree.get_terminals()])
        if tree_depth < 3:
            raise Exception(f'tree depth is less than 3, please check your tree file{tree_file}')    
    return True

def check_tree(tree):
    
    tree_depth = max([len(tree.trace(tree.root,clade)) for clade in tree.get_terminals()])
    if tree_depth < 3:
        raise Exception(f'tree depth is less than 3, please check your tree file{tree}')
    else:   
        return True

def fcheck(*files):
    for file in files:
       if not os.path.exists(file):
           raise FileNotFoundError(f"{file} not exist")


class Gu99():
    def __init__(self,aln_file,*tree_files,cluster_name=None,trees:list=[]) -> None:
        """__init__ 
        :param aln_file: aln_file path
        :param tree_files: tree_files path
        :param cluster_name: clusters name, defaults to None
        :param trees: list of Biopython.Phylo tree object, defaults to []
        """
        fcheck(aln_file,*tree_files)
        self.input = [aln_file]
        if trees!=[]:
            for tree in trees: 
                check_tree(tree)
                self.input.extend([tree.format("newick").strip("\n")])
        else:
            for tree_file in tree_files:
                self.input.extend([load_tree_file(tree_file,check=True)])
        
        if cluster_name is None:
            self.cluster_name = [f"cluster_{i}" for i in range(1,len(self.input))]
        elif isinstance(cluster_name,list):
            #Judgment cluster_name is a list
            self.cluster_name = cluster_name
        else:
            raise Exception("cluster_name must be a list")
        self._calculate()
        self._summary = self.summary()
    def _calculate(self):
        """
        Create a new Gu99 Calculator and Complete the calculation steps.
        """
        self.calculator = _gu99cpp.create_calculator(self.input,self.cluster_name)
        self.calculator.calculate()
    def summary(self):
        calculator = self.calculator
        columns_names = get_colnames(calculator._r_names())
        summary_results = pd.DataFrame(columns=columns_names)
        dicts = calculator._summary()
        for dict in dicts:
          summary_results.loc[dict["name"]] = dict["values"]
        return summary_results
    def results(self):
        calculator = self.calculator
        columns_names = get_colnames(calculator._r_names())
        results = pd.DataFrame(calculator._results(),columns=columns_names,index=[i+1 for i in calculator._kept()])
        return results
    def fundist(self):
        if len(self.input)<=3:
            return None
        else:
            n = len(self.input)-1  
            theta = np.zeros((n, n))
            k = 0
            for i in range(n):
                for j in range(i + 1, n):
                    theta[j, i] = theta[i, j] = self._summary.iloc[0,k]
                    k += 1
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
            columns_names = [f"cluster{i}" for i in range(1,len(self.input))]
            return pd.DataFrame(dist_results,index=columns_names)
class Gu2001(): 
    def __init__(self,aln_file,*tree_files,cluster_name=None,trees:list=[]) -> None:
        """__init__ 
        :param aln_file: aln_file path
        :param tree_files: tree_files path
        :param cluster_name: clusters name, defaults to None
        :param trees: list of Biopython.Phylo tree object, defaults to []
        """
        fcheck(aln_file,*tree_files)
        self.input = [aln_file]
        if trees!=[]:
            for tree in trees: 
                check_tree(tree)
                self.input.extend([tree.format("newick").strip("\n")])
        else:
            for tree_file in tree_files:
                self.input.extend([load_tree_file(tree_file,check=True)])
        if cluster_name is None:
            self.cluster_name = [f"cluster_{i}" for i in range(1,len(self.input))]
        elif isinstance(cluster_name,list):
            #Judgment cluster_name is a list
            self.cluster_name = cluster_name
        else:
            raise Exception("cluster_name must be a list")
        self._calculate()
        self._summary = self.summary()
    def _calculate(self):
        """
        Create a new Gu2001 Calculator and Complete the calculation steps.
        """
        self.calculator = _gu2001cpp.create_calculator(self.input,self.cluster_name)
        self.calculator.calculate()
    def _help(self):
        print(help(_gu2001cpp))
    def summary(self):
        calculator = self.calculator
        columns_names = get_colnames(calculator._r_names())
        summary_results = pd.DataFrame(columns=columns_names)
        dicts = calculator._summary()
        for dict in dicts:
          summary_results.loc[dict["name"]] = dict["values"]
        summary_results.index.name = "Parameters"
        return summary_results
    def results(self):
        calculator = self.calculator
        columns_names = get_colnames(calculator._r_names())
        results = pd.DataFrame(calculator._results(),columns=columns_names,index=[i+1 for i in calculator._kept()])
        results.index.name = "Position"
        return results 
class Type2():
    def __init__(self,aln_file,*tree_files,cluster_name=None,trees:list=[]) -> None:
        """__init__ 
        :param aln_file: aln_file path
        :param tree_files: tree_files path
        :param cluster_name: clusters name, defaults to None
        :param trees: list of Biopython.Phylo tree object, defaults to []
        """
        fcheck(aln_file,*tree_files)
        self.input = [aln_file]
        if trees!=[]:
            for tree in trees: 
                check_tree(tree)
                self.input.extend([tree.format("newick").strip("\n")])
        else:
            for tree_file in tree_files:
                self.input.extend([load_tree_file(tree_file,check=True)])
        if cluster_name is None:
            self.cluster_name = [f"cluster_{i}" for i in range(1,len(self.input))]
        elif isinstance(cluster_name,list):
            #Judgment cluster_name is a list
            self.cluster_name = cluster_name
        else:
            raise Exception("cluster_name must be a list")
        self._calculate()
        self._summary = self.summary()
    def _calculate(self):
        """
        Create a new type2 Calculator and Complete the calculation steps.
        """
        self.calculator = _type2cpp.create_calculator(self.input,self.cluster_name)
        self.calculator.calculate()
    def summary(self):
        calculator = self.calculator
        columns_names = get_colnames(calculator._r_names())
        summary_results = pd.DataFrame(columns=columns_names)
        dicts = calculator._summary()
        for dict in dicts:
          summary_results.loc[dict["name"]] = dict["values"]
        summary_results.index.name = "Parameters"
        return summary_results
    def results(self):
        calculator = self.calculator
        columns_names = get_colnames(calculator._r_names())
        results = pd.DataFrame(calculator._results(),columns=columns_names,index=[i+1 for i in calculator._kept()])
        results.index.name = "Position"
        return results
    def _help(self):
        print(help(_type2cpp))
class Asym():
    def __init__(self,aln_file,*tree_files,cluster_name=None,trees:list=[]) -> None:
        """__init__ 
        :param aln_file: aln_file path
        :param tree_files: tree_files path
        :param cluster_name: clusters name, defaults to None
        :param trees: list of Biopython.Phylo tree object, defaults to []
        """
        fcheck(aln_file,*tree_files)
        self.input = [aln_file]
        if trees!=[]:
            for tree in trees: 
                check_tree(tree)
                self.input.extend([tree.format("newick").strip("\n")])
        else:
            for tree_file in tree_files:
                self.input.extend([load_tree_file(tree_file,check=False)])
        if cluster_name is None:
            self.cluster_name = [f"cluster_{i}" for i in range(1,len(self.input))]
        elif isinstance(cluster_name,list):
            #Judgment cluster_name is a list
            self.cluster_name = cluster_name
        else:
            raise Exception("cluster_name must be a list")
        self._calculate()
    def _calculate(self):
        """
        Create a new asym Calculator and Complete the calculation steps.
        """
        self.calculator = _asymcpp.create_calculator(self.input,self.cluster_name)
        self.calculator.calculate()
    def results(self):
        calculator = self.calculator
        columns_names = get_colnames(calculator._r_names())
        results = pd.DataFrame(calculator._results(),columns=columns_names,index = list(range(1,np.size(calculator._results(),0)+1)))
        results.index.name = "Cluster Number of Outgroup" 
        return results
    def _help(self):
        print(help(_asymcpp))
class Effective():
    def __init__(self,aln_file,*tree_files,cluster_name=None,trees:list=[]) -> None:
        """__init__ 
        :param aln_file: aln_file path
        :param tree_files: tree_files path
        :param cluster_name: clusters name, defaults to None
        :param trees: list of Biopython.Phylo tree object, defaults to []
        """
        fcheck(aln_file,*tree_files)
        self.input = [aln_file]
        if trees!=[]:
            for tree in trees: 
                check_tree(tree)
                self.input.extend([tree.format("newick").strip("\n")])
        else:
            for tree_file in tree_files:
                self.input.extend([load_tree_file(tree_file,check=True)])
        if cluster_name is None:
            self.cluster_name = [f"cluster_{i}" for i in range(1,len(self.input))]
        elif isinstance(cluster_name,list):
            #Judgment cluster_name is a list
            self.cluster_name = cluster_name
        else:
            raise Exception("cluster_name must be a list")
        self._calculate()
        self.type1_effective_number = self.calculator._results1().shape[0]
        self.type2_effective_number = self.calculator._results2().shape[0]
        print(f"Type1 Effective Number of Sites is {self.type1_effective_number}, Type2 Effective Number of Sites is {self.type2_effective_number}")
    def _calculate(self):
        """
        Create a new Effective number of site Calculator and Complete the calculation steps.
        """
        self.calculator = _effectivecpp.create_calculator(self.input,self.cluster_name)
        self.calculator.calculate()
    def type1_results(self):
        calculator = self.calculator
        columns_names = get_colnames(calculator._r_names())
        results = pd.DataFrame(calculator._results1(),columns=columns_names)
        results.index.name = "Number"
        return results
    def type2_results(self):
        calculator = self.calculator
        columns_names = get_colnames(calculator._r_names())
        results = pd.DataFrame(calculator._results2(),columns=columns_names)
        results.index.name = "Number"
        return results
    def _help(self):
        print(help(_effectivecpp))
class Fdr():
    def __init__(self,aln_file,*tree_files,cluster_name=None,trees:list=[]) -> None:
        """__init__ 
        :param aln_file: aln_file path
        :param tree_files: tree_files path
        :param cluster_name: clusters name, defaults to None
        :param trees: list of Biopython.Phylo tree object, defaults to []
        """
        fcheck(aln_file,*tree_files)
        self.input = [aln_file]
        if trees!=[]:
            for tree in trees: 
                check_tree(tree)
                self.input.extend([tree.format("newick").strip("\n")])
        else:
            for tree_file in tree_files:
                self.input.extend([load_tree_file(tree_file,check=True)])
        if cluster_name is None:
            self.cluster_name = [f"cluster_{i}" for i in range(1,len(self.input))]
        elif isinstance(cluster_name,list):
            #Judgment cluster_name is a list
            self.cluster_name = cluster_name
        else:
            raise Exception("cluster_name must be a list")
        self._calculate()
    def _calculate(self):
        """
        Create a new fdr Calculator and Complete the calculation steps.
        """
        self.calculator = _fdrcpp.create_calculator(self.input,self.cluster_name)
        self.calculator.calculate()
    def type1_results(self):
        calculator = self.calculator
        columns_names = get_colnames(calculator._r_names())
        results = pd.DataFrame(calculator._results1(),columns=columns_names)
        results = results.set_index(results.columns[0],drop=True)
        return results
    def type2_results(self):
        calculator = self.calculator
        columns_names = get_colnames(calculator._r_names())
        results = pd.DataFrame(calculator._results2(),columns=columns_names)
        results = results.set_index(results.columns[0],drop=True)
        return results
    def _help(self):
        print(help(_fdrcpp))     
class Rvs():
    def __init__(self,aln_file,*tree_files,cluster_name=None,trees:list=[]) -> None:
        """__init__ 
        :param aln_file: aln_file path
        :param tree_files: tree_files path
        :param cluster_name: clusters name, defaults to None
        :param trees: list of Biopython.Phylo tree object, defaults to []
        """
        fcheck(aln_file,*tree_files)
        self.input = [aln_file]
        if trees!=[]:
            for tree in trees: 
                check_tree(tree)
                self.input.extend([tree.format("newick").strip("\n")])
        else:
            for tree_file in tree_files:
                self.input.extend([load_tree_file(tree_file,check=True)])
        if cluster_name is None:
            self.cluster_name = [f"cluster_{i}" for i in range(1,len(self.input))]
        elif isinstance(cluster_name,list):
            #Judgment cluster_name is a list
            self.cluster_name = cluster_name
        else:
            raise Exception("cluster_name must be a list")
        self._calculate()
    def _calculate(self):
        """
        Create a new rvs Calculator and Complete the calculation steps.
        """
        self.calculator = _rvscpp.create_calculator(self.input,self.cluster_name)
        self.calculator.calculate()
    def summary(self):
        calculator = self.calculator
        columns_names = get_colnames(calculator._s_names())
        summary_results = pd.DataFrame(columns=columns_names)
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
        columns_names = get_colnames(calculator._r_names())
        results = pd.DataFrame(calculator._results(),columns=columns_names,index=[i+1 for i in calculator._kept()])
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
    def __init__(self,aln_file,*tree_files,cluster_name=None,trees=[]) -> None:
        """__init__ 
        :param aln_file: aln_file path
        :param tree_files: tree_files path
        :param cluster_name: clusters name, defaults to None
        :param trees: list of Biopython.Phylo tree object, defaults to []
        """
        fcheck(aln_file,*tree_files)
        self.input = [aln_file]
        if trees!=[]:
            for tree in trees: 
                check_tree(tree)
                self.input.extend([tree.format("newick").strip("\n")])
        else:
            for tree_file in tree_files:
                self.input.extend([load_tree_file(tree_file,check=True)])
        if cluster_name is None:
            self.cluster_name = [f"cluster_{i}" for i in range(1,len(self.input))]
        elif isinstance(cluster_name,list):
            #Judgment cluster_name is a list
            self.cluster_name = cluster_name
        else:
            raise Exception("cluster_name must be a list")
        self._calculate()
    def _calculate(self):
        """
        Create a new TypeOneAnalysis Calculator and Complete the calculation steps.
        """
        self.calculator = _typeOneAnalysiscpp.create_calculator(self.input,self.cluster_name)
        self.calculator.calculate()
    def summary(self):
        calculator = self.calculator
        columns_names = get_colnames(calculator._s_names())
        summary_results = pd.DataFrame(columns=columns_names)
        dicts = calculator._summary()
        for dict in dicts:
          summary_results.loc[dict["name"]] = dict["values"]
        summary_results.index.name = "Parameters"
        return summary_results
    def results(self):
        calculator = self.calculator
        columns_names = get_colnames(calculator._r_names())
        results = pd.DataFrame(calculator._results(),columns=columns_names,index=[i+1 for i in calculator._kept()])
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