from typing import List, Optional
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Phylo
from Bio.Phylo.BaseTree import Tree

from . import _gu2001cpp, _gu99cpp, _type2cpp, _asymcpp, _fdrcpp, _effectivecpp, _rvscpp, _typeOneAnalysiscpp


def read_tree(tree_file: str) -> Tree:
    """Read a tree file in newick or nexus format."""
    try:
        return Phylo.read(tree_file, 'newick')
    except ValueError:
        try:
            return Phylo.read(tree_file, 'nexus')
        except ValueError:
            raise ValueError(f"Tree file '{tree_file}' is not in newick or nexus format.")

def get_colnames(r_names: List[str]) -> List[str]:
    """Get column names from r_names."""
    return r_names.copy()

def load_tree_file(tree_file: str, check: bool = True) -> str:
    """Load a tree file and optionally check its validity."""
    tree = read_tree(tree_file)
    if check:
        check_tree(tree)
    return tree.format("newick").strip("\n")

def check_tree_file(*tree_files: str) -> bool:
    """Check the validity of multiple tree files."""
    for tree_file in tree_files:
        tree = Phylo.read(tree_file, 'newick')
        tree_depth = max(len(tree.trace(tree.root, clade)) for clade in tree.get_terminals())
        if tree_depth < 3:
            raise ValueError(f'Tree depth is less than 3, please check your tree file: {tree_file}')
    return True

def check_tree(tree: Tree) -> bool:
    """Check the validity of a single tree."""
    tree_depth = max(len(tree.trace(tree.root, clade)) for clade in tree.get_terminals())
    if tree_depth < 3:
        raise ValueError(f'Tree depth is less than 3, please check your tree: {tree}')
    return True

def fcheck(*files: str) -> None:
    """Check if files exist."""
    for file in files:
        if not os.path.exists(file):
            raise FileNotFoundError(f"{file} does not exist")

class Gu99:
    def __init__(
        self, 
        aln_file: str, 
        *tree_files: str, 
        cluster_name: Optional[List[str]] = None, 
        trees: List[Tree] = []
    ) -> None:
        """
        Initialize Gu99 analysis.

        Args:
            aln_file: Path to alignment file.
            tree_files: Paths to tree files.
            cluster_name: Names of clusters.
            trees: List of Biopython.Phylo tree objects.
        """
        fcheck(aln_file, *tree_files)
        self.input = [aln_file]
        
        if trees:
            for tree in trees:
                check_tree(tree)
                self.input.append(tree.format("newick").strip("\n"))
        else:
            for tree_file in tree_files:
                self.input.append(load_tree_file(tree_file, check=True))
        
        if cluster_name is None:
            self.cluster_name = [f"cluster_{i}" for i in range(1, len(self.input))]
        elif isinstance(cluster_name, list):
            self.cluster_name = cluster_name
        else:
            raise TypeError("cluster_name must be a list")
        
        self._calculate()
        self._summary = self.summary()

    def _calculate(self) -> None:
        """Create a new Gu99 Calculator and complete the calculation steps."""
        self.calculator = _gu99cpp.create_calculator(self.input, self.cluster_name)
        self.calculator.calculate()

    def summary(self) -> pd.DataFrame:
        """Generate a summary of the results."""
        columns_names = get_colnames(self.calculator._r_names())
        summary_results = pd.DataFrame(columns=columns_names)
        for dict_item in self.calculator._summary():
            summary_results.loc[dict_item["name"]] = dict_item["values"]
        return summary_results

    def results(self) -> pd.DataFrame:
        """Generate detailed results."""
        columns_names = get_colnames(self.calculator._r_names())
        results = pd.DataFrame(
            self.calculator._results(),
            columns=columns_names,
            index=[i+1 for i in self.calculator._kept()]
        )
        return results

    def fundist(self) -> Optional[pd.DataFrame]:
        """Calculate functional distance if applicable."""
        if len(self.input) <= 3:
            return None
        
        n = len(self.input) - 1
        theta = np.zeros((n, n))
        k = 0
        for i in range(n):
            for j in range(i + 1, n):
                theta[j, i] = theta[i, j] = self._summary.iloc[0, k]
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
        
        columns_names = [f"cluster{i}" for i in range(1, len(self.input))]
        return pd.DataFrame(dist_results, index=columns_names)

class Gu2001:
    def __init__(
        self,
        aln_file: str,
        *tree_files: str,
        cluster_name: Optional[List[str]] = None,
        trees: List[Tree] = []
    ) -> None:
        """
        Initialize Gu2001 analysis.

        Args:
            aln_file: Path to alignment file.
            tree_files: Paths to tree files.
            cluster_name: Names of clusters.
            trees: List of Biopython.Phylo tree objects.
        """
        fcheck(aln_file, *tree_files)
        self.input = [aln_file]
        
        if trees:
            for tree in trees:
                check_tree(tree)
                self.input.append(tree.format("newick").strip("\n"))
        else:
            for tree_file in tree_files:
                self.input.append(load_tree_file(tree_file, check=True))
        
        if cluster_name is None:
            self.cluster_name = [f"cluster_{i}" for i in range(1, len(self.input))]
        elif isinstance(cluster_name, list):
            self.cluster_name = cluster_name
        else:
            raise TypeError("cluster_name must be a list")
        
        self._calculate()
        self._summary = self.summary()

    def _calculate(self) -> None:
        """Create a new Gu2001 Calculator and complete the calculation steps."""
        self.calculator = _gu2001cpp.create_calculator(self.input, self.cluster_name)
        self.calculator.calculate()

    def _help(self) -> None:
        """Print help information for _gu2001cpp."""
        print(help(_gu2001cpp))

    def summary(self) -> pd.DataFrame:
        """Generate a summary of the results."""
        columns_names = get_colnames(self.calculator._r_names())
        summary_results = pd.DataFrame(columns=columns_names)
        for dict_item in self.calculator._summary():
            summary_results.loc[dict_item["name"]] = dict_item["values"]
        summary_results.index.name = "Parameters"
        return summary_results

    def results(self) -> pd.DataFrame:
        """Generate detailed results."""
        columns_names = get_colnames(self.calculator._r_names())
        results = pd.DataFrame(
            self.calculator._results(),
            columns=columns_names,
            index=[i+1 for i in self.calculator._kept()]
        )
        results.index.name = "Position"
        return results

class Type2:
    def __init__(
        self,
        aln_file: str,
        *tree_files: str,
        cluster_name: Optional[List[str]] = None,
        trees: List[Tree] = []
    ) -> None:
        """
        Initialize Type2 analysis.

        Args:
            aln_file: Path to alignment file.
            tree_files: Paths to tree files.
            cluster_name: Names of clusters.
            trees: List of Biopython.Phylo tree objects.
        """
        fcheck(aln_file, *tree_files)
        self.input = [aln_file]
        
        if trees:
            for tree in trees:
                check_tree(tree)
                self.input.append(tree.format("newick").strip("\n"))
        else:
            for tree_file in tree_files:
                self.input.append(load_tree_file(tree_file, check=True))
        
        if cluster_name is None:
            self.cluster_name = [f"cluster_{i}" for i in range(1, len(self.input))]
        elif isinstance(cluster_name, list):
            self.cluster_name = cluster_name
        else:
            raise TypeError("cluster_name must be a list")
        
        self._calculate()
        self._summary = self.summary()

    def _calculate(self) -> None:
        """Create a new Type2 Calculator and complete the calculation steps."""
        self.calculator = _type2cpp.create_calculator(self.input, self.cluster_name)
        self.calculator.calculate()

    def summary(self) -> pd.DataFrame:
        """Generate a summary of the results."""
        columns_names = get_colnames(self.calculator._r_names())
        summary_results = pd.DataFrame(columns=columns_names)
        for dict_item in self.calculator._summary():
            summary_results.loc[dict_item["name"]] = dict_item["values"]
        summary_results.index.name = "Parameters"
        return summary_results

    def results(self) -> pd.DataFrame:
        """Generate detailed results."""
        columns_names = get_colnames(self.calculator._r_names())
        results = pd.DataFrame(
            self.calculator._results(),
            columns=columns_names,
            index=[i+1 for i in self.calculator._kept()]
        )
        results.index.name = "Position"
        return results

    def _help(self) -> None:
        """Print help information for _type2cpp."""
        print(help(_type2cpp))

class Asym:
    def __init__(
        self,
        aln_file: str,
        *tree_files: str,
        cluster_name: Optional[List[str]] = None,
        trees: List[Tree] = []
    ) -> None:
        """
        Initialize Asym analysis.

        Args:
            aln_file: Path to alignment file.
            tree_files: Paths to tree files.
            cluster_name: Names of clusters.
            trees: List of Biopython.Phylo tree objects.
        """
        fcheck(aln_file, *tree_files)
        self.input = [aln_file]
        
        if trees:
            for tree in trees:
                check_tree(tree)
                self.input.append(tree.format("newick").strip("\n"))
        else:
            for tree_file in tree_files:
                self.input.append(load_tree_file(tree_file, check=False))
        
        if cluster_name is None:
            self.cluster_name = [f"cluster_{i}" for i in range(1, len(self.input))]
        elif isinstance(cluster_name, list):
            self.cluster_name = cluster_name
        else:
            raise TypeError("cluster_name must be a list")
        
        self._calculate()

    def _calculate(self) -> None:
        """Create a new Asym Calculator and complete the calculation steps."""
        self.calculator = _asymcpp.create_calculator(self.input, self.cluster_name)
        self.calculator.calculate()

    def results(self) -> pd.DataFrame:
        """Generate detailed results."""
        columns_names = get_colnames(self.calculator._r_names())
        results = pd.DataFrame(
            self.calculator._results(),
            columns=columns_names,
            index=range(1, np.size(self.calculator._results(), 0) + 1)
        )
        results.index.name = "Cluster Number of Outgroup"
        return results

    def _help(self) -> None:
        """Print help information for _asymcpp."""
        print(help(_asymcpp))

class Effective:
    def __init__(
        self,
        aln_file: str,
        *tree_files: str,
        cluster_name: Optional[List[str]] = None,
        trees: List[Tree] = []
    ) -> None:
        """
        Initialize Effective analysis.

        Args:
            aln_file: Path to alignment file.
            tree_files: Paths to tree files.
            cluster_name: Names of clusters.
            trees: List of Biopython.Phylo tree objects.
        """
        fcheck(aln_file, *tree_files)
        self.input = [aln_file]
        
        if trees:
            for tree in trees:
                check_tree(tree)
                self.input.append(tree.format("newick").strip("\n"))
        else:
            for tree_file in tree_files:
                self.input.append(load_tree_file(tree_file, check=True))
        
        if cluster_name is None:
            self.cluster_name = [f"cluster_{i}" for i in range(1, len(self.input))]
        elif isinstance(cluster_name, list):
            self.cluster_name = cluster_name
        else:
            raise TypeError("cluster_name must be a list")
        
        self._calculate()
        self.type1_effective_number = self.calculator._results1().shape[0]
        self.type2_effective_number = self.calculator._results2().shape[0]
        print(f"Type1 Effective Number of Sites is {self.type1_effective_number}, Type2 Effective Number of Sites is {self.type2_effective_number}")

    def _calculate(self) -> None:
        """Create a new Effective number of site Calculator and complete the calculation steps."""
        self.calculator = _effectivecpp.create_calculator(self.input, self.cluster_name)
        self.calculator.calculate()

    def type1_results(self) -> pd.DataFrame:
        """Generate Type 1 results."""
        columns_names = get_colnames(self.calculator._r_names())
        results = pd.DataFrame(self.calculator._results1(), columns=columns_names)
        results.index.name = "Number"
        return results

    def type2_results(self) -> pd.DataFrame:
        """Generate Type 2 results."""
        columns_names = get_colnames(self.calculator._r_names())
        results = pd.DataFrame(self.calculator._results2(), columns=columns_names)
        results.index.name = "Number"
        return results

    def _help(self) -> None:
        """Print help information for _effectivecpp."""
        print(help(_effectivecpp))

# Continue with other classes...

class Fdr:
    def __init__(
        self,
        aln_file: str,
        *tree_files: str,
        cluster_name: Optional[List[str]] = None,
        trees: List[Tree] = []
    ) -> None:
        """
        Initialize Fdr analysis.

        Args:
            aln_file: Path to alignment file.
            tree_files: Paths to tree files.
            cluster_name: Names of clusters.
            trees: List of Biopython.Phylo tree objects.
        """
        fcheck(aln_file, *tree_files)
        self.input = [aln_file]
        
        if trees:
            for tree in trees:
                check_tree(tree)
                self.input.append(tree.format("newick").strip("\n"))
        else:
            for tree_file in tree_files:
                self.input.append(load_tree_file(tree_file, check=True))
        
        if cluster_name is None:
            self.cluster_name = [f"cluster_{i}" for i in range(1, len(self.input))]
        elif isinstance(cluster_name, list):
            self.cluster_name = cluster_name
        else:
            raise TypeError("cluster_name must be a list")
        
        self._calculate()

    def _calculate(self) -> None:
        """Create a new fdr Calculator and complete the calculation steps."""
        self.calculator = _fdrcpp.create_calculator(self.input, self.cluster_name)
        self.calculator.calculate()

    def type1_results(self) -> pd.DataFrame:
        """Generate Type 1 results."""
        columns_names = get_colnames(self.calculator._r_names())
        results = pd.DataFrame(self.calculator._results1(), columns=columns_names)
        return results.set_index(results.columns[0], drop=True)

    def type2_results(self) -> pd.DataFrame:
        """Generate Type 2 results."""
        columns_names = get_colnames(self.calculator._r_names())
        results = pd.DataFrame(self.calculator._results2(), columns=columns_names)
        return results.set_index(results.columns[0], drop=True)

    def _help(self) -> None:
        """Print help information for _fdrcpp."""
        print(help(_fdrcpp))

class Rvs:
    def __init__(
        self,
        aln_file: str,
        *tree_files: str,
        cluster_name: Optional[List[str]] = None,
        trees: List[Tree] = []
    ) -> None:
        """
        Initialize Rvs analysis.

        Args:
            aln_file: Path to alignment file.
            tree_files: Paths to tree files.
            cluster_name: Names of clusters.
            trees: List of Biopython.Phylo tree objects.
        """
        fcheck(aln_file, *tree_files)
        self.input = [aln_file]
        
        if trees:
            for tree in trees:
                check_tree(tree)
                self.input.append(tree.format("newick").strip("\n"))
        else:
            for tree_file in tree_files:
                self.input.append(load_tree_file(tree_file, check=True))
        
        if cluster_name is None:
            self.cluster_name = [f"cluster_{i}" for i in range(1, len(self.input))]
        elif isinstance(cluster_name, list):
            self.cluster_name = cluster_name
        else:
            raise TypeError("cluster_name must be a list")
        
        self._calculate()

    def _calculate(self) -> None:
        """Create a new rvs Calculator and complete the calculation steps."""
        self.calculator = _rvscpp.create_calculator(self.input, self.cluster_name)
        self.calculator.calculate()

    def summary(self) -> pd.DataFrame:
        """Generate a summary of the results."""
        columns_names = get_colnames(self.calculator._s_names())
        summary_results = pd.DataFrame(columns=columns_names)
        for dict_item in self.calculator._summary():
            summary_results.loc[dict_item["name"]] = dict_item["values"]
        summary_results.index.name = "Parameters"
        return summary_results

    def results(self) -> pd.DataFrame:
        """
        Generate detailed results.
        
        Parameters:
        - Xk: Number of Changes
        - Rk: Posterior Mean of Evolutionary Rate
        """
        columns_names = get_colnames(self.calculator._r_names())
        results = pd.DataFrame(
            self.calculator._results(),
            columns=columns_names,
            index=[i+1 for i in self.calculator._kept()]
        )
        return results

    def _help(self) -> None:
        """Print help information for _rvscpp."""
        print(help(_rvscpp))

class TypeOneAnalysis:
    """
    A class for type one analysis.
    
    Note: Under the two-state model (functional divergence unrelated F0 or related F1),
    there are eight possible combined states for three duplicate clusters,
    which can be reduced to five nondegenerate patterns.
    S0=(F0, F0, F0) means no type-one divergence occurred in any clusters.
    S1=(F1, F0, F0) means type-one functional divergence occurred only in cluster 1,
    similarly: S2=(F0, F1, F0) and S3=(F0, F0, F1).
    The final pattern S4 is for the rest of four states, each of which has two or three
    clusters that have experienced type-one functional divergence.
    """

    def __init__(
        self,
        aln_file: str,
        *tree_files: str,
        cluster_name: Optional[List[str]] = None,
        trees: List[Tree] = []
    ) -> None:
        """
        Initialize TypeOneAnalysis.

        Args:
            aln_file: Path to alignment file.
            tree_files: Paths to tree files.
            cluster_name: Names of clusters.
            trees: List of Biopython.Phylo tree objects.
        """
        fcheck(aln_file, *tree_files)
        self.input = [aln_file]
        
        if trees:
            for tree in trees:
                check_tree(tree)
                self.input.append(tree.format("newick").strip("\n"))
        else:
            for tree_file in tree_files:
                self.input.append(load_tree_file(tree_file, check=True))
        
        if cluster_name is None:
            self.cluster_name = [f"cluster_{i}" for i in range(1, len(self.input))]
        elif isinstance(cluster_name, list):
            self.cluster_name = cluster_name
        else:
            raise TypeError("cluster_name must be a list")
        
        self._calculate()

    def _calculate(self) -> None:
        """Create a new TypeOneAnalysis Calculator and complete the calculation steps."""
        self.calculator = _typeOneAnalysiscpp.create_calculator(self.input, self.cluster_name)
        self.calculator.calculate()

    def summary(self) -> pd.DataFrame:
        """Generate a summary of the results."""
        columns_names = get_colnames(self.calculator._s_names())
        summary_results = pd.DataFrame(columns=columns_names)
        for dict_item in self.calculator._summary():
            summary_results.loc[dict_item["name"]] = dict_item["values"]
        summary_results.index.name = "Parameters"
        return summary_results

    def results(self) -> pd.DataFrame:
        """Generate detailed results."""
        columns_names = get_colnames(self.calculator._r_names())
        results = pd.DataFrame(
            self.calculator._results(),
            columns=columns_names,
            index=[i+1 for i in self.calculator._kept()]
        )
        results.index.name = "Position"
        return results

    def _help(self) -> None:
        """Print help information for _typeOneAnalysiscpp."""
        print(help(_typeOneAnalysiscpp))

# End of file

if __name__ == "__main__":
    # test gu99
    T = Gu99("../test_data/CASP.aln", "../test_data/cl1.tree", "../test_data/cl2.tree","../test_data/cl3.tree")
    summary = T.summary()
    results = T.results()
    print(summary)
    print("------------")
    print(results)
    # test gu2001
    # T = Gu2001("../test_data/CASP.aln", "../test_data/treelength/cl1.tree", "../test_data/treelength/cl2.tree","../test_data/treelength/cl3.tree")
    # summary = T.summary()
    # results = T.results()
    # print(summary)
    # print("------------")
    # print(results)
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