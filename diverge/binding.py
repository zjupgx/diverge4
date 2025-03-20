from typing import List, Optional, Dict, Any, Tuple, Union, Callable
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Phylo
from Bio.Phylo.BaseTree import Tree

from . import _gu2001cpp, _gu99cpp, _type2cpp, _asymcpp, _fdrcpp, _effectivecpp, _rvscpp, _typeOneAnalysiscpp


def read_tree(tree_file: str) -> Tree:
    """Read a tree file in newick or nexus format."""
    if not os.path.exists(tree_file):
        raise FileNotFoundError(f"Tree file not found: {tree_file}")
        
    formats = ['newick', 'nexus']
    last_error = None
    
    for fmt in formats:
        try:
            return Phylo.read(tree_file, fmt)
        except Exception as e:
            last_error = e
    
    # 如果所有格式都失败
    raise ValueError(f"Tree file '{tree_file}' could not be parsed: {last_error}")

def get_colnames(r_names: List[str]) -> List[str]:
    """Get column names from r_names."""
    return r_names.copy()

def load_tree_file(tree_file: str, check: bool = True) -> str:
    """Load a tree file and optionally check its validity."""
    tree = read_tree(tree_file)
    if check:
        check_tree(tree)
    return tree.format("newick").strip("\n")

def check_tree_depth(tree: Tree, min_depth: int = 3) -> bool:
    """Check if tree has sufficient depth."""
    tree_depth = max(len(tree.trace(tree.root, clade)) for clade in tree.get_terminals())
    return tree_depth >= min_depth

def check_tree(tree: Tree) -> bool:
    """Check the validity of a single tree."""
    if not check_tree_depth(tree):
        raise ValueError(f'Tree depth is less than 3, please check your tree: {tree}')
    return True

def check_tree_file(*tree_files: str) -> bool:
    """Check the validity of multiple tree files."""
    for tree_file in tree_files:
        tree = read_tree(tree_file)
        if not check_tree_depth(tree):
            raise ValueError(f'Tree depth is less than 3, please check your tree file: {tree_file}')
    return True

def fcheck(*files: str) -> None:
    """Check if files exist."""
    for file in files:
        if not os.path.exists(file):
            raise FileNotFoundError(f"{file} does not exist")

class BaseAnalysis:
    """Base class for all analysis types."""
    
    def __init__(
        self,
        aln_file: str,
        *tree_files: str,
        cluster_name: Optional[List[str]] = None,
        trees: List[Tree] = [],
        check_trees: bool = True,
        calculator_module: Any = None
    ) -> None:
        """
        Initialize analysis with alignment and tree files.

        Args:
            aln_file: Path to alignment file.
            tree_files: Paths to tree files.
            cluster_name: Names of clusters.
            trees: List of Biopython.Phylo tree objects.
            check_trees: Whether to check tree validity.
            calculator_module: C++ calculator module to use.
        """
        # 检查文件
        fcheck(aln_file, *tree_files)
        
        # 设置计算模块
        self.calculator_module = calculator_module
        
        # 准备输入
        self.input = [aln_file]
        self._prepare_trees(tree_files, trees, check_trees)
        self._prepare_cluster_names(cluster_name)
        
        # 执行计算
        self._calculate()
        
        # 缓存结果以提高性能
        self._summary_cache = None
        self._results_cache = None
    
    def _prepare_trees(self, tree_files: Tuple[str, ...], trees: List[Tree], check_trees: bool) -> None:
        """Prepare tree inputs."""
        if trees:
            for tree in trees:
                if check_trees:
                    check_tree(tree)
                self.input.append(tree.format("newick").strip("\n"))
        else:
            for tree_file in tree_files:
                self.input.append(load_tree_file(tree_file, check=check_trees))
    
    def _prepare_cluster_names(self, cluster_name: Optional[List[str]]) -> None:
        """Set up cluster names."""
        if cluster_name is None:
            self.cluster_name = [f"cluster_{i}" for i in range(1, len(self.input))]
        elif isinstance(cluster_name, list):
            self.cluster_name = cluster_name
        else:
            raise TypeError("cluster_name must be a list")
    
    def _calculate(self) -> None:
        """Create calculator and perform calculations. To be implemented by subclasses."""
        raise NotImplementedError("Subclasses must implement _calculate method")
    
    def _help(self) -> None:
        """Print help information for the calculator module."""
        if hasattr(self, 'calculator_module'):
            print(help(self.calculator_module))
    
    def summary(self) -> pd.DataFrame:
        """Generate summary of results with caching."""
        if self._summary_cache is not None:
            return self._summary_cache
            
        if not hasattr(self.calculator, '_summary'):
            return None
            
        name_method = '_r_names' if hasattr(self.calculator, '_r_names') else '_s_names'
        columns_names = get_colnames(getattr(self.calculator, name_method)())
        
        summary_results = pd.DataFrame(columns=columns_names)
        for dict_item in self.calculator._summary():
            summary_results.loc[dict_item["name"]] = dict_item["values"]
        summary_results.index.name = "Parameters"
        
        self._summary_cache = summary_results
        return summary_results
    
    def results(self) -> pd.DataFrame:
        """Generate detailed results with caching."""
        if self._results_cache is not None:
            return self._results_cache
            
        if not hasattr(self.calculator, '_results'):
            return None
            
        columns_names = get_colnames(self.calculator._r_names())
        
        if hasattr(self.calculator, '_kept'):
            index = [i+1 for i in self.calculator._kept()]
        else:
            index = range(1, np.size(self.calculator._results(), 0) + 1)
            
        results = pd.DataFrame(
            self.calculator._results(),
            columns=columns_names,
            index=index
        )
        results.index.name = self.get_position_label()
        
        self._results_cache = results
        return results
    
    def get_position_label(self) -> str:
        """Get the label for the position index."""
        return "Position"
    
    def clear_cache(self) -> None:
        """Clear cached results."""
        self._summary_cache = None
        self._results_cache = None
    
    @classmethod
    def from_trees(cls, aln_file: str, trees: List[Tree], **kwargs):
        """Factory method to create analysis from Tree objects."""
        return cls(aln_file, trees=trees, **kwargs)
    
    @classmethod
    def from_files(cls, aln_file: str, *tree_files, **kwargs):
        """Factory method to create analysis from tree files."""
        return cls(aln_file, *tree_files, **kwargs)

    def _set_index_name(self, df: pd.DataFrame, name: str) -> pd.DataFrame:
        """Set index name for dataframe."""
        if df is not None and df.index.name is None:
            df.index.name = name
        return df

class Gu99(BaseAnalysis):
    """Class for Gu99 analysis."""
    
    def __init__(
        self, 
        aln_file: str, 
        *tree_files: str, 
        cluster_name: Optional[List[str]] = None, 
        trees: List[Tree] = []
    ) -> None:
        """Initialize Gu99 analysis."""
        super().__init__(aln_file, *tree_files, cluster_name=cluster_name, trees=trees, calculator_module=_gu99cpp)
        self._summary_result = self.summary()
        
    def summary(self) -> pd.DataFrame:
        """Generate a summary of the results."""
        columns_names = get_colnames(self.calculator._r_names())
        summary_results = pd.DataFrame(columns=columns_names)
        for dict_item in self.calculator._summary():
            summary_results.loc[dict_item["name"]] = dict_item["values"]
        return summary_results

    def _calculate(self) -> None:
        """Create a new Gu99 Calculator and complete the calculation steps."""
        self.calculator = self.calculator_module.create_calculator(self.input, self.cluster_name)
        self.calculator.calculate()
    
    def fundist(self) -> Optional[pd.DataFrame]:
        """Calculate functional distance if applicable."""
        if len(self.input) <= 3:
            return None
        
        summary = self._summary_result
        
        n = len(self.input) - 1
        theta = np.zeros((n, n))
        k = 0
        for i in range(n):
            for j in range(i + 1, n):
                theta[j, i] = theta[i, j] = summary.iloc[0, k]
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

    def plot_distance(self, figsize: Tuple[int, int] = (10, 6), 
                      title: str = "Functional Distance Between Clusters") -> plt.Figure:
        """
        Plot functional distances.
        
        Args:
            figsize: Figure size tuple (width, height)
            title: Plot title
            
        Returns:
            Matplotlib figure object
        """
        dist_df = self.fundist()
        if dist_df is None:
            raise ValueError("Cannot plot distances with less than 3 clusters")
            
        fig, ax = plt.subplots(figsize=figsize)
        dist_df.plot(kind='bar', ax=ax)
        ax.set_title(title)
        ax.set_ylabel("Functional Distance")
        ax.set_xlabel("Cluster")
        plt.tight_layout()
        
        return fig

class Gu2001(BaseAnalysis):
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
        super().__init__(aln_file, *tree_files, cluster_name=cluster_name, trees=trees)

    def _calculate(self) -> None:
        """Create a new Gu2001 Calculator and complete the calculation steps."""
        self.calculator = _gu2001cpp.create_calculator(self.input, self.cluster_name)
        self.calculator.calculate()

    def _help(self) -> None:
        """Print help information for _gu2001cpp."""
        print(help(_gu2001cpp))

class Type2(BaseAnalysis):
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
        super().__init__(aln_file, *tree_files, cluster_name=cluster_name, trees=trees)

    def _calculate(self) -> None:
        """Create a new Type2 Calculator and complete the calculation steps."""
        self.calculator = _type2cpp.create_calculator(self.input, self.cluster_name)
        self.calculator.calculate()

    def _help(self) -> None:
        """Print help information for _type2cpp."""
        print(help(_type2cpp))

class Asym(BaseAnalysis):
    """Class for Asym analysis."""
    
    def __init__(
        self,
        aln_file: str,
        *tree_files: str,
        cluster_name: Optional[List[str]] = None,
        trees: List[Tree] = []
    ) -> None:
        """Initialize Asym analysis."""
        self.calculator_module = _asymcpp
        super().__init__(aln_file, *tree_files, cluster_name=cluster_name, trees=trees, check_trees=False)

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

class Effective(BaseAnalysis):
    """Class for Effective analysis."""
    
    def __init__(
        self,
        aln_file: str,
        *tree_files: str,
        cluster_name: Optional[List[str]] = None,
        trees: List[Tree] = []
    ) -> None:
        """Initialize Effective analysis."""
        self.calculator_module = _effectivecpp
        super().__init__(aln_file, *tree_files, cluster_name=cluster_name, trees=trees)
        
        # 计算并显示有效位点数量
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

class Fdr(BaseAnalysis):
    """Class for Fdr analysis."""
    
    def __init__(
        self,
        aln_file: str,
        *tree_files: str,
        cluster_name: Optional[List[str]] = None,
        trees: List[Tree] = []
    ) -> None:
        """Initialize Fdr analysis."""
        self.calculator_module = _fdrcpp
        super().__init__(aln_file, *tree_files, cluster_name=cluster_name, trees=trees)

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

class Rvs(BaseAnalysis):
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
        super().__init__(aln_file, *tree_files, cluster_name=cluster_name, trees=trees)
    
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

class TypeOneAnalysis(BaseAnalysis):
    """
    A class for type one analysis.
    
    Note: Under the two-state model (functional divergence unrelated F0 or related F1),
    there are eight possible combined states for three duplicate clusters,
    which can be reduced to five nondegenerate patterns.
    S0=(F0, F0, F0) means no type-one divergence occurred in any clusters.
    S1=(F1, F0, F0) means type-one functional divergence occurred only in cluster 1,
    similarly: S2=(F0, F1, F0) and S3=(F0, F0, F1).
    The final pattern S4 is for the rest of four states, each of which has two or three
    clusters that have experienced type-one fCursorSessionTokenunctional divergence.
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
        super().__init__(aln_file, *tree_files, cluster_name=cluster_name, trees=trees)

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