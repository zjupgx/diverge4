from typing import List, Optional, Dict, Any, Tuple, Generator, Union
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Phylo, AlignIO
from Bio.Phylo import BaseTree
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from scipy.cluster.hierarchy import linkage, cut_tree, dendrogram, fcluster
from scipy.spatial.distance import squareform
from pymsaviz import MsaViz
from tqdm import tqdm
import time
import contextlib
import joblib
from tqdm.autonotebook import tqdm
import copy

from .binding import (
    Gu99, Gu2001, Type2, Fdr, Effective, Asym, TypeOneAnalysis, Rvs
)

class CalPipe:
    """
    CalPipe class: Wrap some calculate pipeline of diverge.
    """

    def __init__(self, aln_file: str, *tree_files: str) -> None:
        """
        Initialize CalPipe.

        Args:
            aln_file (str): Alignment file path
            *tree_files (str): Tree file paths
        """
        self._validate_files(aln_file, *tree_files)
        self.func_dict = {
            'Gu99': True, 'Gu2001': False, 'Type2': True, 'Fdr': True,
            'Effective': False, 'FunDist': False, 'Type1Analysis': False, 'Asym': False
        }
        self.tree_files = tree_files
        self.aln_file = aln_file
        self._pipe_select(*self.tree_files)
        self._pipeline()
        self.result_summary = self._result_summary
        self.detail = [key for key, value in self.func_dict.items() if value]

    def _validate_files(self, aln_file: str, *tree_files: str) -> None:
        """
        Validate existence of alignment and tree files.

        Args:
            aln_file (str): Alignment file path
            *tree_files (str): Tree file paths

        Raises:
            ValueError: If any of the files does not exist
        """
        if not os.path.isfile(aln_file):
            raise ValueError(f"Alignment file '{aln_file}' does not exist.")
        for tree_file in tree_files:
            if not os.path.isfile(tree_file):
                raise ValueError(f"Tree file '{tree_file}' does not exist.")

    def _pipe_select(self, *tree_files: str) -> None:
        """Select pipeline based on tree files."""
        cluster_num = len(tree_files)
        has_branch_length = _has_branch_length(*tree_files)
        if has_branch_length:
            self.func_dict['Gu2001'] = True
        if cluster_num >= 3:
            self.func_dict['FunDist'] = True
        if cluster_num == 2:
            self.func_dict['Effective'] = True
        if cluster_num == 3:
            self.func_dict['Asym'] = True
            self.func_dict['Type1Analysis'] = True

    def _pipeline(self) -> None:
        """Execute the calculation pipeline."""
        progress_bar = tqdm(total=sum(self.func_dict.values()), desc="Calculate pipeline")
        
        calculation_methods = {
            'Gu99': self._calculate_gu99,
            'Gu2001': self._calculate_gu2001,
            'Type1Analysis': self._calculate_type1analysis,
            'Type2': self._calculate_type2,
            'Effective': self._calculate_effective,
            'Fdr': self._calculate_fdr,
            'Asym': self._calculate_asym
        }

        for method, func in calculation_methods.items():
            if self.func_dict[method]:
                try:
                    func(progress_bar)
                except Exception as e:
                    printv(f"Error occurred during {method} calculation: {e}")
                    raise

    def _calculate_gu99(self, progress_bar: tqdm) -> None:
        gu99 = Gu99(self.aln_file, *self.tree_files)
        self.gu99_results = gu99.results
        self.gu99_summary = gu99.summary
        progress_bar.update(1)
        progress_bar.set_description("Gu99 calculation running...")
        progress_bar.refresh()
        
        if self.func_dict['FunDist']:
            self.fundist_results = gu99.fundist()
            progress_bar.update(1)
            progress_bar.set_description("Function distance calculation running...")
            progress_bar.refresh()

    def _calculate_gu2001(self, progress_bar: tqdm) -> None:
        gu2001 = Gu2001(self.aln_file, *self.tree_files)
        self.gu2001_results = gu2001.results
        self.gu2001_summary = gu2001.summary
        progress_bar.update(1)
        progress_bar.set_description("Gu2001 distance calculation running...")
        progress_bar.refresh()

    def _calculate_type1analysis(self, progress_bar: tqdm) -> None:
        toa = TypeOneAnalysis(self.aln_file, *self.tree_files)
        self.type1analysis_summary = toa.summary
        self.type1analysis_results = toa.results
        progress_bar.update(1)
        progress_bar.set_description("Type one Analysis distance calculation running...")
        progress_bar.refresh()

    def _calculate_type2(self, progress_bar: tqdm) -> None:
        type2 = Type2(self.aln_file, *self.tree_files)
        self.type2_summary = type2.summary
        self.type2_results = type2.results
        progress_bar.update(1)
        progress_bar.set_description("Type two calculation running...")
        progress_bar.refresh()

    def _calculate_effective(self, progress_bar: tqdm) -> None:
        effective = Effective(self.aln_file, *self.tree_files)
        self.type1_effective = effective.type1_results
        self.type2_effective = effective.type2_results
        progress_bar.update(1)
        progress_bar.set_description("Effective number of site calculation running...")
        progress_bar.refresh()

    def _calculate_fdr(self, progress_bar: tqdm) -> None:
        fdr = Fdr(self.aln_file, *self.tree_files)
        self.type1_fdr = fdr.type1_results
        self.type2_fdr = fdr.type2_results
        progress_bar.update(1)
        progress_bar.set_description("Fdr calculation running...")
        progress_bar.refresh()

    def _calculate_asym(self, progress_bar: tqdm) -> None:
        asym = Asym(self.aln_file, *self.tree_files)
        self.asym_results = asym.results
        progress_bar.update(1)
        progress_bar.set_description("Asym test running...")
        progress_bar.refresh()

    def _result_summary(self) -> Dict[str, Any]:
        """Generate result summary."""
        _result_summary = {}
        if self.func_dict['Gu99']:
            _result_summary['gu99_summary'] = self.gu99_summary
            _result_summary['gu99_results'] = self.gu99_results
        if self.func_dict['Gu2001']:
            _result_summary['gu2001_summary'] = self.gu2001_summary
            _result_summary['gu2001_results'] = self.gu2001_results
        if self.func_dict['Type1Analysis']:
            _result_summary['type1analysis_summary'] = self.type1analysis_summary
            _result_summary['type1analysis_results'] = self.type1analysis_results
        if self.func_dict['Type2']:
            _result_summary['type2_summary'] = self.type2_summary
            _result_summary['type2_results'] = self.type2_results
        if self.func_dict['Effective']:
            _result_summary['type1_effective'] = self.type1_effective
            _result_summary['type2_effective'] = self.type2_effective
        if self.func_dict['Fdr']:
            _result_summary['type1_fdr'] = self.type1_fdr
            _result_summary['type2_fdr'] = self.type2_fdr
        if self.func_dict['FunDist']:
            _result_summary['fundist_results'] = self.fundist_results
        if self.func_dict['Asym']:
            _result_summary['asym_results'] = self.asym_results
        return _result_summary

    def __str__(self) -> str:
        """
        String representation of CalPipe.

        Returns:
            str: A string describing the calculation pipeline and available results
        """
        strtext = "Diverge calculation pipeline\n"
        for i, j in enumerate(self.detail):
            strtext += f"step{i+1}: {j}\n"
        strtext += "#####################\n"
        strtext += f"You can get the result by calling the result_summary attribute or the specific attribute as follow:\n {self.result_summary.keys()}"
        return strtext

    def __repr__(self) -> str:
        """Representation of CalPipe."""
        return self.__str__()

def _has_branch_length(*tree_files: str) -> bool:
    """
    Check if all tree files have branch lengths.

    Args:
        *tree_files (str): Tree file paths

    Returns:
        bool: True if all trees have branch lengths, False otherwise

    Raises:
        Exception: If there's an inconsistency in branch length presence across trees
    """
    k = sum(_check_tree(tree_file) for tree_file in tree_files)
    if k not in (0, len(tree_files)):
        raise Exception('Tree file error, please check your tree files')
    return k == len(tree_files)

def _check_tree(tree_file: str) -> bool:
    """
    Check if a tree file has branch lengths.

    Args:
        tree_file (str): Path to the tree file

    Returns:
        bool: True if the tree has branch lengths, False otherwise
    """
    tree = Phylo.read(tree_file, 'newick')
    return any(clade.branch_length for clade in tree.find_clades())

def plot_phylogenetic_tree(tree: BaseTree.Tree, font_size: int = 5, figsize: Tuple[int, int] = (5, 6)) -> None:
    """
    Plot a phylogenetic tree.

    Args:
        tree (BaseTree.Tree): The phylogenetic tree to plot
        font_size (int, optional): Font size for the plot. Defaults to 5.
        figsize (Tuple[int, int], optional): Figure size. Defaults to (5, 6).
    """
    plt.rc('font', size=font_size)
    fig = plt.figure(figsize=figsize, dpi=400)
    axes = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, axes=axes)
    plt.show()

def read_tree(tree_file: str) -> BaseTree.Tree:
    """
    Read a tree file in newick or nexus format.

    Args:
        tree_file (str): Path to the tree file

    Returns:
        BaseTree.Tree: The parsed phylogenetic tree

    Raises:
        ValueError: If the tree file is not in newick or nexus format
    """
    try:
        return Phylo.read(tree_file, 'newick')
    except ValueError:
        try:
            return Phylo.read(tree_file, 'nexus')
        except ValueError:
            raise ValueError(f"Tree file '{tree_file}' is not in newick or nexus format.")

def printv(*text: Any, show_time: bool = True, verbose: bool = True) -> None:
    """
    Print text with time and verbose.

    Args:
        *text (Any): Text to print
        show_time (bool, optional): Whether to show time. Defaults to True.
        verbose (bool, optional): Whether to print. Defaults to True.
    """    
    if verbose:
        if show_time:
            print(time.ctime() + "\t", *text)
        else:
            print(*text)

@contextlib.contextmanager
def tqdm_joblib(*args: Any, **kwargs: Any) -> Generator[tqdm, None, None]:
    """
    Context manager to patch joblib to report into tqdm progress bar.

    Args:
        *args (Any): Arguments for tqdm
        **kwargs (Any): Keyword arguments for tqdm

    Yields:
        Generator[tqdm, None, None]: A tqdm progress bar object
    """
    tqdm_object = tqdm(*args, **kwargs)
    class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):
        def __init__(self, *args: Any, **kwargs: Any) -> None:
            super().__init__(*args, **kwargs)

        def __call__(self, *args: Any, **kwargs: Any) -> Any:
            tqdm_object.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)

    old_batch_callback = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
    try:
        yield tqdm_object
    finally:
        joblib.parallel.BatchCompletionCallBack = old_batch_callback
        tqdm_object.close()
        
        

def tree_construct(aln,method='nj',dist_calc='identity'):
    calculator = DistanceCalculator(dist_calc)
    constructor = DistanceTreeConstructor()
    printv("Calculating distance matrix from msa using",dist_calc,"method...")
    dm = calculator.get_distance(aln)
    printv("Constructing phylogenetic tree using",method,"method...")
    if method == 'nj':
        tree = constructor.nj(dm)
    else:
        tree = constructor.upgma(dm)
    for clade in tree.get_nonterminals():
        clade.name = None
    return tree

def aln_to_tree(aln_file,method='nj'):
    try:
        aln = AlignIO.read(aln_file,'clustal')
    except ValueError:
        aln = AlignIO.read(aln_file,'fasta')
    return tree_construct(aln,method)

# from a Bio.Phylo.BaseTree.Tree object, return a distance matrix
# distance is defined as the mini number of nodes between two leaves  
def tree_to_dist(tree):
    if tree.total_branch_length() == 0:
        for clade in tree.find_clades():
            clade.branch_length = 1
    clades_name = [clade.name for clade in tree.get_terminals() if clade.name != '' ]
    dist = pd.DataFrame(index=clades_name, columns=clades_name)
    for clade in tree.get_terminals():
        for other in tree.get_terminals():
            if clade.name != other.name:
                dist.loc[clade.name, other.name] = tree.distance(clade, other)
    return dist.fillna(0)

def auto_split(aln_file,plot=False,exclude_level=1,dist_calc='blosum62',tree_construct_method='nj',cluster_method='ward',n_clusters=None,root_at_midpoint=True):
    printv("Running subtree auto split process")
    try:
        aln = AlignIO.read(aln_file,'clustal')
    except ValueError:
        aln = AlignIO.read(aln_file,'fasta')
    
    
    tree = tree_construct(aln,method=tree_construct_method,dist_calc=dist_calc)
    if root_at_midpoint==True:
        tree.root_at_midpoint()
    if plot:
        printv("Phylogenetic tree of the alingnment file:")
        plot_phylogenetic_tree(tree)
    sub_trees,subtree_clades,exclude_list = split_tree2(tree,n_clusters=n_clusters,exclude_level=exclude_level,cluster_method=cluster_method,plot=plot)
    if plot:
        printv("Subtrees:")
        for subtree in sub_trees:
            plot_phylogenetic_tree(subtree)
    return sub_trees,subtree_clades,exclude_list

def split_tree2(tree,n_clusters=None,exclude_level=0,cluster_method='ward',plot=True):
    tree_ = copy.deepcopy(tree)
    exclude_list = []
    for clade in tree_.get_terminals():
        trace = tree_.trace(tree_.root,clade)
        if len(trace) <= exclude_level+1:
            tree_.prune(clade)
            exclude_list.extend([clade.name])
    printv("Prune terminals:",exclude_list)
    tree_dist = tree_to_dist(tree_)
    condensed_tree_dist = tree_dist.values[np.triu_indices(len(tree_dist), k=1)]
    Z = linkage(condensed_tree_dist, method=cluster_method)
    if n_clusters is None:
        last = Z[-10:, 2]
        last_rev = last[::-1]
        idxs = np.arange(1, len(last) + 1)
        if plot==True:
            plt.plot(idxs, last_rev)
        acceleration = np.diff(last, 2)  
        acceleration_rev = acceleration[::-1]
        n_clusters = acceleration_rev.argmax() + 2  
    printv("recommend n_cluster:",n_clusters)
    clusters = cut_tree(Z, n_clusters=n_clusters).squeeze()
    sub_trees = []
    subtree_clades = []
    for i in range(n_clusters):
        subtree_clades.extend([tree_dist.index[np.where(clusters == i)].values])
        root = tree.common_ancestor(tree_dist.index[np.where(clusters == i)].values)
        sub_tree = Phylo.BaseTree.Tree.from_clade(root)
        sub_tree.rooted = False
        sub_trees.extend([sub_tree])
    return sub_trees,subtree_clades,exclude_list

def save_subtrees(subtrees, directory, format='newick'):
    if not os.path.exists(directory):
        os.makedirs(directory)
    for i, subtree in enumerate(subtrees):
        filename = f"{directory}/subtree_{i}.{format}"
        Phylo.write(subtree, filename, format)
        printv(f"Subtree {i} saved as {filename}")

# split tree into clusters by distance_mtx
def split_tree(tree_file,n_clusters=2,height=None):
    tree = read_tree(tree_file)
    # construct the distance matrix 
    tree_dist = tree_to_dist(tree)
    ## transform distance matrix into a condensed matrix
    condensed_tree_dist = tree_dist.values[np.triu_indices(len(tree_dist), k=1)]
    Z = linkage(condensed_tree_dist, method='complete')
    clusters = cut_tree(Z,n_clusters=n_clusters).reshape(-1)
    ## from clusters to split tree object and return a list of split tree
    sub_trees = []
    root_distance = []
    for i in range(n_clusters):
        root = tree.common_ancestor(tree_dist.index[np.where(clusters == i)].values)
        # function = lambda x:x-1 if x>i else x
        # index = function(np.where(clusters == i)[0])
        root_distance.append(max([root.distance(i) for i in root.get_terminals()]))
        sub_trees.append(root)

    return sub_trees,root_distance

def plot_msa(aln,aln_format='clustal',color_scheme="None",marker_list=None,wrap_length=60,show_grid=True,show_plot=True,save_plot=False,path="./msa.png",marker='v',marker_color="red"):
    mv=MsaViz(aln,color_scheme=color_scheme,wrap_length=wrap_length,show_grid=show_grid, format =aln_format)
    if marker_list is not None:
        mv.add_markers(marker_list,marker=marker,color=marker_color)
    if show_plot:
        mv.plotfig()
    if save_plot:
        mv.savefig(path,pad_inches=0.5)
    return mv

def view_cutoff_msa(aln,results,colname,color_scheme="None",cutoff=0.5,show_grid=True,show_plot=True,save_plot=False,path="./msa.png",marker='v',marker_color="red",aln_format='clustal'):
    """view msa with cutoff value
    color_scheme : str | None, optional
            Color scheme. If None, `Zappo`(AA) or `Nucleotide`(NT) is set.
            [`Clustal`|`Zappo`|`Taylor`|`Flower`|`Blossom`|`Sunset`|`Ocean`|
            `Hydrophobicity`|`HelixPropensity`|`StrandPropensity`|`TurnPropensity`|
            `BuriedIndex`|`Nucleotide`|`Purine/Pyrimidine`|`None`]
    marker : str | Marker type of matplotlib.
            See https://matplotlib.org/stable/api/markers_api.html for details.
    """
    try:
        results.loc[:,colname]
    except(KeyError):
        raise Exception(f'column {colname} not found in results')
    marker_list = results[results.loc[:,colname] > cutoff].index.tolist()
    mv = plot_msa(aln,aln_format=aln_format,marker_list=marker_list,show_grid=show_grid,show_plot=show_plot,save_plot=save_plot,path=path,marker=marker,color_scheme=color_scheme,marker_color=marker_color)
    return mv

def view_cutoff_msa2(aln_path,results,colname,color_scheme="None",cutoff=0.5,show_grid=True,show_plot=True,save_plot=False,path="./msa.png",marker='v',marker_color="red",aln_format='clustal'):
    """view msa with cutoff value
    color_scheme : str | None, optional
            Color scheme. If None, `Zappo`(AA) or `Nucleotide`(NT) is set.
            [`Clustal`|`Zappo`|`Taylor`|`Flower`|`Blossom`|`Sunset`|`Ocean`|
            `Hydrophobicity`|`HelixPropensity`|`StrandPropensity`|`TurnPropensity`|
            `BuriedIndex`|`Nucleotide`|`Purine/Pyrimidine`|`None`]
    marker : str | Marker type of matplotlib.
            See https://matplotlib.org/stable/api/markers_api.html for details.
    """
    try :
        aln = AlignIO.read(aln_path,'clustal')
    except ValueError:
        try:
            aln = AlignIO.read(aln_path,'fasta')
        except ValueError:
            raise ValueError("The alignment file is not in fasta or clustal format")
        
    try:
        results.loc[:,colname]
    except(KeyError):
        raise Exception(f'column {colname} not found in results')
    marker_list = results[results.loc[:,colname] > cutoff].index.tolist()
    mv = plot_msa(aln,aln_format=aln_format,marker_list=marker_list,show_grid=show_grid,show_plot=show_plot,save_plot=save_plot,path=path,marker=marker,color_scheme=color_scheme,marker_color=marker_color)
    return mv

def pre_check_tree(*tree_files):
    for tree_file in tree_files:
        tree = Phylo.read(tree_file, 'newick')
        # for clade in tree.find_clades():
        #     clade.branch_length = 1
        # get tree depth
        # for clade in tree.get_terminals():
        #     trace = tree.trace(tree.root,clade)
        tree_depth = max([len(tree.trace(tree.root,clade)) for clade in tree.get_terminals()])
        if tree_depth < 3:
            raise Exception(f'tree depth is less than 3, please check your tree file{tree_file}')
            return False    
    return True

def get_genefam_member(tree):
    tree_terminal_names = [clade.name for clade in tree.get_terminals()]
    # calde_name : species_name + "_" + gene_name
    # get gene list
    gene_list = [i.split('_')[-1] for i in tree_terminal_names]
    gene_list = list(set(gene_list))
    return gene_list

def get_genefam_clusters(tree_obj,gene_fam,exclude_list=[],exclude_level=1):
    tree_list = []
    for clade in tree_obj.get_terminals():
        trace = tree_obj.trace(tree_obj.root,clade)
        if len(trace) <= exclude_level+1:
            exclude_list.extend([clade.name])
    exclude_list = list(set(exclude_list))
    if isinstance(gene_fam,list):
        for gene in gene_fam:
            gene_terminals = [clade for clade in tree_obj.get_terminals() if gene in clade.name and clade.name not in exclude_list]
            subtree_clade = tree_obj.common_ancestor(gene_terminals)
            subtree = Phylo.BaseTree.Tree.from_clade(subtree_clade)
            tree_list.extend([subtree])
        return tree_list,gene_fam
    elif isinstance(gene_fam,dict):
        for cluster_name,gene_list in gene_fam.items():
            gene_terminals = [clade for clade in tree_obj.get_terminals() if any(gene in clade.name for gene in gene_list) and clade.name not in exclude_list]
            subtree_clade = tree_obj.common_ancestor(gene_terminals)
            subtree = Phylo.BaseTree.Tree.from_clade(subtree_clade)
            tree_list.extend([subtree])
        return tree_list,gene_fam.keys()