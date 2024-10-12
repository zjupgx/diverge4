import os
import time
from itertools import combinations
from typing import Any, List, Optional, Tuple

import numpy as np
import pandas as pd
from Bio import Phylo, AlignIO
from Bio.Phylo.BaseTree import Tree, Clade
from joblib import Parallel, delayed
from scipy.cluster.hierarchy import linkage
from skbio import TabularMSA, Protein, DistanceMatrix
from skbio.sequence.distance import hamming
from skbio.tree import nj
from tqdm import tqdm

from diverge import Gu99, Type2
from diverge.utils import printv, tqdm_joblib

def timeit(func):
    """
    Decorator to measure the execution time of a function.

    Args:
        func (callable): The function to be timed.

    Returns:
        callable: The wrapped function.
    """
    def timed(*args: Any, **kwargs: Any) -> Any:
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        print(f"{func.__name__} took {end_time - start_time:.4f} seconds")
        return result
    return timed

def pre_check_tree(trees: List[Tree]) -> bool:
    """
    Check if all trees have a depth of at least 3.

    Args:
        trees (List[Tree]): List of trees to check.

    Returns:
        bool: True if all trees have depth >= 3, False otherwise.
    """
    for tree in trees:
        tree_depth = max([len(tree.trace(tree.root, clade)) for clade in tree.get_terminals()])
        if tree_depth < 3:
            return False    
    return True

def aln_read(aln_file: str) -> AlignIO.MultipleSeqAlignment:
    """
    Read an alignment file in clustal or fasta format.

    Args:
        aln_file (str): Path to the alignment file.

    Returns:
        AlignIO.MultipleSeqAlignment: The parsed alignment.

    Raises:
        Exception: If the alignment file is not in fasta or clustal format.
    """
    try:
        aln = AlignIO.read(aln_file, 'clustal')
    except ValueError:
        try:
            aln = AlignIO.read(aln_file, 'fasta')
        except ValueError:
            print("==============\n", aln_file, "==============\n")
            raise Exception("The alignment file is not in fasta or clustal format")
    return aln

def sk_aln_read(aln_file: str) -> TabularMSA:
    """
    Read an alignment file using scikit-bio.

    Args:
        aln_file (str): Path to the alignment file.

    Returns:
        TabularMSA: The parsed alignment.
    """
    aln = TabularMSA.read(aln_file, constructor=Protein)
    try:
        aln.reassign_index(minter='id')
    except KeyError:
        pass
    return aln

def skbio_to_biopython_tree(skbio_tree: Any) -> Tree:
    """
    Convert a scikit-bio tree to a Biopython tree.

    Args:
        skbio_tree (Any): A scikit-bio tree object.

    Returns:
        Tree: The converted Biopython tree.
    """
    def convert_node(skbio_node: Any) -> Clade:
        if skbio_node.is_tip():
            return Clade(branch_length=skbio_node.length, name=skbio_node.name)
        else:
            clade = Clade(branch_length=skbio_node.length)
            clade.clades = [convert_node(child) for child in skbio_node.children]
            return clade
    root_clade = convert_node(skbio_tree.root())
    biopython_tree = Tree(root=root_clade, rooted=True)
    return biopython_tree

@timeit
def tree_construct(aln: TabularMSA) -> Tuple[DistanceMatrix, Tree]:
    """
    Construct a tree from an alignment.

    Args:
        aln (TabularMSA): The alignment.

    Returns:
        Tuple[DistanceMatrix, Tree]: The distance matrix and constructed tree.
    """
    dm = DistanceMatrix.from_iterable(aln, metric=hamming, keys=aln.index)
    tree = nj(dm)
    biopython_tree = skbio_to_biopython_tree(tree)
    return dm, biopython_tree

def re_clean_tree(tree: Tree) -> Tree:
    """
    Remove names from non-terminal nodes of a tree.

    Args:
        tree (Tree): The tree to clean.

    Returns:
        Tree: The cleaned tree.
    """
    for clade in tree.get_nonterminals():
        clade.name = None
    return tree

def get_cluster(aln: TabularMSA, *tree_files: str, trees: List[Tree] = []) -> Tuple[List[int], int, dict]:
    """
    Get cluster information from alignment and trees.

    Args:
        aln (TabularMSA): The alignment.
        *tree_files (str): Paths to tree files.
        trees (List[Tree], optional): List of pre-loaded trees. Defaults to [].

    Returns:
        Tuple[List[int], int, dict]: Cluster assignments, number of clusters, and tree dictionary.
    """
    names = aln.index.to_list()
    tree_cluster = [0] * len(names)
    i = 1
    tree_dict = {}
    if trees:
        for tree in trees:
            tree_id = f"subtree_{i}"
            tree_dict[i] = tree_id
            tree_terminal = [i.name for i in tree.get_terminals()]
            t_list = [names.index(j) for j in tree_terminal]
            for k in t_list:
                tree_cluster[k] = i
            i += 1
    else:
        def read_tree(tree_file: str) -> Tuple[str, List[int]]:
            tree = Phylo.read(tree_file, 'newick')
            tree_id = os.path.basename(tree_file).split('.')[0]
            tree_terminal = [i.name for i in tree.get_terminals()]
            t_list = [names.index(j) for j in tree_terminal]
            return tree_id, t_list
        for tree_file in tree_files:
            tree_id, t_list = read_tree(tree_file)
            tree_dict[i] = tree_id
            for k in t_list:
                tree_cluster[k] = i
            i += 1
    return tree_cluster, i - 1, tree_dict

def sub_dm(dm: DistanceMatrix, c_list: np.ndarray) -> DistanceMatrix:
    """
    Create a sub-distance matrix based on a boolean mask.

    Args:
        dm (DistanceMatrix): The original distance matrix.
        c_list (np.ndarray): Boolean mask for selecting entries.

    Returns:
        DistanceMatrix: The sub-distance matrix.
    """
    sub_ids = np.array(dm.ids)[c_list]
    sub_data = dm.filter(sub_ids).data
    return DistanceMatrix(sub_data, ids=sub_ids)

def sep_cluster(tree_cluster: List[int], cluster_num: int) -> Tuple[List[np.ndarray], List[Tuple[List[int], List[int]]]]:
    """
    Separate clusters into groups.

    Args:
        tree_cluster (List[int]): Cluster assignments.
        cluster_num (int): Number of clusters.

    Returns:
        Tuple[List[np.ndarray], List[Tuple[List[int], List[int]]]]: Cluster list and group list.
    """
    group_list = get_group_list(cluster_num)
    cluster_list = []
    for group in group_list:
        cluster = np.zeros(len(tree_cluster))
        group1, group2 = group[0], group[1]
        cluster[np.isin(tree_cluster, group1)] = 1
        cluster[np.isin(tree_cluster, group2)] = 2
        cluster_list.extend([cluster])
    return cluster_list, group_list

def get_group_list(group_num: int) -> List[Tuple[List[int], List[int]]]:
    """
    Generate all possible group combinations.

    Args:
        group_num (int): Number of groups.

    Returns:
        List[Tuple[List[int], List[int]]]: List of group combinations.
    """
    nums = frozenset(range(1, group_num + 1))
    group_list = []
    for num in range(1, group_num // 2 + 1):
        for combo in combinations(nums, num):
            group1 = frozenset(combo)
            group2 = nums - group1
            if len(group1) <= len(group2):
                group_list.append((list(group1), list(group2)))
    return group_list

@timeit
def tree_reconstruct(dm: DistanceMatrix, cluster: np.ndarray) -> Tuple[Tree, Tree]:
    """
    Reconstruct trees based on a distance matrix and cluster assignments.

    Args:
        dm (DistanceMatrix): The distance matrix.
        cluster (np.ndarray): Cluster assignments.

    Returns:
        Tuple[Tree, Tree]: Two reconstructed trees.
    """
    def preserve_original_names(tree: Any, original_names: List[str]) -> Any:
        name_map = {name.replace('_', ' '): name for name in original_names}
        for tip in tree.tips():
            if tip.name in name_map:
                tip.name = name_map[tip.name]
        return tree
    
    cluster1_list = np.isin(cluster, 1)
    cluster2_list = np.isin(cluster, 2)
    sub_dm1 = sub_dm(dm, cluster1_list)
    sub_dm2 = sub_dm(dm, cluster2_list)
    
    original_names1 = sub_dm1.ids
    original_names2 = sub_dm2.ids
    
    tree1 = preserve_original_names(nj(sub_dm1), original_names1)
    tree2 = preserve_original_names(nj(sub_dm2), original_names2)
    
    return skbio_to_biopython_tree(tree1), skbio_to_biopython_tree(tree2)

@timeit
def process_tree(aln_file: str, super_cluster: List[Tree], sp_type: int) -> Tuple[Optional[List[float]], Optional[List[str]], Optional[dict]]:
    """
    Process trees for functional divergence analysis.

    Args:
        aln_file (str): Path to the alignment file.
        super_cluster (List[Tree]): List of trees to process.
        sp_type (int): Type of functional divergence analysis (1 or 2).

    Returns:
        Tuple[Optional[List[float]], Optional[List[str]], Optional[dict]]: 
        Results, positions, and summary of the analysis.
    """
    if pre_check_tree(super_cluster):
        if sp_type == 1:
            calculator = Gu99(aln_file, trees=super_cluster)
        else:
            calculator = Type2(aln_file, trees=super_cluster)
        summary = calculator.summary()
        position = calculator.results().index.values.tolist()
        results = calculator.results().values.tolist()
        return results, position, summary
    else:
        print(f'Tree depth is less than 3, please check your tree file')
        return None, None, None

@timeit
def get_super_cluster_pp(aln_file: str, *tree_files: str, sp_type: int = 1, trees: List[Tree] = [], verbose: bool = True) -> Tuple[np.ndarray, List[List[Tree]], List[Tuple[List[int], List[int]]], dict, List[str], List[dict]]:
    """
    Perform super cluster analysis for functional divergence.

    Args:
        aln_file (str): Path to the alignment file.
        *tree_files (str): Paths to tree files.
        sp_type (int, optional): Type of functional divergence analysis. Defaults to 1.
        trees (List[Tree], optional): Pre-loaded trees. Defaults to [].
        verbose (bool, optional): Whether to print progress. Defaults to True.

    Returns:
        Tuple[np.ndarray, List[List[Tree]], List[Tuple[List[int], List[int]]], dict, List[str], List[dict]]:
        Results array, super cluster list, group list, tree dictionary, position list, and summary list.
    """
    aln = sk_aln_read(aln_file)
    tree_cluster, cluster_num, tree_dict = get_cluster(aln, *tree_files, trees=trees)
    cluster_list, group_list = sep_cluster(tree_cluster, cluster_num)
    dm, _ = tree_construct(aln)

    def process_cluster_and_tree(cluster: np.ndarray) -> Tuple[List[Tree], Tuple[Optional[List[float]], Optional[List[str]], Optional[dict]]]:
        tree1, tree2 = tree_reconstruct(dm, cluster)
        tree1, tree2 = re_clean_tree(tree1), re_clean_tree(tree2)
        super_cluster = [tree1, tree2]
        return super_cluster, process_tree(aln_file, super_cluster, sp_type)

    n_jobs = 8  
    if verbose:
        printv(f"Running type{sp_type} functional divergence super cluster analysis...")
        with tqdm_joblib(desc="Processing super cluster groups", total=len(cluster_list)) as progress_bar:
            results = Parallel(n_jobs=n_jobs)(delayed(process_cluster_and_tree)(cluster) for cluster in cluster_list)
        printv('Finish super cluster analysis.')
    else:
        results = Parallel(n_jobs=n_jobs)(delayed(process_cluster_and_tree)(cluster) for cluster in cluster_list)

    super_cluster_list = []
    results_list = []
    position_list = []
    summary_list = []
    new_group_list = []
    for i, (super_cluster, (results, position, summary)) in enumerate(results):
        if results is not None:
            super_cluster_list.append(super_cluster)
            new_group_list.append(group_list[i])
            results_list.append(results)
            summary_list.append(summary)
            if not position_list:
                position_list = position
    results_array = np.reshape(np.array(results_list), (len(super_cluster_list), -1))
    return results_array, super_cluster_list, new_group_list, tree_dict, position_list, summary_list

class SuperCluster:
    """Class for performing super cluster analysis of functional divergence."""

    def __init__(self, aln_file: str, *tree_files: str, trees: List[Tree] = [], sp_type: int = 1, verbose: bool = True):
        """
        Initialize SuperCluster.

        Args:
            aln_file (str): Path to the alignment file.
            *tree_files (str): Paths to tree files.
            trees (List[Tree], optional): Pre-loaded trees. Defaults to [].
            sp_type (int, optional): Type of functional divergence analysis. Defaults to 1.
            verbose (bool, optional): Whether to print progress. Defaults to True.
        """
        self.pp_list, self.tree_list, self.group_list, self.tree_dict, self.position_list, self.summary_list = get_super_cluster_pp(aln_file, *tree_files, sp_type=sp_type, trees=trees, verbose=verbose)
        self.results = pd.DataFrame(self.pp_list, index=[f"{i}" for i in self.group_list], columns=self.position_list)
        self.summary = self.get_summary()
    def get_summary(self):
        col_name = [f"{i}" for i in self.group_list]
        summary_df = pd.DataFrame(index=list(self.summary_list[0].index.values)+["Qk>0.5","Qk>0.67","Qk>0.9"],columns=col_name)
        for i in range(len(self.group_list)):
          # summary_df.iloc[0,i] = f"{round(self.param_list[i][0][0],2)}±{round(self.param_list[i][1][0],2)}"
          # summary_df.iloc[1,i] = round(self.param_list[i][2][0],2)
          index_s = self.summary_list[i].shape[0]
          summary_df.iloc[:index_s,i] = self.summary_list[i].iloc[:,0]
          summary_df.iloc[index_s,i] = sum(self.pp_list[i,:]>0.5)
          summary_df.iloc[index_s+1,i] = sum(self.pp_list[i,:]>0.67)
          summary_df.iloc[index_s+2,i] = sum(self.pp_list[i,:]>0.9)
        return summary_df
    
# Main execution
if __name__ == "__main__":
    aln_file = "E:/verysync/diverge_pybind/web/statics/example_data/ERBB_family.fas"
    tree_files = [
        "E:/verysync/diverge_pybind/web/statics/example_data/EGFR.tree",
        "E:/verysync/diverge_pybind/web/statics/example_data/ERBB2.tree",
        "E:/verysync/diverge_pybind/web/statics/example_data/ERBB3.tree",
        "E:/verysync/diverge_pybind/web/statics/example_data/ERBB4.tree"
    ]
    super_cluster = SuperCluster(aln_file, *tree_files)
    print(super_cluster.results)