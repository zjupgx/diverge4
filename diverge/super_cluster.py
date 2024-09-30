import itertools
import os
import time
from itertools import combinations
from typing import Any

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

# Utility functions
def timeit(func):
    """Decorator to measure the execution time of a function."""
    def timed(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        print(f"{func.__name__} took {end_time - start_time:.4f} seconds")
        return result
    return timed

def pre_check_tree(trees):
    for tree in trees:
        tree_depth = max([len(tree.trace(tree.root, clade)) for clade in tree.get_terminals()])
        if tree_depth < 3:
            return False    
    return True

# File reading functions

def aln_read(aln_file):
    try:
        aln = AlignIO.read(aln_file, 'clustal')
    except ValueError:
        try:
            aln = AlignIO.read(aln_file, 'fasta')
        except ValueError:
            print("==============\n", aln_file, "==============\n")
            raise Exception("The alignment file is not in fasta or clustal format")
    return aln

def sk_aln_read(aln_file):
    aln = TabularMSA.read(aln_file, constructor=Protein)
    aln.reassign_index(minter='id')
    return aln

# Tree manipulation functions
def skbio_to_biopython_tree(skbio_tree):
    def convert_node(skbio_node):
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
def tree_construct(aln):
    dm = DistanceMatrix.from_iterable(aln, metric=hamming, keys=aln.index)
    tree = nj(dm)
    biopython_tree = skbio_to_biopython_tree(tree)
    return dm, biopython_tree

def re_clean_tree(tree):
    for clade in tree.get_nonterminals():
        clade.name = None
    return tree

# Cluster and group functions

def get_cluster(aln, *tree_files, trees: list = []):
    names = [record.metadata['id'] for record in aln]
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
        def read_tree(tree_file):
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

def sub_dm(dm, c_list):
    sub_ids = np.array(dm.ids)[c_list]
    sub_data = dm.filter(sub_ids).data
    return DistanceMatrix(sub_data, ids=sub_ids)

def sep_cluster(tree_cluster, cluster_num):
    group_list = get_group_list(cluster_num)
    cluster_list = []
    for group in group_list:
        cluster = np.zeros(len(tree_cluster))
        group1, group2 = group[0], group[1]
        cluster[np.isin(tree_cluster, group1)] = 1
        cluster[np.isin(tree_cluster, group2)] = 2
        cluster_list.extend([cluster])
    return cluster_list, group_list

def get_group_list(group_num):
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
def tree_reconstruct(dm, cluster):
    def preserve_original_names(tree, original_names):
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
def process_tree(aln_file, super_cluster, sp_type):
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

# Main function for super cluster analysis
@timeit
def get_super_cluster_pp(aln_file, *tree_files, sp_type: int = 1, trees: list = [], verbose=True):
    aln = sk_aln_read(aln_file)
    tree_cluster, cluster_num, tree_dict = get_cluster(aln, *tree_files, trees=trees)
    cluster_list, group_list = sep_cluster(tree_cluster, cluster_num)
    dm, _ = tree_construct(aln)

    def process_cluster_and_tree(cluster):
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

    for super_cluster, (results, position, summary) in results:
        if results is not None:
            super_cluster_list.append(super_cluster)
            results_list.append(results)
            summary_list.append(summary)
            if not position_list:
                position_list = position
    results_array = np.reshape(np.array(results_list), (len(super_cluster_list), -1))
    return results_array, super_cluster_list, group_list, tree_dict, position_list, summary_list

# Main class
class SuperCluster():
    def __init__(self, aln_file, *tree_files, trees: list = [], sp_type: int = 1, verbose=True):
        """
        Attributes:
          sp_type (int): Using SuperCluster to calculate sp_type I or sp_type II functional divergence.
          aln_file (str): The path to the alignment file.
          tree_files (str): A list of paths to the phylogenetic tree files.
          result (pandas.DataFrame): A dataframe that stores the posterior probability of every site in the alignment for each group in the supercluster.
        """
        self.pp_list, self.tree_list, self.group_list, self.tree_dict, self.position_list, self.summary_list = get_super_cluster_pp(aln_file, *tree_files, sp_type=sp_type, trees=trees, verbose=verbose)
        self.result = pd.DataFrame(self.pp_list, index=[f"{i}" for i in self.group_list], columns=self.position_list)

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
    print(super_cluster.result)