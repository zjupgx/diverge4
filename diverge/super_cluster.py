import os
import time
from itertools import combinations
from typing import Any, List, Optional, Tuple, Dict
from functools import lru_cache

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
from Bio.Phylo.TreeConstruction import DistanceCalculator
from diverge import Gu99, Type2
from diverge.binding import Gu99Batch
from diverge.utils import printv, tqdm_joblib
import math
# BLOSUM62 background distribution for Jensen-Shannon divergence
BLOSUM62_BACKGROUND = [0.074, 0.052, 0.045, 0.054, 0.025, 0.034, 0.054, 0.074,
                       0.026, 0.068, 0.099, 0.058, 0.025, 0.047, 0.039, 0.057,
                       0.051, 0.013, 0.032, 0.073]

AMINO_ACIDS_NO_GAP = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
                      'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']


def gap_percentage(sequence: str) -> float:
    """Calculate the percentage of gaps in a sequence.
    
    Computes the proportion of gap characters ('-') in a protein or DNA sequence.
    Used for quality assessment and filtering of alignment columns.
    
    Args:
        sequence: Input sequence string that may contain gap characters.
        
    Returns:
        Percentage of gaps as a float between 0.0 and 1.0.
        Returns 0.0 for empty sequences.
        
    Example:
        >>> gap_percentage("ATCG--A-")
        0.375  # 3 gaps out of 8 positions
    """
    if not sequence:
        return 0.0
    num_gaps = sequence.count('-')
    return num_gaps / len(sequence)

def frequency_count(column: str, symbols: List[str], pseudocount: float = 0.000001) -> List[float]:
    """Calculate symbol frequency counts for a column.
    
    Computes normalized frequency counts for each symbol in the provided list,
    with optional pseudocount addition to avoid zero probabilities.
    
    Args:
        column: String representing an alignment column.
        symbols: List of symbols to count (e.g., amino acids).
        pseudocount: Small value added to each count to avoid zeros (default: 1e-6).
                   Set to 0 to disable pseudocount addition.
                   
    Returns:
        List of normalized frequencies (probabilities) for each symbol.
        Returns list of zeros if total count is zero.
        
    Example:
        >>> frequency_count("AAC", ["A", "C", "G"], 0.0)
        [0.667, 0.333, 0.0]  # 2 A's, 1 C, 0 G's out of 3 total
    """
    if pseudocount != 0:
        freq_counts = [pseudocount] * len(symbols)
    else:
        freq_counts = [0.0] * len(symbols)
        
    for i, symbol in enumerate(symbols):
        freq_counts[i] += column.count(symbol)
        
    total_counts = sum(freq_counts)
    if total_counts == 0:
        return [0.0] * len(symbols)
        
    return [count / total_counts for count in freq_counts]

def js_divergence(col1: str, col2: Optional[List[float]] = None) -> float:
    """Calculate Jensen-Shannon divergence for conservation scoring.
    
    Computes the Jensen-Shannon divergence between the amino acid composition
    of an alignment column and either another distribution or the BLOSUM62
    background distribution. The result is weighted by (1 - gap_percentage)
    to account for alignment quality.
    
    The Jensen-Shannon divergence is a symmetric measure of similarity between
    two probability distributions, ranging from 0 (identical) to 1 (completely different).
    
    Args:
        col1: Alignment column as string.
        col2: Optional probability distribution to compare against.
             If None, uses BLOSUM62 background distribution.
             
    Returns:
        Conservation score as float, weighted by alignment quality.
        Returns -1.0 if distributions have different lengths.
        Higher values indicate lower conservation (more divergent from background).
        
    Note:
        Uses base-2 logarithm for divergence calculation. Gap characters
        reduce the final score proportionally to their frequency.
    """
    fc1 = frequency_count(col1, AMINO_ACIDS_NO_GAP, 0.001)
    
    if col2 is not None:
        fc2 = col2
    else:
        fc2 = BLOSUM62_BACKGROUND
        
    if len(fc1) != len(fc2):
        return -1.0
        
    # Create R distribution
    r = [0.5 * fc1[i] + 0.5 * fc2[i] for i in range(len(fc1))]
    
    divergence = 0.0
    for i in range(len(fc1)):
        if r[i] != 0.0:
            if fc1[i] == 0.0:
                divergence += fc2[i] * math.log(fc2[i] / r[i], 2)
            elif fc2[i] == 0.0:
                divergence += fc1[i] * math.log(fc1[i] / r[i], 2)
            else:
                divergence += (fc1[i] * math.log(fc1[i] / r[i], 2) + 
                             fc2[i] * math.log(fc2[i] / r[i], 2))
                             
    divergence /= 2
    
    return (1 - gap_percentage(col1)) * divergence

def normalize_to_01(score_list: List[Optional[float]]) -> List[Optional[float]]:
    """Map values to range [0,1], ignoring None values.
    
    Performs min-max normalization on a list of scores, preserving None values.
    Useful for standardizing different types of scores to a common scale.
    
    Args:
        score_list: List of numerical scores, may contain None values.
        
    Returns:
        List of normalized scores in [0,1] range. None values are preserved.
        If all values are None or identical, handles edge cases appropriately.
        
    Note:
        If all valid scores are identical, they are mapped to 1.0.
        None values indicate positions that should be excluded from analysis.
    """
    new_scores = []
    valid_scores = [s for s in score_list if s is not None]
    
    if not valid_scores:
        return score_list
        
    min_val = min(valid_scores)
    max_val = max(valid_scores)
    score_range = max_val - min_val
    
    for score in score_list:
        if score is None:
            new_scores.append(None)
        elif score_range == 0:
            new_scores.append(1.0)
        else:
            new_scores.append((score - min_val) / score_range)
            
    return new_scores

def conservation_window_score(sdp_scores: List[Optional[float]], 
                            alignment: List[str], window_len: int, 
                            lambda_param: float = 0.7) -> List[Optional[float]]:
    """Apply conservation window heuristic to scores.
    
    Implements a sliding window approach to incorporate local sequence conservation
    into functional divergence scores. This method combines the original SDP
    (Site-specific Divergence Pattern) scores with conservation information
    from neighboring alignment positions.
    
    The algorithm:
    1. Calculates conservation scores for each alignment position using JS divergence
    2. Computes windowed conservation averages around each position
    3. Normalizes windowed conservation scores to [0,1]
    4. Combines original and conservation scores using linear interpolation
    
    Args:
        sdp_scores: Original site-specific divergence pattern scores.
        alignment: List of aligned sequences as strings.
        window_len: Half-window size for conservation calculation.
                   Total window = 2*window_len + 1.
        lambda_param: Weight for combining scores (default: 0.7).
                     Final score = lambda*original + (1-lambda)*conservation.
                     
    Returns:
        List of weighted scores combining original SDP and conservation information.
        Positions with None in original scores remain None.
        
    Note:
        Higher lambda_param values give more weight to original SDP scores.
        Window boundary effects are handled by using available positions only.
    """
    win_scores = sdp_scores[:]
    cons_scores = []
    
    # Calculate conservation scores for all columns
    for i in range(len(alignment[0])):
        column = ""
        for seq in alignment:
            if i < len(seq):
                column += seq[i]
        
        if sdp_scores[i] is None:
            cons_scores.append(None)
        else:
            cons_scores.append(js_divergence(column))
            
    # Calculate windowed conservation scores
    wincon_scores = []
    for i, sdp_score in enumerate(sdp_scores):
        if sdp_score is None:
            wincon_scores.append(None)
            continue
            
        window_sum = 0.0
        num_terms = 0.0
        
        start_pos = max(i - window_len, 0)
        end_pos = min(i + window_len + 1, len(sdp_scores))
        
        for j in range(start_pos, end_pos):
            if i != j and cons_scores[j] is not None and cons_scores[j] >= 0:
                num_terms += 1
                window_sum += cons_scores[j]
                
        if num_terms != 0:
            wincon_scores.append(window_sum / num_terms)
        else:
            wincon_scores.append(0.0)
            
    wincon_scores = normalize_to_01(wincon_scores)
    
    # Combine SDP and conservation scores
    for i, sdp_score in enumerate(sdp_scores):
        if sdp_score is not None and wincon_scores[i] is not None:
            win_scores[i] = ((1 - lambda_param) * wincon_scores[i] + lambda_param * sdp_scores[i])
                           
    return win_scores



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
        # Timing info removed for production
        return result
    return timed

def pre_check_tree(trees: List[Tree]) -> bool:
    """Check if all trees have a depth of at least 3.
    
    Validates that each phylogenetic tree in the list has sufficient depth
    for meaningful functional divergence analysis. Trees with insufficient
    depth may not provide reliable evolutionary signal.
    
    Args:
        trees: List of Biopython Tree objects to validate.
        
    Returns:
        True if all trees have depth >= 3, False if any tree fails the check.
        
    Note:
        Tree depth is calculated as the maximum path length from root
        to any terminal (leaf) node in the tree.
    """
    for tree in trees:
        tree_depth = max([len(tree.trace(tree.root, clade)) for clade in tree.get_terminals()])
        if tree_depth < 3:
            return False    
    return True

def pre_check(aln_file: str,trees: List[Tree]) -> bool:
    """Check if all trees have a depth of at least 3 and match alignment.
    
    Performs comprehensive validation of trees and alignment compatibility:
    1. Verifies that all tree terminal names exist in the alignment
    2. Ensures all trees have sufficient depth (>= 3) for analysis
    
    Args:
        aln_file: Path to the alignment file.
        trees: List of Biopython Tree objects to validate.
        
    Returns:
        True if all validation checks pass.
        
    Raises:
        ValueError: If any tree terminal is not found in alignment,
                   or if any tree has insufficient depth (<3).
                   
    Note:
        This function reads the alignment file to extract sequence IDs
        for comparison with tree terminal names.
    """
    aln = aln_read(aln_file)
    seq_ids = [record.id for record in aln]
    for tree in trees:
        for t in tree.get_terminals():
            if t.name not in seq_ids:
                raise ValueError(f"Tree {tree.name} has a terminal {t.name} that is not in the alignment")
        tree_depth = max([len(tree.trace(tree.root, clade)) for clade in tree.get_terminals()])
        if tree_depth < 3:
            raise ValueError('Tree depth is less than 3, please check your tree file')    
    return True

def aln_read(aln_file: str) -> AlignIO.MultipleSeqAlignment:
    """Read an alignment file in clustal or fasta format.
    
    Attempts to read an alignment file, trying clustal format first,
    then falling back to fasta format if that fails.
    
    Args:
        aln_file: Path to the alignment file.
        
    Returns:
        Parsed multiple sequence alignment object.
        
    Raises:
        Exception: If the file cannot be parsed as either clustal or fasta format.
        
    Note:
        The function automatically detects the format by attempting
        to parse with each format sequentially.
    """
    try:
        aln = AlignIO.read(aln_file, 'clustal')
    except ValueError:
        try:
            aln = AlignIO.read(aln_file, 'fasta')
        except ValueError:
            raise Exception("The alignment file is not in fasta or clustal format")
    return aln

def sk_aln_read(aln_file: str) -> TabularMSA:
    """Read an alignment file using scikit-bio.
    
    Reads a FASTA alignment file and converts it to a scikit-bio TabularMSA
    object. Attempts to reassign sequence indices using sequence IDs.
    
    Args:
        aln_file: Path to the alignment file (must be FASTA format).
        
    Returns:
        TabularMSA object containing protein sequences.
        
    Note:
        If sequence ID reassignment fails (due to missing or duplicate IDs),
        the function continues with default numeric indices.
    """
    from skbio import io
    # Read the alignment and convert the generator to TabularMSA
    sequences = list(io.read(aln_file, format='fasta', constructor=Protein))
    aln = TabularMSA(sequences)
    try:
        aln.reassign_index(minter='id')
    except KeyError:
        pass
    return aln


def skbio_to_biopython_tree(dm,skbio_tree: Any) -> Tree:
    """Convert a scikit-bio tree to a Biopython tree.
    
    Recursively converts a scikit-bio tree structure to the equivalent
    Biopython Tree object, preserving branch lengths and terminal names.
    Handles name formatting issues (spaces converted to underscores).
    
    Args:
        dm: Distance matrix containing valid sequence IDs for name validation.
        skbio_tree: A scikit-bio tree object to convert.
        
    Returns:
        Converted Biopython Tree object with equivalent structure.
        
    Raises:
        ValueError: If a terminal name cannot be found in the distance matrix IDs,
                   even after attempting space-to-underscore conversion.
                   
    Note:
        Internal nodes in the converted tree will have no names (None).
        Only terminal (leaf) nodes retain their original names.
    """
    def convert_node(skbio_node: Any) -> Clade:
        if skbio_node.is_tip():
            save_name = skbio_node.name
            if " " in save_name:
                # check if save_name in dm.ids
                if save_name not in dm.ids:
                    if save_name.replace(" ","_") in dm.ids:
                        save_name = save_name.replace(" ","_")
                    else:
                        raise ValueError(f"The sequence name {save_name} is not in the alignment")
            return Clade(branch_length=skbio_node.length, name=save_name)
        else:
            clade = Clade(branch_length=skbio_node.length)
            clade.clades = [convert_node(child) for child in skbio_node.children]
            return clade
    root_clade = convert_node(skbio_tree.root())
    biopython_tree = Tree(root=root_clade, rooted=True)
    return biopython_tree

def tree_construct(aln: TabularMSA) -> Tuple[DistanceMatrix, Tree]:
    """Construct a phylogenetic tree from an alignment.
    
    Uses Hamming distance to compute pairwise distances between sequences
    and applies neighbor-joining algorithm to construct the tree.
    
    Args:
        aln: Multiple sequence alignment as TabularMSA object.
        
    Returns:
        Tuple containing:
        - DistanceMatrix: Pairwise distances between sequences
        - Tree: Phylogenetic tree constructed using neighbor-joining
        
    Note:
        The function is decorated with @timeit for performance monitoring.
        Hamming distance counts the number of differing positions between sequences.
    """
    dm = DistanceMatrix.from_iterable(aln, metric=hamming, keys=aln.index)
    tree = nj(dm)
    biopython_tree = skbio_to_biopython_tree(dm,tree)
    return dm, biopython_tree

def re_clean_tree(tree: Tree) -> Tree:
    """Remove names from non-terminal nodes of a tree.
    
    Cleans up internal node names in a phylogenetic tree, leaving only
    terminal (leaf) nodes with their original names. This is often necessary
    after tree construction to ensure proper formatting.
    
    Args:
        tree: Biopython Tree object to clean.
        
    Returns:
        The same tree object with internal node names set to None.
        
    Note:
        Modifies the tree in-place by setting the name attribute of
        all non-terminal clades to None.
    """
    for clade in tree.get_nonterminals():
        clade.name = None
    return tree

def get_cluster_bio(aln, *tree_files: str, trees: List[Tree] = []) -> Tuple[List[int], int, dict]:
    """Get cluster information from alignment and trees using Biopython alignment.
    
    Assigns cluster labels to sequences based on their presence in the provided trees.
    Sequences not present in any tree are assigned to cluster 0.
    
    Args:
        aln: Biopython multiple sequence alignment object.
        *tree_files: Variable number of paths to tree files in Newick format.
        trees: Optional list of pre-loaded Tree objects (alternative to tree_files).
        
    Returns:
        Tuple containing:
        - List[int]: Cluster assignments for each sequence (0 = no cluster)
        - int: Total number of clusters found
        - dict: Mapping from cluster number to tree identifier
        
    Note:
        Uses Biopython alignment objects with .id attribute for sequence names.
        Tree files are read in Newick format and identified by their base filename.
    """
    names = [record.id for record in aln]
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

@lru_cache(maxsize=128)
def _get_indices_from_boolean_array(bool_tuple: Tuple[bool, ...]) -> Tuple[int, ...]:
    """Cached function to get indices from boolean array."""
    return tuple(np.where(np.array(bool_tuple))[0])

def sub_dm(dm: DistanceMatrix, c_list: np.ndarray) -> DistanceMatrix:
    """Optimized distance matrix subsetting with caching.
    
    Extracts a subset of the distance matrix based on boolean or integer indexing.
    Implements caching for boolean arrays to improve performance when the same
    subsetting patterns are used repeatedly.
    
    Args:
        dm: Original distance matrix to subset.
        c_list: Boolean array (True for positions to include) or
               integer array (indices of positions to include).
               
    Returns:
        New DistanceMatrix containing only the selected rows and columns.
        
    Note:
        Boolean indexing uses LRU caching to speed up repeated operations.
        The function preserves sequence IDs in the subset matrix.
    """
    if c_list.dtype == bool:
        # Use caching for frequently used boolean patterns
        bool_tuple = tuple(c_list)
        idx = np.array(_get_indices_from_boolean_array(bool_tuple))
    else:
        idx = np.asarray(c_list, dtype=np.intp) 
    
    # More efficient matrix subsetting
    sub_data = dm.data[np.ix_(idx, idx)]
    sub_ids = [dm.ids[i] for i in idx]
    
    return DistanceMatrix(sub_data, ids=sub_ids)

def sep_cluster(tree_cluster: List[int], cluster_num: int, mode: str = "ALL") -> Tuple[List[np.ndarray], List[Tuple[List[int], List[int]]]]:
    """Separate clusters into groups for comparative analysis.
    
    Creates all possible pairwise groupings of clusters for functional divergence analysis.
    Each grouping divides clusters into two sets for comparison.
    
    Args:
        tree_cluster: List of cluster assignments for each sequence (0 = unassigned).
        cluster_num: Total number of clusters to consider.
        mode: Grouping strategy:
             - "ALL": Generate all possible cluster groupings
             - "Simple": Only generate groupings where one group has size 1
             
    Returns:
        Tuple containing:
        - List[np.ndarray]: Binary cluster assignments (1 and 2) for each grouping
        - List[Tuple[List[int], List[int]]]: Original cluster groupings as integer lists
        
    Raises:
        ValueError: If mode is not "ALL" or "Simple".
        
    Note:
        The "Simple" mode is useful for one-vs-rest comparisons,
        while "ALL" mode generates all possible pairwise comparisons.
    """
    if mode not in ["ALL", "Simple"]:
        raise ValueError(f"Mode must be 'ALL' or 'Simple', got '{mode}'")
    
    group_list = get_group_list(cluster_num)
    
    # Filter group_list based on mode
    if mode == "Simple":
        # Only keep groups where one of the groups has size 1
        group_list = [group for group in group_list if len(group[0]) == 1 or len(group[1]) == 1]
    
    cluster_list = []
    for group in group_list:
        cluster = np.zeros(len(tree_cluster))
        group1, group2 = group[0], group[1]
        cluster[np.isin(tree_cluster, group1)] = 1
        cluster[np.isin(tree_cluster, group2)] = 2
        cluster_list.extend([cluster])
    return cluster_list, group_list

def get_group_list(group_num: int) -> List[Tuple[List[int], List[int]]]:
    """Generate all possible group combinations for cluster comparison.
    
    Creates all possible ways to divide clusters into two non-empty groups,
    ensuring each combination appears only once by ordering groups by size.
    
    Args:
        group_num: Total number of clusters to divide.
        
    Returns:
        List of tuples, each containing two lists representing a group division.
        The first group is always smaller than or equal in size to the second.
        
    Example:
        >>> get_group_list(3)
        [([1], [2, 3]), ([2], [1, 3]), ([3], [1, 2])]
        # For 3 clusters: three ways to divide into groups
        
    Note:
        Uses frozenset operations for efficient set arithmetic.
        Clusters are numbered starting from 1, not 0.
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

def _build_single_tree(sub_dm: DistanceMatrix, original_names: List[str]) -> Tree:
    """
    Build a single tree from distance matrix and preserve names.
    
    Args:
        sub_dm: Distance matrix subset
        original_names: Original sequence names
        
    Returns:
        Converted Biopython tree
    """
    
    # Build tree using neighbor-joining
    skbio_tree = nj(sub_dm)
    
    
    return skbio_to_biopython_tree(sub_dm,skbio_tree)

@timeit
def tree_reconstruct(dm: DistanceMatrix, cluster: np.ndarray) -> Tuple[Tree, Tree]:
    """
    Optimized tree reconstruction with parallel processing and caching.

    Args:
        dm (DistanceMatrix): The distance matrix.
        cluster (np.ndarray): Cluster assignments.

    Returns:
        Tuple[Tree, Tree]: Two reconstructed trees.

    """
    # Vectorized boolean operations - much faster than element-wise comparison
    cluster1_mask = (cluster == 1)
    cluster2_mask = (cluster == 2)
    
    # Early exit if clusters are empty
    if not np.any(cluster1_mask) or not np.any(cluster2_mask):
        raise ValueError("One or both clusters are empty")
    
    # Extract distance matrix subsets
    sub_dm1 = sub_dm(dm, cluster1_mask)
    sub_dm2 = sub_dm(dm, cluster2_mask)
    
    # Get original names lists
    original_names1 = list(sub_dm1.ids)
    original_names2 = list(sub_dm2.ids)
    tree1 = _build_single_tree(sub_dm1, original_names1)
    tree2 = _build_single_tree(sub_dm2, original_names2)
    
    return tree1, tree2


@timeit
def diverge_cal(aln_file: str, super_cluster: List[Tree], sp_type: int) -> Tuple[Optional[List[Any]], Optional[List[str]], Optional[Any]]:
    """Process trees for functional divergence analysis.
    
    Performs functional divergence analysis on a set of phylogenetic trees
    using either Gu99 (Type I) or Type2 analysis methods.
    
    Args:
        aln_file: Path to the multiple sequence alignment file.
        super_cluster: List of phylogenetic Tree objects for analysis.
        sp_type: Type of functional divergence analysis:
                - 1: Gu99 analysis (Type I functional divergence)
                - 2: Type2 analysis (Type II functional divergence)
                
    Returns:
        Tuple containing:
        - Optional[List[Any]]: Results matrix as list of lists, None if analysis failed
        - Optional[List[str]]: Position indices as strings, None if analysis failed
        - Optional[Any]: Summary DataFrame with analysis parameters, None if failed
        
    Note:
        Analysis will fail if trees don't pass pre-validation checks
        (insufficient depth or mismatched sequence names).
    """
    if pre_check(aln_file, super_cluster):
        if sp_type == 1:
            calculator = Gu99(aln_file, trees=super_cluster)
        else:
            calculator = Type2(aln_file, trees=super_cluster)
        summary = calculator.summary  # Call the method, not reference it
        position = calculator.results.index.values.tolist()
        results = calculator.results.values.tolist()
        return results, position, summary
    else:
        return None, None, None


@timeit
def diverge_cal_batch(aln_file: str, super_clusters: List[List[Tree]], sp_type: int, max_threads: Optional[int] = None) -> List[Tuple[Optional[List[Any]], Optional[List[str]], Optional[Any]]]:
    """
    Process multiple tree clusters for functional divergence analysis using parallel batch processing.
    
    This function provides significant performance improvements over sequential processing,
    especially for large numbers of clusters on multi-core systems.

    Args:
        aln_file (str): Path to the alignment file.
        super_clusters (List[List[Tree]]): List of tree clusters to process.
        sp_type (int): Type of functional divergence analysis (1 or 2).
        max_threads (Optional[int]): Maximum number of threads to use. If None, uses all available cores.

    Returns:
        List[Tuple[Optional[List[Any]], Optional[List[str]], Optional[Any]]]: 
        List of results for each cluster, where each result contains:
        - Results matrix (or None if failed)
        - Position list (or None if failed) 
        - Summary dataframe (or None if failed)
    """
    if sp_type != 1:
        # For non-Gu99 analyses, fall back to sequential processing
        print("Warning: Batch processing only available for Gu99 (sp_type=1). Using sequential processing.")
        results = []
        for super_cluster in tqdm(super_clusters, desc="Processing super cluster groups"):
            calc_result = diverge_cal(aln_file, super_cluster, sp_type)
            results.append(calc_result)
        return results
    
    # Filter valid clusters first
    valid_clusters = []
    valid_indices = []
    
    for i, super_cluster in enumerate(super_clusters):
        if pre_check(aln_file, super_cluster):
            valid_clusters.append(super_cluster)
            valid_indices.append(i)
    
    if not valid_clusters:
        print("Warning: No valid clusters found for batch processing.")
        return [(None, None, None)] * len(super_clusters)
    
    print(f"Starting parallel batch processing of {len(valid_clusters)} valid clusters "
          f"(out of {len(super_clusters)} total) using Gu99Batch...")
    
    # Create batch calculator
    batch = Gu99Batch(max_threads=max_threads)
    
    # Add tasks for valid clusters
    for i, super_cluster in enumerate(valid_clusters):
        batch.add_task(
            aln_file,
            trees=super_cluster,
            task_name=f"SuperCluster_{valid_indices[i]}"
        )
    
    # Execute batch calculation
    batch.calculate_batch()
    
    # Get results
    batch_results = batch.get_results()
    
    # Map results back to original order
    final_results = [(None, None, None)] * len(super_clusters)
    
    for i, result in enumerate(batch_results):
        original_index = valid_indices[i]
        
        if result['success']:
            # Extract data from batch result
            summary_df = result.get('summary')
            results_df = result.get('results')
            
            if summary_df is not None and results_df is not None:
                position = results_df.index.values.tolist()
                results_values = results_df.values.tolist()
                final_results[original_index] = (results_values, position, summary_df)
            else:
                print(f"Warning: Missing data in successful result for cluster {original_index}")
                final_results[original_index] = (None, None, None)
        else:
            print(f"Warning: Task failed for cluster {original_index}: {result['error_message']}")
            final_results[original_index] = (None, None, None)
    
    print(f"Batch processing completed: {len([r for r in final_results if r[0] is not None])}/{len(super_clusters)} clusters successful")
    
    return final_results

# Helper to enable process-parallel NJ + cleanup
# using raw NumPy arrays to leverage joblib memmapping
from typing import cast

def _reconstruct_and_clean_from_array(dm_data: np.ndarray, dm_ids: List[str], cluster: np.ndarray) -> List[Tree]:
    """Reconstruct phylogenetic trees from distance matrix data and cluster assignments.
    
    Helper function that builds two phylogenetic trees from distance matrix subsets
    corresponding to cluster assignments. Handles edge cases where clusters might be empty.
    
    Args:
        dm_data: 2D numpy array containing pairwise distances.
        dm_ids: List of sequence identifiers corresponding to matrix rows/columns.
        cluster: Array with cluster assignments (1 and 2) for each sequence.
        
    Returns:
        List containing two Tree objects:
        - First tree: sequences assigned to cluster 1
        - Second tree: sequences assigned to cluster 2
        Empty trees are returned if any cluster is empty.
        
    Note:
        Uses neighbor-joining algorithm for tree construction.
        Trees are cleaned to remove internal node names.
        Designed for parallel processing with joblib.
    """
    # Boolean masks for two clusters
    cluster1_mask = (cluster == 1)
    cluster2_mask = (cluster == 2)

    # Indices for submatrices
    idx1 = np.where(cluster1_mask)[0]
    idx2 = np.where(cluster2_mask)[0]

    # Early exit if any subgroup is empty
    if idx1.size == 0 or idx2.size == 0:
        # Build empty trees to keep pipeline robust
        empty_tree1 = Tree(root=Clade())
        empty_tree2 = Tree(root=Clade())
        return [empty_tree1, empty_tree2]

    sub1 = dm_data[np.ix_(idx1, idx1)]
    sub2 = dm_data[np.ix_(idx2, idx2)]
    ids1 = [dm_ids[i] for i in idx1]
    ids2 = [dm_ids[i] for i in idx2]

    dm1 = DistanceMatrix(sub1, ids=ids1)
    dm2 = DistanceMatrix(sub2, ids=ids2)

    # Run NJ
    t1 = nj(dm1)
    t2 = nj(dm2)

    # Preserve names and convert to Biopython trees
    t1 = tree_reconstruct.__defaults__[0](t1, ids1) if False else t1 
    bi1 = skbio_to_biopython_tree(dm1,t1)
    bi2 = skbio_to_biopython_tree(dm2,t2)
    bi1 = re_clean_tree(bi1)
    bi2 = re_clean_tree(bi2)
    return [bi1, bi2]

@timeit
def get_super_cluster_pp(aln_file: str, *tree_files: str, sp_type: int = 1, trees: List[Tree] = [], verbose: bool = True, mode: str = "ALL", metric: str = 'hamming', parallel: bool = True, max_threads: Optional[int] = None, conswins: Optional[Dict] = None) -> Tuple[np.ndarray, List[List[Tree]], List[Tuple[List[int], List[int]]], dict, List[str], List[Any]]:
    """
    Perform super cluster analysis for functional divergence with optional parallel processing.

    Args:
        aln_file (str): Path to the alignment file.
        *tree_files (str): Paths to tree files.
        sp_type (int, optional): Type of functional divergence analysis. Defaults to 1.
        trees (List[Tree], optional): Pre-loaded trees. Defaults to [].
        verbose (bool, optional): Whether to print progress. Defaults to True.
        mode (str, optional): Mode for group selection. "ALL" returns all groups, 
                             "Simple" returns only groups where one group has size 1. Defaults to "ALL".
        metric (str, optional): Distance metric for tree construction. Defaults to 'hamming'.
        parallel (bool, optional): Whether to use parallel batch processing for Gu99 analyses. Defaults to True.
        max_threads (Optional[int], optional): Maximum number of threads for parallel processing. 
                                              If None, uses all available cores. Defaults to None.
        conswins (Optional[Dict], optional): conswins like groupsim.
                                           Expected keys: 'cons_win_len' (int), 'lambda_param' (float).
                                           Defaults to None (no conservation weighting).

    Returns:
        Tuple[np.ndarray, List[List[Tree]], List[Tuple[List[int], List[int]]], dict, List[str], List[Any]]:
        Results array, super cluster list, group list, tree dictionary, position list, and summary list.
    """
    if metric == 'hamming':
    
        aln = sk_aln_read(aln_file)
        tree_cluster, cluster_num, tree_dict = get_cluster(aln, *tree_files, trees=trees)
        cluster_list, group_list = sep_cluster(tree_cluster, cluster_num, mode=mode)
        dm, _ = tree_construct(aln)
    else:
        aln = AlignIO.read(aln_file, 'fasta')
        tree_cluster,cluster_num,tree_dict = get_cluster_bio(aln,*tree_files,trees=trees)
        cluster_list, group_list = sep_cluster(tree_cluster, cluster_num, mode=mode)
        
        calculator = DistanceCalculator(metric)
        dm_bio = calculator.get_distance(aln)
        n = len(dm_bio.matrix)
        matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1): 
                val = dm_bio.matrix[i][j] if j < len(dm_bio.matrix[i]) else dm_bio.matrix[j][i]
                matrix[i, j] = matrix[j, i] = val
        dm = DistanceMatrix(matrix, ids=[record.id for record in aln])
    dm_data: np.ndarray = dm.data
    dm_ids: List[str] = list(dm.ids)


    super_clusters=[]
    for cluster in cluster_list:
        super_cluster = _reconstruct_and_clean_from_array(dm_data,dm_ids,cluster)
        super_clusters.append(super_cluster)

    # Choose processing strategy based on parameters
    if parallel and sp_type == 1:
        # Use parallel batch processing for Gu99 analyses
        calc_results_all = diverge_cal_batch(aln_file, super_clusters, sp_type, max_threads)
    else:
        # Use sequential processing
        if parallel and sp_type != 1:
            print(f"Warning: Parallel processing not available for sp_type={sp_type}. Using sequential processing.")
        
        calc_results_all = []
        for super_cluster in tqdm(super_clusters, desc="Processing super cluster groups"):
            calc_results = diverge_cal(aln_file, super_cluster, sp_type)
            calc_results_all.append(calc_results)
    super_cluster_list = []
    results_list = []
    position_list: List[str] = []
    summary_list = []
    new_group_list: List[Tuple[List[int], List[int]]] = []
    
    # Aggregate
    if calc_results_all:
        for i, (super_cluster, (calc_results, position, summary)) in enumerate(zip(super_clusters, calc_results_all)):
            if calc_results is not None:
                super_cluster_list.append(super_cluster)
                new_group_list.append(group_list[i])
                results_list.append(calc_results)
                summary_list.append(summary)
                if not position_list and position is not None:
                    position_list = position
    
    # Handle empty results
    if not results_list:
        results_array = np.array([])
    else:
        results_array = np.reshape(np.array(results_list), (len(super_cluster_list), -1))
    
    # Apply conservation window weighting if requested
    if conswins is not None and results_array.size > 0:
        # Read alignment for conservation analysis
        aln = AlignIO.read(aln_file, 'fasta')
        alignment_seqs = [str(record.seq) for record in aln]
        
        cons_win_len = conswins.get('cons_win_len', 3)
        lambda_param = conswins.get('lambda_param', 0.7)
        
        # Apply conservation weighting to each row (cluster) of results
        weighted_results = []
        for i, result_row in enumerate(results_array):
            # Convert to list of Optional[float] for conservation_window_score
            scores = [float(score) if not np.isnan(score) else None for score in result_row]
            
            # Apply conservation window scoring
            weighted_scores = conservation_window_score(
                scores, alignment_seqs, cons_win_len, lambda_param
            )
            
            # Convert back to numpy array, replacing None with NaN
            weighted_row = np.array([score if score is not None else np.nan 
                                   for score in weighted_scores])
            weighted_results.append(weighted_row)
        
        results_array = np.array(weighted_results)
    
    return results_array, super_cluster_list, new_group_list, tree_dict, position_list, summary_list

class SuperCluster:
    """Class for performing super cluster analysis of functional divergence."""

    def __init__(self, aln_file: str, *tree_files: str, trees: List[Tree] = [], sp_type: int = 1, verbose: bool = True, mode: str = "ALL", filter: bool = False, metric: str = 'hamming', parallel: bool = True, max_threads: Optional[int] = None, conswins: Optional[Dict] = None):
        """
        Initialize SuperCluster with optional parallel processing support.

        Args:
            aln_file (str): Path to the alignment file.
            *tree_files (str): Paths to tree files.
            trees (List[Tree], optional): Pre-loaded trees. Defaults to [].
            sp_type (int, optional): Type of functional divergence analysis. Defaults to 1.
            verbose (bool, optional): Whether to print progress. Defaults to True.
            mode (str, optional): Mode for group selection. "ALL" returns all groups, 
                                 "Simple" returns only groups where one group has size 1. Defaults to "ALL".
            filter (bool, optional): Whether to apply filtering to results. Defaults to False.
            metric (str, optional): Distance metric for tree construction. Defaults to 'hamming'.
            parallel (bool, optional): Whether to use parallel batch processing for Gu99 analyses. Defaults to True.
            max_threads (Optional[int], optional): Maximum number of threads for parallel processing. 
                                                  If None, uses all available cores. Defaults to None.
            conswins (Optional[Dict], optional): ConservationWindow parameters for weighting scores.
                                               Expected keys: 'cons_win_len' (int), 'lambda_param' (float).
                                               Defaults to None (no conservation weighting).
        """
        self.aln_file = aln_file
        self.tree_files = tree_files
        self.filter = filter
        self.parallel = parallel
        self.max_threads = max_threads
        self.conswins = conswins
        
        # Perform super cluster analysis with optional parallel processing
        self.pp_list, self.tree_list, self.group_list, self.tree_dict, self.position_list, self.summary_list = get_super_cluster_pp(
            aln_file, *tree_files, 
            sp_type=sp_type, 
            trees=trees, 
            verbose=verbose, 
            mode=mode, 
            metric=metric,
            parallel=parallel,
            max_threads=max_threads,
            conswins=conswins
        )
        self.summary = self.get_summary()
        
    def get_summary(self):
        """Generate summary statistics for super cluster analysis.
        
        Creates a comprehensive summary DataFrame containing:
        - Original summary statistics from each cluster group analysis
        - Count statistics for different Qk threshold levels (>0.5, >0.67, >0.9)
        
        Returns:
            DataFrame with cluster groups as columns and parameters/statistics as rows.
            Includes both original analysis parameters and threshold-based counts.
            
        Note:
            Qk represents posterior probability of functional divergence at each site.
            Higher Qk values indicate stronger evidence for functional divergence.
        """
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
    
    @property
    def results(self):
        results = pd.DataFrame(self.pp_list, index=[f"{i}" for i in self.group_list], columns=self.position_list)
        if self.filter:
            from .filter import filter_results
            aln = AlignIO.read(self.aln_file, 'fasta')
            trees = [Phylo.read(tree_file, 'newick') for tree_file in self.tree_files]
            results, _ = filter_results(results, aln, trees)
            return results
        else:
            return results
    
    def benchmark_parallel_vs_sequential(self, runs: int = 1) -> Dict[str, float]:
        """Benchmark parallel vs sequential processing performance.
        
        Compares the execution time of parallel batch processing against
        sequential processing for functional divergence analysis.
        
        Args:
            runs: Number of benchmark runs to average (default: 1).
                 More runs provide more reliable timing estimates.
                 
        Returns:
            Dictionary containing timing results and performance metrics:
            - 'sequential_time': Average time for sequential processing (seconds)
            - 'parallel_time': Average time for parallel processing (seconds) 
            - 'speedup': Ratio of sequential/parallel times
            - 'sequential_times': List of individual sequential run times
            - 'parallel_times': List of individual parallel run times
            
        Raises:
            ValueError: If no tree files are available for benchmarking.
            
        Note:
            Only benchmarks Gu99 analyses (sp_type=1) as parallel processing
            is currently only implemented for this analysis type.
            Prints detailed timing results and speedup information.
        """
        if not hasattr(self, 'tree_files') or not self.tree_files:
            raise ValueError("Cannot benchmark: no tree files available")
            
        print(f"Benchmarking parallel vs sequential processing ({runs} run{'s' if runs > 1 else ''})...")
        
        # Sequential timing
        sequential_times = []
        for run in range(runs):
            print(f"Sequential run {run + 1}/{runs}...")
            start_time = time.time()
            
            # Run sequential analysis
            get_super_cluster_pp(
                self.aln_file, *self.tree_files,
                sp_type=1,  # Only benchmark Gu99
                verbose=False,
                parallel=False
            )
            
            sequential_time = time.time() - start_time
            sequential_times.append(sequential_time)
        
        # Parallel timing
        parallel_times = []
        for run in range(runs):
            print(f"Parallel run {run + 1}/{runs}...")
            start_time = time.time()
            
            # Run parallel analysis
            get_super_cluster_pp(
                self.aln_file, *self.tree_files,
                sp_type=1,  # Only benchmark Gu99
                verbose=False,
                parallel=True,
                max_threads=self.max_threads
            )
            
            parallel_time = time.time() - start_time
            parallel_times.append(parallel_time)
        
        # Calculate averages
        avg_sequential = np.mean(sequential_times)
        avg_parallel = np.mean(parallel_times)
        speedup = avg_sequential / avg_parallel
        
        results = {
            'sequential_time': avg_sequential,
            'parallel_time': avg_parallel,
            'speedup': speedup,
            'sequential_times': sequential_times,
            'parallel_times': parallel_times
        }
        
        print(f"\n=== Benchmark Results ===")
        print(f"Average Sequential Time: {avg_sequential:.2f} seconds")
        print(f"Average Parallel Time: {avg_parallel:.2f} seconds")
        print(f"Speedup: {speedup:.2f}x")
        
        if speedup > 1:
            print(f"Parallel processing is {speedup:.2f}x faster!")
        else:
            print(f"Sequential processing is {1/speedup:.2f}x faster (parallel overhead)")
            
        return results
        
# Main execution
if __name__ == "__main__":
    aln_file = "E:/verysync/diverge_pybind/web/statics/example_data/ERBB_family.fas"
    tree_files = [
        "E:/verysync/diverge_pybind/web/statics/example_data/EGFR.tree",
        "E:/verysync/diverge_pybind/web/statics/example_data/ERBB2.tree",
        "E:/verysync/diverge_pybind/web/statics/example_data/ERBB3.tree",
        "E:/verysync/diverge_pybind/web/statics/example_data/ERBB4.tree"
    ]
    # Example without conservation weighting
    super_cluster = SuperCluster(aln_file, *tree_files)
    
    # Example with conservation weighting (ConservationWindow parameters)
    conswins_params = {
        'cons_win_len': 3,      # Conservation window size 
        'lambda_param': 0.7     # Lambda for linear combination (0.7 = 70% original score + 30% conservation)
    }
    # super_cluster_with_conswins = SuperCluster(aln_file, *tree_files, conswins=conswins_params)