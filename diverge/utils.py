from typing import List, Optional, Dict, Any, Tuple, Generator
import os
import numpy as np
from Bio import Phylo
from tqdm import tqdm
import time
import contextlib
import joblib
from tqdm.autonotebook import tqdm

from .binding import (
    Gu99, Gu2001, Type2, Fdr, Effective, Asym, TypeOneAnalysis, Rvs
)

# Analysis registry: centralized configuration for all analyses
ANALYSIS_REGISTRY = {
    'Gu99': {
        'class': Gu99,
        'enabled': lambda ctx: True,
        'result_attrs': [('results', 'gu99_results'), ('summary', 'gu99_summary')],
        'description': 'Gu99 calculation',
        'run_fundist': True,
    },
    'Gu2001': {
        'class': Gu2001,
        'enabled': lambda ctx: ctx['has_branch_length'],
        'result_attrs': [('results', 'gu2001_results'), ('summary', 'gu2001_summary')],
        'description': 'Gu2001 distance calculation',
        'run_fundist': False,
    },
    'Type1Analysis': {
        'class': TypeOneAnalysis,
        'enabled': lambda ctx: ctx['cluster_num'] == 3,
        'result_attrs': [('summary', 'type1analysis_summary'), ('results', 'type1analysis_results')],
        'description': 'Type one Analysis distance calculation',
        'run_fundist': False,
    },
    'Type2': {
        'class': Type2,
        'enabled': lambda ctx: True,
        'result_attrs': [('summary', 'type2_summary'), ('results', 'type2_results')],
        'description': 'Type two calculation',
        'run_fundist': False,
    },
    'Effective': {
        'class': Effective,
        'enabled': lambda ctx: ctx['cluster_num'] == 2,
        'result_attrs': [('type1_results', 'type1_effective'), ('type2_results', 'type2_effective')],
        'description': 'Effective number of site calculation',
        'run_fundist': False,
    },
    'Fdr': {
        'class': Fdr,
        'enabled': lambda ctx: True,
        'result_attrs': [('type1_results', 'type1_fdr'), ('type2_results', 'type2_fdr')],
        'description': 'Fdr calculation',
        'run_fundist': False,
    },
    'FunDist': {
        'class': None,  # Special case: called on Gu99 instance
        'enabled': lambda ctx: ctx['cluster_num'] >= 3,
        'result_attrs': [],
        'description': 'Function distance calculation',
        'run_fundist': False,
    },
    'Asym': {
        'class': Asym,
        'enabled': lambda ctx: ctx['cluster_num'] == 3,
        'result_attrs': [('results', 'asym_results')],
        'description': 'Asym test',
        'run_fundist': False,
    },
}

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
        self.tree_files = tree_files
        self.aln_file = aln_file

        # Build execution context for configuration-driven pipeline
        self._context = {
            'cluster_num': len(tree_files),
            'has_branch_length': _has_branch_length(*tree_files)
        }

        # Determine which analyses are enabled based on context
        self.func_dict = {
            name: config['enabled'](self._context)
            for name, config in ANALYSIS_REGISTRY.items()
        }

        # Store results in unified dict
        self._results = {}

        # Execute pipeline
        self._pipeline()

        # Maintain backward compatibility
        self.detail = [name for name, enabled in self.func_dict.items() if enabled]

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

    def _pipeline(self) -> None:
        """Execute the calculation pipeline using configuration-driven approach."""
        progress_bar = tqdm(total=sum(self.func_dict.values()), desc="Calculate pipeline")

        # Execute analyses in order defined by ANALYSIS_REGISTRY
        for name, config in ANALYSIS_REGISTRY.items():
            if self.func_dict.get(name, False):
                try:
                    self._execute_analysis(name, config, progress_bar)
                except Exception as e:
                    printv(f"Error occurred during {name} calculation: {e}")
                    raise

        progress_bar.close()

    def _execute_analysis(self, name: str, config: Dict[str, Any], progress_bar: tqdm) -> None:
        """
        Execute a single analysis and store results.

        Args:
            name (str): Analysis name
            config (dict): Analysis configuration from ANALYSIS_REGISTRY
            progress_bar (tqdm): Progress bar instance
        """
        # Special case: FunDist is called on Gu99 instance
        if name == 'FunDist':
            if hasattr(self, '_gu99_instance'):
                self.fundist_results = self._gu99_instance.fundist()
                self._results['fundist_results'] = self.fundist_results
            progress_bar.update(1)
            progress_bar.set_description(f"{config['description']} running...")
            progress_bar.refresh()
            return

        # Standard analysis execution
        analyzer = config['class'](self.aln_file, *self.tree_files)

        # Store instance for FunDist (Gu99 special case)
        if name == 'Gu99':
            self._gu99_instance = analyzer

        # Map results from analyzer to instance attributes
        for source_attr, target_attr in config['result_attrs']:
            value = getattr(analyzer, source_attr)
            setattr(self, target_attr, value)
            self._results[target_attr] = value

        # Update progress
        progress_bar.update(1)
        progress_bar.set_description(f"{config['description']} running...")
        progress_bar.refresh()

    @property
    def result_summary(self) -> Dict[str, Any]:
        """
        Get unified result summary.

        Returns:
            Dict[str, Any]: Dictionary containing all analysis results
        """
        return self._results

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