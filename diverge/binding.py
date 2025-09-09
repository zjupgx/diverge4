from typing import List, Optional, Dict, Any, Tuple, Union, Callable
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Phylo,AlignIO
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
    
    # If all formats fail
    raise ValueError(f"Tree file '{tree_file}' could not be parsed: {last_error}")

def get_colnames(r_names: List[str]) -> List[str]:
    """Get column names from r_names.
    
    Creates a copy of the input list to avoid modifying the original.
    
    Args:
        r_names: List of column names to copy.
        
    Returns:
        A new list containing the same column names.
    """
    return r_names.copy()

def load_tree_file(tree_file: str, check: bool = True) -> str:
    """Load a tree file and optionally check its validity.
    
    Reads a tree file, validates it if requested, and returns it in Newick format.
    
    Args:
        tree_file: Path to the tree file to load.
        check: Whether to validate the tree structure (depth >= 3).
        
    Returns:
        The tree in Newick format string, with trailing newlines removed.
        
    Raises:
        FileNotFoundError: If the tree file doesn't exist.
        ValueError: If the tree is invalid (when check=True).
    """
    tree = read_tree(tree_file)
    if check:
        check_tree(tree)
    return tree.format("newick").strip("\n")

def check_tree_depth(tree: Tree, min_depth: int = 3) -> bool:
    """Check if tree has sufficient depth.
    
    Calculates the maximum path length from root to any terminal node
    and compares it against the minimum required depth.
    
    Args:
        tree: Biopython Tree object to check.
        min_depth: Minimum required depth (default: 3).
        
    Returns:
        True if tree depth exceeds min_depth, False otherwise.
    """
    tree_depth = max(len(tree.trace(tree.root, clade)) for clade in tree.get_terminals())
    return tree_depth > min_depth

def check_tree(tree: Tree) -> bool:
    """Check the validity of a single tree.
    
    Validates that a tree has sufficient depth for phylogenetic analysis.
    Trees with depth <= 3 are considered invalid for functional divergence analysis.
    
    Args:
        tree: Biopython Tree object to validate.
        
    Returns:
        True if the tree is valid.
        
    Raises:
        ValueError: If the tree depth is less than or equal to 3.
    """
    if not check_tree_depth(tree):
        raise ValueError(f'Tree depth is less than 3, please check your tree: {tree}')
    return True

def check_tree_file(*tree_files: str) -> bool:
    """Check the validity of multiple tree files.
    
    Validates that all provided tree files contain trees with sufficient depth
    for phylogenetic analysis.
    
    Args:
        *tree_files: Variable number of paths to tree files.
        
    Returns:
        True if all trees are valid.
        
    Raises:
        ValueError: If any tree has insufficient depth (<=3).
        FileNotFoundError: If any tree file doesn't exist.
    """
    for tree_file in tree_files:
        tree = read_tree(tree_file)
        if not check_tree_depth(tree):
            raise ValueError(f'Tree depth is less than 3, please check your tree file: {tree_file}')
    return True

def fcheck(*files: str) -> None:
    """Check if files exist.
    
    Validates that all specified files exist on the filesystem.
    Useful for batch validation of input files before processing.
    
    Args:
        *files: Variable number of file paths to check.
        
    Raises:
        FileNotFoundError: If any of the specified files don't exist.
    """
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
        fcheck(aln_file, *tree_files)
        
        self.calculator_module = calculator_module
        
        self.input = [aln_file]
        self._prepare_trees(tree_files, trees, check_trees)
        self._prepare_cluster_names(cluster_name)
        
        self._calculate()
        
        self._summary_cache = None
        self._results_cache = None
    
    def _prepare_trees(self, tree_files: Tuple[str, ...], trees: List[Tree], check_trees: bool) -> None:
        """Prepare tree inputs for analysis.
        
        Processes either tree files or pre-loaded Tree objects and converts them
        to Newick format strings for the calculator. Validates trees if requested.
        
        Args:
            tree_files: Tuple of paths to tree files.
            trees: List of pre-loaded Biopython Tree objects.
            check_trees: Whether to validate tree structure before processing.
            
        Side Effects:
            Appends Newick format strings to self.input list.
        """
        if trees:
            for tree in trees:
                if check_trees:
                    check_tree(tree)
                self.input.append(tree.format("newick").strip("\n"))
        else:
            for tree_file in tree_files:
                self.input.append(load_tree_file(tree_file, check=check_trees))
    
    def _prepare_cluster_names(self, cluster_name: Optional[List[str]]) -> None:
        """Set up cluster names for analysis.
        
        Creates cluster names based on input or generates default names.
        Default names follow the pattern 'cluster_1', 'cluster_2', etc.
        
        Args:
            cluster_name: Optional list of custom cluster names. If None,
                         default names will be generated.
                         
        Side Effects:
            Sets self.cluster_name attribute.
            
        Raises:
            TypeError: If cluster_name is not None or a list.
        """
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
        """Print help information for the calculator module.
        
        Displays the built-in help documentation for the C++ calculator module
        associated with this analysis type. Useful for understanding available
        methods and parameters.
        """
        if hasattr(self, 'calculator_module'):
            print(help(self.calculator_module))
    @property
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
    
    @property
    def results(self) -> pd.DataFrame:
        """Generate detailed results with caching.
        
        Returns:
            DataFrame with detailed analysis results. The index represents 
            sequence positions and is 0-based (starting from 0).
        """
        if self._results_cache is not None:
            return self._results_cache
            
        if not hasattr(self.calculator, '_results'):
            return None
            
        columns_names = get_colnames(self.calculator._r_names())
        
        if hasattr(self.calculator, '_kept'):
            index = [i for i in self.calculator._kept()]
        else:
            index = range(0, np.size(self.calculator._results(), 0))
            
        results = pd.DataFrame(
            self.calculator._results(),
            columns=columns_names,
            index=index
        )
        results.index.name = self.get_position_label()
        
        self._results_cache = results
        return results
    
    def get_position_label(self) -> str:
        """Get the label for the position index.
        
        Returns:
            Label for position index. Positions are 0-based.
        """
        return "Position"
    
    def clear_cache(self) -> None:
        """Clear cached results.
        
        Removes cached summary and results DataFrames to free memory
        or force recalculation. Useful when underlying data has changed.
        
        Side Effects:
            Resets _summary_cache and _results_cache to None.
        """
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
        """Set index name for dataframe.
        
        Helper method to set the index name of a DataFrame if it doesn't already have one.
        
        Args:
            df: DataFrame to modify (can be None).
            name: Name to set for the index.
            
        Returns:
            The same DataFrame with index name set, or None if input was None.
        """
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
        filter:bool = False,
        trees: List[Tree] = []
    ) -> None:
        """Initialize Gu99 analysis."""
        super().__init__(aln_file, *tree_files, cluster_name=cluster_name, trees=trees, calculator_module=_gu99cpp)
        self._summary_result = self.summary
        self.filter = filter
        self.aln_file = aln_file
        self.tree_files = tree_files
    @property
    def summary(self) -> pd.DataFrame:
        """Generate a summary of the results."""
        columns_names = get_colnames(self.calculator._r_names())
        summary_results = pd.DataFrame(columns=columns_names)
        for dict_item in self.calculator._summary():
            summary_results.loc[dict_item["name"]] = dict_item["values"]
        return summary_results
    @property
    def results(self) -> pd.DataFrame:
        results = super().results
        if self.filter:
            from .filter import filter_results
            aln = AlignIO.read(self.aln_file,"fasta")
            trees = [Phylo.read(tree_file,"newick") for tree_file in self.tree_files]
            results,_ = filter_results(results.T,aln,trees)
            results = results.T
        return results
    def _calculate(self) -> None:
        """Create a new Gu99 Calculator and complete the calculation steps."""
        self.calculator = self.calculator_module.create_calculator(self.input, self.cluster_name)
        self.calculator.calculate()
    
    def fundist(self) -> Optional[pd.DataFrame]:
        """Calculate functional distance if applicable.
        
        Computes functional distances between clusters based on theta values
        from the Gu99 analysis. Requires at least 3 clusters to perform the calculation.
        
        The method calculates pairwise distances using -log(1 - theta) transformation
        and applies a branch-length correction algorithm to estimate functional
        distances between clusters.
        
        Returns:
            DataFrame with functional distances for each cluster, or None if
            insufficient clusters (<=2) are available.
            
        Note:
            This calculation is only meaningful when there are at least 3 clusters
            in the analysis.
        """
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

class Gu99Batch:
    """
    Parallel batch processing class for multiple Gu99 analyses.
    
    This class enables true parallel processing of multiple Gu99 tasks
    without being limited by Python's GIL, providing near-linear speedup
    with the number of CPU cores.
    
    Features:
    - True parallel processing in C++
    - Memory efficient with isolated task spaces
    - Fault tolerance - single task failures don't affect others
    - Scalable to multi-core systems
    """
    
    def __init__(self, max_threads: Optional[int] = None) -> None:
        """
        Initialize batch calculator.
        
        Args:
            max_threads: Maximum number of threads to use. 
                        If None, uses all available CPU cores.
        """
        self.max_threads = max_threads if max_threads is not None else 0
        self.batch_calculator = _gu99cpp.create_batch_calculator(self.max_threads)
        self.tasks: List[Dict[str, Any]] = []
        self._results_cache: Optional[List[Dict[str, Any]]] = None
        
    def add_task(
        self,
        aln_file: str,
        *tree_files: str,
        cluster_name: Optional[List[str]] = None,
        trees: List[Tree] = [],
        check_trees: bool = True,
        task_name: Optional[str] = None
    ) -> int:
        """
        Add a task to the batch queue.
        
        Args:
            aln_file: Path to alignment file
            tree_files: Paths to tree files
            cluster_name: Names of clusters
            trees: List of Biopython.Phylo tree objects (alternative to tree_files)
            check_trees: Whether to validate trees
            task_name: Optional name for this task
            
        Returns:
            Task ID (0-based index)
        """
        # Validate inputs
        fcheck(aln_file, *tree_files)
        
        if trees and tree_files:
            raise ValueError("Cannot specify both tree_files and trees")
        
        # Prepare input arguments
        input_args = [aln_file]
        
        # Handle tree inputs
        if trees:
            for tree in trees:
                if check_trees:
                    check_tree(tree)
                input_args.append(tree.format("newick").strip("\n"))
        else:
            for tree_file in tree_files:
                input_args.append(load_tree_file(tree_file, check=check_trees))
        
        # Prepare cluster names
        if cluster_name is None:
            cluster_names = [f"cluster_{i}" for i in range(1, len(input_args))]
        elif isinstance(cluster_name, list):
            cluster_names = cluster_name
        else:
            raise TypeError("cluster_name must be a list")
        
        # Store task metadata
        task_id = len(self.tasks)
        task_info = {
            'task_id': task_id,
            'task_name': task_name or f"Task_{task_id}",
            'aln_file': aln_file,
            'tree_files': tree_files,
            'cluster_name': cluster_names,
            'input_args': input_args
        }
        self.tasks.append(task_info)
        
        # Add to C++ batch calculator
        self.batch_calculator.add_task(input_args, cluster_names)
        
        return task_id
    
    def calculate_batch(self) -> None:
        """
        Execute all queued tasks in parallel.
        
        This method performs true parallel computation in C++,
        bypassing Python's GIL limitations.
        """
        if not self.tasks:
            raise ValueError("No tasks have been added to the batch")
        
        print(f"Starting batch calculation of {len(self.tasks)} tasks "
              f"using up to {self.batch_calculator.get_max_threads()} threads...")
        
        # Clear any cached results
        self._results_cache = None
        
        # Execute batch calculation
        self.batch_calculator.calculate_batch()
        
        # Check for any failures
        success_status = self.batch_calculator.get_success_status()
        error_messages = self.batch_calculator.get_error_messages()
        
        failed_tasks = []
        for i, (success, error_msg) in enumerate(zip(success_status, error_messages)):
            if not success:
                failed_tasks.append((i, self.tasks[i]['task_name'], error_msg))
        
        if failed_tasks:
            print(f"Warning: {len(failed_tasks)} tasks failed:")
            for task_id, task_name, error_msg in failed_tasks:
                print(f"  - {task_name} (ID: {task_id}): {error_msg}")
        
        successful_tasks = sum(success_status)
        print(f"Batch calculation completed: {successful_tasks}/{len(self.tasks)} tasks successful")
    
    def get_results(self) -> List[Dict[str, Any]]:
        """
        Get results from all completed tasks.
        
        Returns:
            List of dictionaries containing results for each task.
            Each dictionary contains:
            - task_id: Task identifier
            - task_name: Task name
            - success: Whether the task completed successfully
            - error_message: Error message if failed
            - summary: Summary DataFrame (if successful)
            - results: Results DataFrame (if successful)
        """
        if self._results_cache is not None:
            return self._results_cache
        
        raw_results = self.batch_calculator.get_all_results()
        processed_results = []
        
        for i, raw_result in enumerate(raw_results):
            result_dict = {
                'task_id': raw_result['task_id'],
                'task_name': self.tasks[i]['task_name'],
                'success': raw_result['success'],
                'error_message': raw_result['error_message']
            }
            
            if raw_result['success']:
                # Process summary
                if 'summary' in raw_result:
                    r_names = raw_result['r_names'] if 'r_names' in raw_result else []
                    
                    # Collect summary data first
                    summary_data = {}
                    for summary_item in raw_result['summary']:
                        summary_data[summary_item["name"]] = summary_item["values"]
                    
                    if summary_data:
                        # Determine number of columns from first item
                        first_values = next(iter(summary_data.values()))
                        n_cols = len(first_values) if hasattr(first_values, '__len__') else 1
                        
                        # Generate column names if not provided
                        if not r_names:
                            r_names = [f"cluster_{j}" for j in range(1, n_cols + 1)]
                        
                        # Create DataFrame with proper structure
                        summary_df = pd.DataFrame(index=summary_data.keys(), columns=r_names)
                        
                        # Fill the DataFrame
                        for param_name, values in summary_data.items():
                            if hasattr(values, '__len__') and len(values) > 1:
                                summary_df.loc[param_name] = values
                            else:
                                # Handle single values or arrays with one element
                                val = values[0] if hasattr(values, '__len__') and len(values) == 1 else values
                                summary_df.loc[param_name] = [val] * n_cols
                        
                        summary_df.index.name = "Parameters"
                        result_dict['summary'] = summary_df
                
                # Process detailed results
                if 'results' in raw_result and 'kept' in raw_result:
                    results_array = raw_result['results']
                    kept_positions = raw_result['kept']
                    r_names = raw_result['r_names']
                    
                    results_df = pd.DataFrame(
                        results_array,
                        columns=r_names,
                        index=kept_positions
                    )
                    results_df.index.name = "Position"
                    result_dict['results'] = results_df
            
            processed_results.append(result_dict)
        
        self._results_cache = processed_results
        return processed_results
    
    def get_successful_results(self) -> List[Dict[str, Any]]:
        """Get results only from successfully completed tasks."""
        all_results = self.get_results()
        return [result for result in all_results if result['success']]
    
    def get_failed_tasks(self) -> List[Dict[str, Any]]:
        """Get information about failed tasks."""
        all_results = self.get_results()
        return [result for result in all_results if not result['success']]
    
    def get_task_result(self, task_id: int) -> Dict[str, Any]:
        """
        Get result for a specific task.
        
        Args:
            task_id: Task identifier
            
        Returns:
            Dictionary containing task result
        """
        results = self.get_results()
        if 0 <= task_id < len(results):
            return results[task_id]
        else:
            raise IndexError(f"Task ID {task_id} out of range")
    
    def clear_tasks(self) -> None:
        """Clear all tasks and results."""
        self.batch_calculator.clear_tasks()
        self.tasks.clear()
        self._results_cache = None
    
    @property
    def num_tasks(self) -> int:
        """Get number of tasks in the batch."""
        return len(self.tasks)
    
    @property
    def max_threads_available(self) -> int:
        """Get maximum number of threads available."""
        return self.batch_calculator.get_max_threads()
    
    def print_summary(self) -> None:
        """Print a summary of batch processing results."""
        if not self.tasks:
            print("No tasks in batch")
            return
        
        results = self.get_results()
        successful = sum(1 for r in results if r['success'])
        failed = len(results) - successful
        
        print(f"\n=== Batch Processing Summary ===")
        print(f"Total tasks: {len(results)}")
        print(f"Successful: {successful}")
        print(f"Failed: {failed}")
        print(f"Max threads used: {self.max_threads_available}")
        
        if failed > 0:
            print(f"\nFailed tasks:")
            for result in results:
                if not result['success']:
                    print(f"  - {result['task_name']}: {result['error_message']}")
    
    def save_results(self, output_dir: str, format: str = 'csv') -> None:
        """
        Save all successful results to files.
        
        Args:
            output_dir: Directory to save results
            format: Output format ('csv', 'xlsx', 'pickle')
        """
        import os
        
        os.makedirs(output_dir, exist_ok=True)
        successful_results = self.get_successful_results()
        
        if not successful_results:
            print("No successful results to save")
            return
        
        for result in successful_results:
            task_name = result['task_name'].replace(' ', '_')
            
            if 'summary' in result:
                summary_file = os.path.join(output_dir, f"{task_name}_summary.{format}")
                if format == 'csv':
                    result['summary'].to_csv(summary_file)
                elif format == 'xlsx':
                    result['summary'].to_excel(summary_file)
                elif format == 'pickle':
                    result['summary'].to_pickle(summary_file)
            
            if 'results' in result:
                results_file = os.path.join(output_dir, f"{task_name}_results.{format}")
                if format == 'csv':
                    result['results'].to_csv(results_file)
                elif format == 'xlsx':
                    result['results'].to_excel(results_file)
                elif format == 'pickle':
                    result['results'].to_pickle(results_file)
        
        print(f"Results saved to {output_dir} in {format} format")
    
    @classmethod
    def from_task_list(
        cls, 
        task_configs: List[Dict[str, Any]], 
        max_threads: Optional[int] = None
    ) -> 'Gu99Batch':
        """
        Create batch calculator from a list of task configurations.
        
        Args:
            task_configs: List of dictionaries with task parameters
            max_threads: Maximum threads to use
            
        Returns:
            Configured Gu99Batch instance
        """
        batch = cls(max_threads=max_threads)
        
        for i, config in enumerate(task_configs):
            task_name = config.get('task_name', f'Task_{i}')
            batch.add_task(
                config['aln_file'],
                *config.get('tree_files', []),
                cluster_name=config.get('cluster_name'),
                trees=config.get('trees', []),
                check_trees=config.get('check_trees', True),
                task_name=task_name
            )
        
        return batch

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
        if check_tree_file(*tree_files):
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
            index=range(0, np.size(self.calculator._results(), 0) )
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
        
        # Calculate and display effective site numbers
        self.type1_effective_number = self.calculator._results1().shape[0]
        self.type2_effective_number = self.calculator._results2().shape[0]
        print(f"Type1 Effective Number of Sites is {self.type1_effective_number}, Type2 Effective Number of Sites is {self.type2_effective_number}")

    def _calculate(self) -> None:
        """Create a new Effective number of site Calculator and complete the calculation steps."""
        self.calculator = _effectivecpp.create_calculator(self.input, self.cluster_name)
        self.calculator.calculate()

    def type1_results(self) -> pd.DataFrame:
        """Generate Type 1 results.
        
        Retrieves and formats Type 1 functional divergence results from the calculator.
        Type 1 results typically represent sites under positive selection or
        functional constraints in specific lineages.
        
        Returns:
            DataFrame containing Type 1 analysis results with proper column names
            and numbered index.
        """
        columns_names = get_colnames(self.calculator._r_names())
        results = pd.DataFrame(self.calculator._results1(), columns=columns_names)
        results.index.name = "Number"
        return results

    def type2_results(self) -> pd.DataFrame:
        """Generate Type 2 results.
        
        Retrieves and formats Type 2 functional divergence results from the calculator.
        Type 2 results typically represent sites with altered evolutionary rates
        between different lineages or functional classes.
        
        Returns:
            DataFrame containing Type 2 analysis results with proper column names
            and numbered index.
        """
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
    @property
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
            index=[i for i in self.calculator._kept()]
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
    @property
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
            index=[i for i in self.calculator._kept()]
        )
        results.index.name = "Position"
        return results

    def _help(self) -> None:
        """Print help information for _typeOneAnalysiscpp."""
        print(help(_typeOneAnalysiscpp))

# End of file