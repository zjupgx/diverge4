"""
DIVERGE - Detection of Functional Divergence in Protein Families

Main module for DIVERGE analysis tools.
"""

__version__ = '4.0.7'
__author__ = 'DIVERGE Team'

# Export main functionality
from .binding import (
    # Core analysis classes
    Gu99,
    Gu2001,
    Type2,
    Asym,
    Effective,
    Fdr,
    Rvs,
    TypeOneAnalysis,
    
    # Utility functions
    read_tree,
    check_tree,
    check_tree_file,
    load_tree_file,
)

from .utils import (
    CalPipe
)

from .filter import filter_results

from .super_cluster import SuperCluster

# Export visualization functionality
from .treeplot import (
    draw_phylogenetic_tree,
)






