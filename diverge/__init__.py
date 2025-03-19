"""
DIVERGE - Detection of Functional Divergence in Protein Families

Main module for DIVERGE analysis tools.
"""

__version__ = '4.0.0'
__author__ = 'DIVERGE Team'

# 导出主要功能
from .binding import (
    # 核心分析类
    Gu99,
    Gu2001,
    Type2,
    Asym,
    Effective,
    Fdr,
    Rvs,
    TypeOneAnalysis,
    
    # 实用函数
    read_tree,
    check_tree,
    check_tree_file,
    load_tree_file,
)

# 导出可视化功能
from .treeplot import (
    plot_tree,
    highlight_nodes,
)

# 提供便捷函数
def analyze_type1(aln_file, *tree_files, **kwargs):
    """Convenience function to run Type 1 analysis."""
    return Gu99(aln_file, *tree_files, **kwargs)

def analyze_type2(aln_file, *tree_files, **kwargs):
    """Convenience function to run Type 2 analysis."""
    return Type2(aln_file, *tree_files, **kwargs)

def analyze_effective(aln_file, *tree_files, **kwargs):
    """Convenience function to calculate effective sites."""
    return Effective(aln_file, *tree_files, **kwargs)

__version__ = '4.0.1'


