"""
DIVERGE - Detection of Functional Divergence in Protein Families

Main module for DIVERGE analysis tools.
"""

__version__ = '4.1.0'
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
)

from .utils import (
    CalPipe
)

from .super_cluster import SuperCluster






