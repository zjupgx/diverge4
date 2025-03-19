"""
Tests for base functionality.
"""
import unittest
import tempfile
import os
import numpy as np
from Bio.Phylo.BaseTree import Tree, Clade

from diverge.binding import (
    read_tree,
    check_tree,
    check_tree_depth,
    BaseAnalysis
)

class TestTreeFunctions(unittest.TestCase):
    """Test tree handling functions."""
    
    def test_check_tree_depth(self):
        """Test the tree depth checker."""
        # 创建一个深度为3的树
        tree = Tree(
            root=Clade(
                branch_length=1.0,
                clades=[
                    Clade(branch_length=1.0, name="A", clades=[
                        Clade(branch_length=1.0, name="A1")
                    ]),
                    Clade(branch_length=1.0, name="B")
                ]
            )
        )
        
        # 应该通过检查
        self.assertTrue(check_tree_depth(tree))
        
        # 测试深度不足的树
        shallow_tree = Tree(
            root=Clade(
                branch_length=1.0,
                clades=[
                    Clade(branch_length=1.0, name="A"),
                    Clade(branch_length=1.0, name="B")
                ]
            )
        )
        
        # 应该不通过检查
        self.assertFalse(check_tree_depth(shallow_tree))

# 更多测试... 