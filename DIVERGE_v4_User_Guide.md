# DIVERGE v4 Complete User Guide

## Table of Contents
1. [Introduction](#introduction)
2. [Installation](#installation)  
3. [Input File Requirements](#input-file-requirements)
4. [Core Analysis Methods](#core-analysis-methods)
5. [SuperCluster Algorithm](#supercluster-algorithm)
6. [Batch Processing](#batch-processing)
7. [Advanced Features](#advanced-features)
8. [Performance Optimization](#performance-optimization)
9. [Examples and Use Cases](#examples-and-use-cases)
10. [Troubleshooting](#troubleshooting)

## Introduction

DIVERGE v4 is a comprehensive bioinformatics package for analyzing functional divergence in protein family evolution. It implements multiple statistical methods to detect sites under positive selection and functional constraints across different evolutionary lineages.

### Key Features
- **Multiple Analysis Methods**: Gu99, Gu2001, Type2, Asym, Effective, FDR, RVS, TypeOneAnalysis
- **SuperCluster Algorithm**: Novel approach for large multi-gene family analysis
- **Parallel Processing**: High-performance batch processing with C++ backend
- **Conservation Weighting**: Advanced scoring with local sequence conservation
- **Flexible Input**: Support for multiple alignment and tree formats

## Installation

### From PyPI (Recommended)
```bash
pip install diverge4
```

### From Source
```bash
git clone https://github.com/zjupgx/diverge4.git
cd diverge4
pip install pybind11
python setup.py install
```

### Dependencies
- Python 3.9+
- NumPy, SciPy, Pandas
- Biopython, scikit-bio
- PyBind11 (for C++ extensions)
- Matplotlib, Plotly (visualization)

## Input File Requirements

### Alignment Files
- **Formats**: FASTA or CLUSTAL
- **Type**: Amino acid sequences only
- **Quality**: Remove gaps in terminal regions
- **Naming**: Sequence names must match exactly between alignment and trees

```python
# Example alignment format (FASTA)
>Sequence1
MKLLVLGLFI--PSLVLG
>Sequence2
MKLLVLSIFI--PSLVLG
```

### Phylogenetic Trees
- **Format**: Newick format
- **Requirements**: 
  - Minimum depth of 3 levels
  - **CRITICAL**: Remove all internal node names
  - Branch lengths optional but recommended
- **Naming**: Terminal node names must match alignment sequence names

```python
# Example tree format (Newick)
((Sequence1:0.1,Sequence2:0.1):0.05,Sequence3:0.15);
```

## Core Analysis Methods

### 1. Gu99 Analysis (Type I Functional Divergence)

Tests for altered evolutionary rates after gene duplication events.

```python
from diverge import Gu99

# Basic usage
gu99 = Gu99("alignment.aln", "cluster1.tree", "cluster2.tree")
print(gu99.summary)
print(gu99.results.head())

# With custom cluster names
gu99 = Gu99("alignment.aln", "cluster1.tree", "cluster2.tree", 
            cluster_name=["EGFR", "ERBB2"])

# Access results
theta_values = gu99.summary  # Theta coefficients and statistics
site_scores = gu99.results   # Site-specific posterior probabilities

# Functional distance calculation (3+ clusters)
if gu99.fundist() is not None:
    distances = gu99.fundist()
    print("Functional distances between clusters:")
    print(distances)
```

### 2. Gu2001 Analysis (Updated Type I Method)

Improved method for Type I functional divergence detection.

```python
from diverge import Gu2001

gu2001 = Gu2001("alignment.aln", "cluster1.tree", "cluster2.tree")
print("Theta coefficient:", gu2001.summary.iloc[0, 0])
print("Sites with Qk > 0.9:", sum(gu2001.results.iloc[:, 0] > 0.9))
```

### 3. Type2 Analysis (Type II Functional Divergence)

Detects changes in amino acid properties between clusters.

```python
from diverge import Type2

type2 = Type2("alignment.aln", "cluster1.tree", "cluster2.tree")
print("Alpha coefficient:", type2.summary.iloc[0, 0])
```

### 4. Asymmetric Analysis

Tests for asymmetric evolutionary patterns using outgroup information.

```python
from diverge import Asym

asym = Asym("alignment.aln", "ingroup1.tree", "ingroup2.tree", "outgroup.tree")
results = asym.results()  # Results by outgroup cluster
```

### 5. Effective Sites Calculation

Determines effective number of sites for Type I and Type II analyses.

```python
from diverge import Effective

effective = Effective("alignment.aln", "cluster1.tree", "cluster2.tree")
print(f"Type I effective sites: {effective.type1_effective_number}")
print(f"Type II effective sites: {effective.type2_effective_number}")

# Detailed results
type1_sites = effective.type1_results()
type2_sites = effective.type2_results()
```

### 6. False Discovery Rate (FDR)

Applies multiple testing correction to control false discovery rate.

```python
from diverge import Fdr

fdr = Fdr("alignment.aln", "cluster1.tree", "cluster2.tree")
type1_corrected = fdr.type1_results()  # FDR-corrected Type I results
type2_corrected = fdr.type2_results()  # FDR-corrected Type II results
```

### 7. Rate Variation Among Sites (RVS)

Estimates site-specific evolutionary rates.

```python
from diverge import Rvs

rvs = Rvs("alignment.aln", "cluster1.tree", "cluster2.tree")
print("Rate parameters:", rvs.summary)
rates = rvs.results()  # Site-specific rate estimates
```

### 8. TypeOneAnalysis (Multi-cluster Extension)

Analyzes complex patterns across multiple clusters (requires exactly 3 clusters).

```python
from diverge import TypeOneAnalysis

analysis = TypeOneAnalysis("alignment.aln", "cluster1.tree", 
                          "cluster2.tree", "cluster3.tree")
print("State probabilities:", analysis.summary)
patterns = analysis.results()  # Site-specific pattern probabilities
```

## SuperCluster Algorithm

The SuperCluster algorithm is designed for large-scale analysis of multi-gene families, reducing computational complexity from exponential to linear.

### Basic Usage

```python
from diverge import SuperCluster

# Simple analysis with all possible cluster groupings
super_cluster = SuperCluster("alignment.aln", "cluster1.tree", "cluster2.tree", 
                           "cluster3.tree", "cluster4.tree")

# Access results
print("Summary statistics:")
print(super_cluster.summary)

# Site-specific scores for each cluster comparison
results = super_cluster.results
print(f"Analysis performed on {len(results.columns)} positions")
print(f"Number of cluster groupings: {len(results.index)}")
```

### Advanced SuperCluster Options

```python
# Mode selection
super_cluster_simple = SuperCluster("alignment.aln", *tree_files, mode="Simple")
# Simple mode: only one-vs-rest comparisons

super_cluster_all = SuperCluster("alignment.aln", *tree_files, mode="ALL")  
# ALL mode: all possible pairwise comparisons (default)

# Different distance metrics
super_cluster_jc = SuperCluster("alignment.aln", *tree_files, metric="jc69")
# Options: 'hamming' (default), 'jc69', 'k2p', etc.

# Parallel processing control
super_cluster = SuperCluster("alignment.aln", *tree_files, 
                           parallel=True, max_threads=8)

# Conservation weighting
conswins_params = {
    'cons_win_len': 3,      # Window size for conservation calculation
    'lambda_param': 0.7     # Weight: 70% original + 30% conservation  
}
super_cluster_weighted = SuperCluster("alignment.aln", *tree_files, 
                                    conswins=conswins_params)
```

### Working with SuperCluster Results

```python
# Interpret results
results = super_cluster.results
summary = super_cluster.summary

# Find highly divergent sites (Qk > 0.9)  
for group_idx, row in results.iterrows():
    high_sites = row[row > 0.9]
    if len(high_sites) > 0:
        print(f"Group {group_idx}: {len(high_sites)} sites with Qk > 0.9")
        print(f"Positions: {high_sites.index.tolist()}")

# Group information
group_list = super_cluster.group_list
tree_dict = super_cluster.tree_dict
for i, (group1, group2) in enumerate(group_list):
    print(f"Comparison {i}: Clusters {group1} vs {group2}")
```

## Batch Processing

### Gu99Batch for High-Performance Analysis

The Gu99Batch class enables parallel processing of multiple Gu99 analyses, bypassing Python's GIL limitations.

```python
from diverge import Gu99Batch

# Create batch calculator
batch = Gu99Batch(max_threads=8)  # Use 8 CPU cores

# Add multiple tasks
task_configs = [
    {
        'aln_file': 'dataset1.aln',
        'tree_files': ['d1_cluster1.tree', 'd1_cluster2.tree'],
        'task_name': 'Dataset_1'
    },
    {
        'aln_file': 'dataset2.aln', 
        'tree_files': ['d2_cluster1.tree', 'd2_cluster2.tree'],
        'task_name': 'Dataset_2'
    }
]

# Add tasks to batch
for config in task_configs:
    batch.add_task(
        config['aln_file'],
        *config['tree_files'],
        task_name=config['task_name']
    )

# Execute batch processing
batch.calculate_batch()

# Retrieve results
results = batch.get_results()
successful_results = batch.get_successful_results()

# Print summary
batch.print_summary()

# Save results to files
batch.save_results('output_directory', format='csv')
```

### Performance Benchmarking

```python
# Benchmark parallel vs sequential performance
benchmark_results = super_cluster.benchmark_parallel_vs_sequential(runs=3)

print(f"Speedup: {benchmark_results['speedup']:.2f}x")
print(f"Sequential time: {benchmark_results['sequential_time']:.2f}s")
print(f"Parallel time: {benchmark_results['parallel_time']:.2f}s")
```

## Advanced Features

### Conservation Window Scoring

Incorporate local sequence conservation into functional divergence scores:

```python
from diverge.super_cluster import conservation_window_score

# Manual conservation scoring
scores = [0.8, 0.6, 0.9, 0.4, 0.7]  # Original SDP scores
alignment_seqs = ["MKLLV", "MKLIV", "MKLLV"]  # Alignment sequences
window_len = 2
lambda_param = 0.7

weighted_scores = conservation_window_score(scores, alignment_seqs, 
                                          window_len, lambda_param)
```

### Custom Tree Construction

Build phylogenetic trees from distance matrices:

```python
from diverge.super_cluster import tree_construct, sk_aln_read

# Read alignment and construct tree
aln = sk_aln_read("alignment.aln")
distance_matrix, tree = tree_construct(aln)

# Use constructed tree in analysis
gu99 = Gu99("alignment.aln", trees=[tree])
```

### Result Filtering

Apply quality filters to remove unreliable predictions:

```python
# Enable filtering in analysis
gu99_filtered = Gu99("alignment.aln", "cluster1.tree", "cluster2.tree", 
                     filter=True)

super_cluster_filtered = SuperCluster("alignment.aln", *tree_files, 
                                    filter=True)
```

## Performance Optimization

### Memory Management

```python
# Clear cached results to free memory
analysis.clear_cache()

# For large datasets, process in chunks
def process_large_dataset(alignment_chunks, trees):
    results = []
    for chunk in alignment_chunks:
        chunk_result = SuperCluster(chunk, *trees)
        results.append(chunk_result.results)
        chunk_result.clear_cache()  # Free memory
    return results
```

### Parallel Processing Guidelines

```python
# Optimal thread configuration
import os
optimal_threads = min(os.cpu_count(), len(cluster_groups))

# For SuperCluster analysis
super_cluster = SuperCluster("alignment.aln", *tree_files,
                           parallel=True, 
                           max_threads=optimal_threads)

# For batch processing  
batch = Gu99Batch(max_threads=optimal_threads)
```

## Examples and Use Cases

### Complete Analysis Pipeline

```python
from diverge import SuperCluster, CalPipe

# Method 1: Using CalPipe for automated analysis
pipeline = CalPipe("alignment.aln", "tree1.tree", "tree2.tree", "tree3.tree")
comprehensive_results = pipeline.run_all_analyses()

# Method 2: Manual SuperCluster analysis with conservation
conswins = {'cons_win_len': 3, 'lambda_param': 0.8}
super_cluster = SuperCluster("ERBB_family.fas", "EGFR.tree", "ERBB2.tree", 
                           "ERBB3.tree", "ERBB4.tree", 
                           parallel=True, conswins=conswins)

# Analyze results
results = super_cluster.results
summary = super_cluster.summary

# Find significant sites
significant_sites = {}
for group_idx, row in results.iterrows():
    sites_09 = row[row > 0.9]
    sites_067 = row[(row > 0.67) & (row <= 0.9)]
    
    significant_sites[f"Group_{group_idx}"] = {
        'high_confidence': sites_09.index.tolist(),
        'medium_confidence': sites_067.index.tolist()
    }

print("Significant functional divergence sites:")
for group, sites in significant_sites.items():
    print(f"{group}: {len(sites['high_confidence'])} high, "
          f"{len(sites['medium_confidence'])} medium confidence sites")
```

### Large-Scale Multi-Dataset Analysis

```python
import os
from diverge import Gu99Batch

# Prepare multiple datasets
datasets = []
for dataset_dir in os.listdir("datasets/"):
    aln_file = os.path.join(dataset_dir, "alignment.aln")
    tree_files = [os.path.join(dataset_dir, f) for f in os.listdir(dataset_dir) 
                  if f.endswith('.tree')]
    
    datasets.append({
        'aln_file': aln_file,
        'tree_files': tree_files,
        'task_name': dataset_dir
    })

# Create and run batch analysis
batch = Gu99Batch.from_task_list(datasets, max_threads=16)
batch.calculate_batch()

# Process results
successful_results = batch.get_successful_results()
for result in successful_results:
    task_name = result['task_name']
    summary = result['summary']
    theta = summary.iloc[0, 0]  # Theta coefficient
    
    print(f"{task_name}: θ = {theta:.4f}")
    
    # Find highly divergent sites
    results_df = result['results']
    high_sites = results_df[results_df.iloc[:, 0] > 0.9]
    print(f"  High confidence sites: {len(high_sites)}")
```

## Troubleshooting

### Common Issues and Solutions

#### 1. Tree Format Problems
**Error**: "Tree depth is less than 3"
```python
# Solution: Check and fix tree files
from diverge.binding import check_tree_file, read_tree

try:
    check_tree_file("problematic.tree")
except ValueError as e:
    print(f"Tree validation failed: {e}")
    
# Manually inspect tree
tree = read_tree("problematic.tree")
print(f"Tree depth: {max(len(tree.trace(tree.root, clade)) for clade in tree.get_terminals())}")
```

#### 2. Sequence Name Mismatches
**Error**: "Terminal X not found in alignment"
```python
# Solution: Check sequence name consistency
from Bio import AlignIO, Phylo

aln = AlignIO.read("alignment.aln", "fasta")
tree = Phylo.read("tree.tree", "newick")

aln_names = [record.id for record in aln]
tree_names = [terminal.name for terminal in tree.get_terminals()]

print("Alignment sequences:", aln_names)
print("Tree terminals:", tree_names)
print("Mismatched names:", set(tree_names) - set(aln_names))
```

#### 3. Memory Issues with Large Datasets
```python
# Solution: Use chunked processing or reduce thread count
super_cluster = SuperCluster("large_alignment.aln", *tree_files,
                           max_threads=4,  # Reduce threads
                           parallel=True)

# Or process smaller subsets
def process_alignment_subset(aln_file, start_pos, end_pos, tree_files):
    # Implementation for processing alignment subsets
    pass
```

#### 4. Performance Issues
```python
# Solution: Optimize parameters and use batch processing
# Check available cores
import os
print(f"Available CPU cores: {os.cpu_count()}")

# Use appropriate thread count
batch = Gu99Batch(max_threads=min(8, os.cpu_count()))

# Monitor performance
import time
start_time = time.time()
super_cluster = SuperCluster("alignment.aln", *tree_files, parallel=True)
execution_time = time.time() - start_time
print(f"Analysis completed in {execution_time:.2f} seconds")
```

### Validation and Quality Control

```python
# Validate input files before analysis
def validate_inputs(aln_file, tree_files):
    from diverge.super_cluster import pre_check, aln_read
    from diverge.binding import read_tree
    
    try:
        # Check alignment
        aln = aln_read(aln_file)
        print(f"Alignment: {len(aln)} sequences, {len(aln[0])} positions")
        
        # Check trees
        trees = [read_tree(tf) for tf in tree_files]
        if pre_check(aln_file, trees):
            print("All input files validated successfully")
            return True
    except Exception as e:
        print(f"Validation failed: {e}")
        return False

# Use before analysis
if validate_inputs("alignment.aln", tree_files):
    super_cluster = SuperCluster("alignment.aln", *tree_files)
```

---

## Summary

DIVERGE v4 provides a comprehensive suite of tools for functional divergence analysis, from simple pairwise comparisons to complex multi-cluster analyses. The SuperCluster algorithm and parallel processing capabilities make it suitable for large-scale evolutionary studies, while the variety of analysis methods ensures appropriate statistical approaches for different biological questions.

Key recommendations:
- Start with SuperCluster for multi-gene family analysis
- Use batch processing for multiple datasets
- Apply conservation weighting for improved accuracy  
- Validate input files before analysis
- Monitor performance and adjust parameters as needed

For additional support and examples, refer to the Tutorial directory and case studies included with the package.