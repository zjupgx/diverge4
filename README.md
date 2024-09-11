# DIVERGE v4

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

DIVERGE v4 is a Python package designed for large-scale analysis of functional divergence across multi-gene families. It is an upgrade of the widely used DIVERGE software, incorporating a modular Python structure and a user-friendly web server. This package allows for the identification of amino acid sites undergoing significant evolutionary shifts, helping to uncover functional divergence after gene duplication.

![framework](.\statics\framework.jpg)

## Features

- **Modular Python package**: Built for scalability, with 10 customizable modules for functional divergence analysis.
- **Web interface**: A user-friendly web server developed using the Streamlit framework, making the package accessible even without programming knowledge.
- **Super-Cluster Algorithm**: Designed for large-scale multi-gene families, efficiently predicting functional divergence sites.
- **Comprehensive Database**: A database containing analysis results of human gene families, with search functions based on UniProtKB, Ensembl, HGNC IDs, or gene names.

## Installation

To install the DIVERGE v4 Python package, use pip:

```bash
pip install diverge4
```
If you prefer to compile the package from source using setup.py, you will need to install the pybind11 library, which provides the C++ bindings for Python used in this package. You can install it via pip:
```bash
pip install pybind11
```
Once `pybind11` is installed, you can compile DIVERGE v4 by running the following commands:
```bash
git clone https://github.com/zjupgx/diverge4.git
cd diverge4
python setup.py install
```
`pybind` is necessary because `DIVERGE v4` uses C++ for its core data structures and computationally intensive tasks, which are exposed to Python via `pybind11`.
## Quick Start

### Example Usage

Here’s a simple example to get started with DIVERGE v4. This example demonstrates how to load a multiple sequence alignment (MSA) file and perform a basic analysis:

```python
from diverge import Gu99
# Perform DIVERGE's Gu99 analysis
gu99 = Gu99("./test_data/CASP.aln","./test_data/cl1.tree","./test_data/cl2.tree","./test_data/cl3.tree")
# View the results
result = gu99.results()
print(result)
```
For more detailed examples, check the [tutorials](./Tutorial.ipynb) and [casestudy](./CaseStudy/ERBB_analysis.ipynb).
## Main Functional Modules

DIVERGE v4 provides a variety of independent computing processes to create custom pipelines for functional divergence analysis. Below is a list of the main functions available in the package:

| **Function** | **Description** |
|--------------|-----------------|
| **Type-I Divergence (Gu99 method)** | Detect type-I functional divergence using the Gu (1999) method. |
| **Type-I Divergence (Gu2001 method)** | Detect type-I functional divergence using the Gu (2001) method. Requires a phylogenetic tree file that contains branch length data. |
| **Type-II Divergence** | Detect type-II functional divergence of gene families. |
| **Rate Variation Among Sites (RVS)** | Estimate rate variations among sites for a given cluster, as described in Gu and Zhang (1997). **Only one cluster is allowed per run.** |
| **Functional Distance Analysis** | Estimate type-I functional distance between pairs of clusters and compute type-I functional branch lengths. **Requires at least three homologous gene clusters.** |
| **FDR for Predictions** | Calculate the false discovery rate of functionally diverging sites. |
| **Asymmetric Test for Type-I Functional Divergence** | Statistically test whether the degree of type-I functional divergence differs between two duplicate genes. **Requires three clusters.** |
| **Effective Number of Sites Related to Functional Divergence** | Estimate the effective number of sites related to type-I or type-II functional divergence. **Requires two clusters.** |
| **Gene-Specific Type-I Analysis** | Site-specific posterior profile for predicting gene-specific type-I functional divergence-related sites. **Requires three clusters.** |
| **Super-Cluster** | Perform large-scale functional divergence analysis using the Super-Cluster method, designed for multi-gene families with a large number of subfamilies. |

Each module can be used independently, allowing for a customized analysis pipeline based on your research needs.

## Web Server and Database
If you prefer not to use Python, the web interface provides a simple yet powerful tool for performing functional divergence analysis.

The web server also includes a comprehensive database of human protein families, offering access to the results of precomputed functional divergence analyses. Users can search the database using identifiers such as UniProtKB ID, Ensembl ID, HGNC ID, or gene name.

### Web Server Features:

- Functional Divergence Analysis: Upload your multiple sequence alignment and phylogeny files, and receive results on functional divergence directly from the web server.
- Data Visualization: Explore protein families, visualize amino acid sites undergoing functional divergence, and interact with the results.
- Searchable Database: Search the database of human protein families and access detailed functional annotations, including Gene Ontology terms and pathways.

To access the web server:

1. Visit [DIVERGE Web Server](https://pgx.zju.edu.cn/diverge).
2. Upload your multiple sequence alignment and phylogeny files.
3. Receive and download the analysis results from the results page.

The web server is designed to handle large datasets, enabling researchers to run complex analyses without local computational constraints.

## Citation

If you use DIVERGE in your research, please cite:

```

```
## License

DIVERGE v4 is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.
