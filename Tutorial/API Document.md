# DIVERGE 包 API 文档

## 介绍

DIVERGE（**D**etection of **I**nter-specific sequence con**VERGE**nce）是一个用于检测蛋白质家族中功能分化的 Python 包。本包提供了多种分析方法来识别在进化过程中可能与功能分化相关的氨基酸位点。

## 安装

```bash
pip install diverge
```

## 主要功能

DIVERGE 包含以下主要分析类型：

- **Gu99**: 基于 Gu (1999) 方法的功能分化分析
- **Gu2001**: 基于 Gu (2001) 方法的功能分化分析
- **Type2**: Type-II 功能分化分析
- **Asym**: 非对称分化分析
- **Effective**: 有效位点数量分析
- **Fdr**: 错误发现率控制分析
- **Rvs**: 进化速率变异分析
- **TypeOneAnalysis**: Type-I 功能分化精细分析

## 快速开始

```python
from diverge import Gu99, Type2, Effective

# Type-I 功能分化分析
type1 = Gu99("alignment.aln", "cluster1.tree", "cluster2.tree")
print(type1.summary())
print(type1.results())

# Type-II 功能分化分析
type2 = Type2("alignment.aln", "cluster1.tree", "cluster2.tree")
print(type2.summary())
print(type2.results())

# 有效位点分析
eff = Effective("alignment.aln", "cluster1.tree", "cluster2.tree")
print(eff.type1_results())
print(eff.type2_results())
```

## 核心 API

### BaseAnalysis

所有分析类的基类，提供通用功能。

```python
BaseAnalysis(aln_file, *tree_files, cluster_name=None, trees=[], 
            check_trees=True, calculator_module=None)
```

**参数：**
- `aln_file`: 字符串，序列比对文件路径
- `*tree_files`: 字符串，树文件路径
- `cluster_name`: 列表或None，集群名称
- `trees`: 列表，Biopython.Phylo Tree 对象
- `check_trees`: 布尔值，是否检查树的有效性
- `calculator_module`: C++ 计算模块

**方法：**
- `summary()`: 返回分析结果的摘要 DataFrame
- `results()`: 返回详细结果 DataFrame
- `clear_cache()`: 清除缓存的结果
- `_help()`: 打印计算模块的帮助信息

### Gu99

Type-I 功能分化分析（Gu 1999 方法）。

```python
Gu99(aln_file, *tree_files, cluster_name=None, trees=[])
```

**特有方法：**
- `fundist()`: 计算功能距离，返回 DataFrame
- `plot_distance(figsize=(10, 6), title="Functional Distance Between Clusters")`: 绘制功能距离图

### Gu2001

Type-I 功能分化分析（Gu 2001 方法），考虑了树长度影响。

```python
Gu2001(aln_file, *tree_files, cluster_name=None, trees=[])
```

### Type2

Type-II 功能分化分析。

```python
Type2(aln_file, *tree_files, cluster_name=None, trees=[])
```

### Asym

非对称性功能分化分析。

```python
Asym(aln_file, *tree_files, cluster_name=None, trees=[])
```

### Effective

有效位点数量分析。

```python
Effective(aln_file, *tree_files, cluster_name=None, trees=[])
```

**特有属性：**
- `type1_effective_number`: Type-I 有效位点数量
- `type2_effective_number`: Type-II 有效位点数量

**特有方法：**
- `type1_results()`: 返回 Type-I 有效位点 DataFrame
- `type2_results()`: 返回 Type-II 有效位点 DataFrame

### Fdr

错误发现率控制分析。

```python
Fdr(aln_file, *tree_files, cluster_name=None, trees=[])
```

**特有方法：**
- `type1_results()`: 返回 Type-I FDR 控制结果
- `type2_results()`: 返回 Type-II FDR 控制结果

### Rvs

进化速率变异分析。

```python
Rvs(aln_file, *tree_files, cluster_name=None, trees=[])
```

### TypeOneAnalysis

Type-I 功能分化精细分析（用于3个或以上集群）。

```python
TypeOneAnalysis(aln_file, *tree_files, cluster_name=None, trees=[])
```

## 实用函数

### 树处理函数

```python
# 读取树文件
tree = read_tree("path/to/tree.nwk")

# 检查树的有效性
is_valid = check_tree(tree)

# 检查多个树文件的有效性
are_valid = check_tree_file("tree1.nwk", "tree2.nwk")

# 加载树文件，可选检查有效性
tree_str = load_tree_file("tree.nwk", check=True)
```

### 便捷分析函数

```python
# Type-I 分析快捷方式
result = analyze_type1("alignment.aln", "cluster1.tree", "cluster2.tree")

# Type-II 分析快捷方式
result = analyze_type2("alignment.aln", "cluster1.tree", "cluster2.tree")

# 有效位点分析快捷方式
result = analyze_effective("alignment.aln", "cluster1.tree", "cluster2.tree")
```

## 输出处理

DIVERGE 提供了将结果保存为不同格式的函数。

```python
from diverge.io.writers import save_results_to_csv, save_results_to_excel

# 保存为 CSV
save_results_to_csv(analysis.results(), "results.csv")

# 保存为 Excel (多个表格)
save_results_to_excel({
    "Summary": analysis.summary(),
    "Results": analysis.results()
}, "analysis_results.xlsx")
```

## 并行处理

对于大型数据集，DIVERGE 提供了并行处理功能。

```python
from diverge.core.parallel import parallel_map, chunked_parallel

# 并行处理简单任务
results = parallel_map(my_function, items, n_jobs=4)

# 并行处理分块数据
results = chunked_parallel(process_batch, large_dataset, chunk_size=100)
```

## 实例说明

### 分析 Type-I 功能分化并导出结果

```python
from diverge import Gu99
from diverge.io.writers import save_results_to_excel

# 执行分析
analysis = Gu99("protein.aln", "clade1.tree", "clade2.tree", "clade3.tree",
               cluster_name=["Vertebrates", "Plants", "Fungi"])

# 查看结果摘要
summary = analysis.summary()
print(summary)

# 获取详细结果
results = analysis.results()
print(results)

# 计算功能距离并绘图
analysis.plot_distance(title="Functional Distance Between Clades")
plt.savefig("functional_distance.png")

# 导出所有结果
save_results_to_excel({
    "Summary": summary,
    "Results": results,
    "FunctionalDistance": analysis.fundist()
}, "type1_analysis.xlsx")
```

### 使用 Tree 对象直接分析

```python
from Bio import Phylo
from diverge import Type2

# 加载树对象
tree1 = Phylo.read("tree1.nwk", "newick")
tree2 = Phylo.read("tree2.nwk", "newick")

# 直接使用树对象进行分析
analysis = Type2.from_trees("alignment.aln", trees=[tree1, tree2],
                           cluster_name=["Group1", "Group2"])

# 或使用工厂方法
analysis = Type2("alignment.aln", trees=[tree1, tree2])
```

## 注意事项

1. 树深度需要至少为3，否则会抛出 ValueError。
2. 分析结果通常包括：
   - 统计参数（theta 等）
   - 位点特异性分数
   - 显著性检验结果
3. Gu99 和 Type2 分析通常用于鉴定可能与功能分化相关的位点。
4. 使用 Effective 分析可以评估功能分化位点的有效数量。
5. Fdr 分析有助于控制多重检验中的假阳性率。

## 错误处理

常见错误类型和解决方法：

- `FileNotFoundError`: 确保文件路径正确。
- `ValueError: Tree depth is less than 3`: 确保树的深度足够。
- `TypeError: cluster_name must be a list`: 集群名称必须是列表类型。

## 版本兼容性

- 此 API 文档适用于 DIVERGE 4.0.0 及以上版本。
- 依赖 Python 3.6+ 和 Biopython 1.76+。

## 参考文献

1. Gu, X. (1999). Statistical methods for testing functional divergence after gene duplication. Molecular Biology and Evolution, 16(12), 1664-1674.
2. Gu, X. (2001). Maximum-likelihood approach for gene family evolution under functional divergence. Molecular Biology and Evolution, 18(4), 453-464.
3. Gu, X. (2006). A simple statistical method for estimating type-II (cluster-specific) functional divergence of protein sequences. Molecular Biology and Evolution, 23(10), 1937-1945.
