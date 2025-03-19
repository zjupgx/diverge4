
# DIVERGE API 调用文档

## 安装与设置

### 安装

```bash
pip install diverge
```

### 导入

```python
# 导入主要分析类
from diverge import Gu99, Gu2001, Type2, Asym, Effective, Fdr, Rvs, TypeOneAnalysis

# 导入实用函数
from diverge import read_tree, check_tree, load_tree_file
```

## 基本工作流

DIVERGE 的典型工作流包括：
1. 准备序列比对文件和系统树文件
2. 实例化相应的分析类
3. 检索结果摘要
4. 获取详细结果
5. 可选：执行额外分析如功能距离计算
6. 可选：导出或可视化结果

## 详细 API 调用

### 1. Gu99（Type-I 功能分化分析）

```python
# 基本用法，使用文件路径
analysis = Gu99(
    aln_file="path/to/alignment.aln",     # 序列比对文件
    "path/to/tree1.tree",                 # 第一个树文件
    "path/to/tree2.tree",                 # 第二个树文件
    cluster_name=["CladeA", "CladeB"]     # 自定义集群名称(可选)
)

# 使用 Tree 对象
from Bio import Phylo
tree1 = Phylo.read("path/to/tree1.tree", "newick")
tree2 = Phylo.read("path/to/tree2.tree", "newick")

analysis = Gu99(
    aln_file="path/to/alignment.aln",
    trees=[tree1, tree2],
    cluster_name=["CladeA", "CladeB"]
)

# 获取结果
summary = analysis.summary()      # 获取参数摘要
results = analysis.results()      # 获取位点详细结果
fundist = analysis.fundist()      # 计算功能距离（仅当有≥3个集群时）

# 打印摘要
print(summary)
# 示例输出:
#           CladeA-CladeB
# Theta       0.123456
# SE          0.056789
# MLE         0.134567
# ...

# 查看详细结果
print(results.head())
# 示例输出:
#    Position   AlphaC   BetaC   Score   P-value
# 1        1    0.123    0.456   0.789    0.032
# 2        2    0.234    0.567   0.890    0.045
# ...

# 绘制功能距离图形（仅当fundist不为None时）
if fundist is not None:
    import matplotlib.pyplot as plt
    fig = analysis.plot_distance()
    plt.show()
```

### 2. Gu2001（改进的 Type-I 分析）

```python
# 基本用法
analysis = Gu2001(
    aln_file="path/to/alignment.aln",
    "path/to/tree1.tree",
    "path/to/tree2.tree",
    cluster_name=["CladeA", "CladeB"]
)

# 获取结果
summary = analysis.summary()
results = analysis.results()

# 打印帮助信息
analysis._help()  # 打印C++模块相关信息
```

### 3. Type2（Type-II 功能分化分析）

```python
# 基本调用
type2_analysis = Type2(
    aln_file="path/to/alignment.aln", 
    "path/to/tree1.tree", 
    "path/to/tree2.tree",
    cluster_name=["Mammals", "Birds"]
)

# 获取结果
summary = type2_analysis.summary()
results = type2_analysis.results()

print(summary)
# 示例输出:
#           Mammals-Birds
# Theta       0.345678
# SE          0.078912
# ...

# 筛选显著结果
significant_sites = results[results['P-value'] < 0.05]
print(f"发现 {len(significant_sites)} 个显著的 Type-II 功能分化位点")
```

### 4. Asym（非对称分化分析）

```python
# 基本调用
asym_analysis = Asym(
    aln_file="path/to/alignment.aln",
    "path/to/tree1.tree", 
    "path/to/tree2.tree",
    "path/to/tree3.tree"
)

# 获取结果 - 注意Asym只有results方法，没有summary
results = asym_analysis.results()
print(results)

# 示例输出:
#   Cluster Number of Outgroup   Clade1-Clade2   Clade1-Clade3   Clade2-Clade3
# 1                          1          0.123          0.234          0.345
# 2                          2          0.456          0.567          0.678
# 3                          3          0.789          0.890          0.901
```

### 5. Effective（有效位点分析）

```python
# 基本调用
eff_analysis = Effective(
    aln_file="path/to/alignment.aln",
    "path/to/tree1.tree", 
    "path/to/tree2.tree"
)

# 实例化时会自动打印有效位点数量
# 示例输出: "Type1 Effective Number of Sites is 15, Type2 Effective Number of Sites is 8"

# 获取详细结果
type1_results = eff_analysis.type1_results()
type2_results = eff_analysis.type2_results()

print(type1_results.head())
# 示例输出：
#    Number   Position   Score   ...
# 0       1         45   0.987   ...
# 1       2         67   0.954   ...
```

### 6. Fdr（假发现率控制）

```python
# 基本调用
fdr_analysis = Fdr(
    aln_file="path/to/alignment.aln",
    "path/to/tree1.tree", 
    "path/to/tree2.tree",
    "path/to/tree3.tree"
)

# 获取 Type-I 和 Type-II FDR 控制结果
type1_fdr = fdr_analysis.type1_results()
type2_fdr = fdr_analysis.type2_results()

print(type1_fdr.head())
# 示例输出:
#   alpha   pi0    FDR     Position   q-value
# 1  0.01   0.85   0.025        134     0.015
# ...

# 根据 q-value 筛选位点
significant_sites = type1_fdr[type1_fdr['q-value'] < 0.05]
```

### 7. Rvs（进化速率变异分析）

```python
# 基本调用（注意：通常只需要一个树）
rvs_analysis = Rvs(
    aln_file="path/to/alignment.aln",
    "path/to/tree.tree"
)

# 获取结果
summary = rvs_analysis.summary()
results = rvs_analysis.results()

print(summary)
# 示例输出:
#                 Value
# Mean             0.781
# Variance         0.324
# ...

print(results.head())
# 示例输出:
#    Position   Xk    Rk
# 1         1  2.34  0.78
# 2         2  1.56  0.45
# ...
```

### 8. TypeOneAnalysis（多集群 Type-I 分析）

```python
# 基本调用（适用于3个或更多集群）
type1_analysis = TypeOneAnalysis(
    aln_file="path/to/alignment.aln",
    "path/to/tree1.tree", 
    "path/to/tree2.tree",
    "path/to/tree3.tree",
    cluster_name=["Vertebrates", "Insects", "Plants"]
)

# 获取结果
summary = type1_analysis.summary()
results = type1_analysis.results()

print(summary)
# 示例输出包含S0-S4状态的相关统计

print(results.head())
# 详细位点分析结果
```

## 树文件处理工具

```python
# 读取树文件
tree = read_tree("path/to/tree.nwk")  # 自动检测newick或nexus格式

# 检查单个树的有效性
is_valid = check_tree(tree)

# 加载树文件并检查有效性
tree_str = load_tree_file("path/to/tree.nwk", check=True)
```

## 完整工作流示例

### 示例1：Type-I 与 Type-II 功能分化分析并导出结果

```python
import pandas as pd
import matplotlib.pyplot as plt
from diverge import Gu99, Type2

# 1. 执行 Type-I 分析
type1 = Gu99(
    "protein_family.aln",
    "clade1.tree", 
    "clade2.tree",
    cluster_name=["Group1", "Group2"]
)

# 2. 执行 Type-II 分析
type2 = Type2(
    "protein_family.aln",
    "clade1.tree", 
    "clade2.tree",
    cluster_name=["Group1", "Group2"]
)

# 3. 获取结果
type1_summary = type1.summary()
type1_results = type1.results()
type2_summary = type2.summary()
type2_results = type2.results()

# 4. 筛选显著位点
alpha = 0.05  # 显著性水平
type1_significant = type1_results[type1_results['P-value'] < alpha]
type2_significant = type2_results[type2_results['P-value'] < alpha]

# 5. 找出特异性和共享位点
type1_positions = set(type1_significant.index)
type2_positions = set(type2_significant.index)
common_positions = type1_positions.intersection(type2_positions)
type1_unique = type1_positions - type2_positions
type2_unique = type2_positions - type1_positions

print(f"Type-I 特异性位点: {len(type1_unique)}")
print(f"Type-II 特异性位点: {len(type2_unique)}")
print(f"共享显著位点: {len(common_positions)}")

# 6. 导出结果
with pd.ExcelWriter('functional_divergence_results.xlsx') as writer:
    type1_summary.to_excel(writer, sheet_name='Type1_Summary')
    type1_results.to_excel(writer, sheet_name='Type1_Results')
    type2_summary.to_excel(writer, sheet_name='Type2_Summary')
    type2_results.to_excel(writer, sheet_name='Type2_Results')
    
    # 创建显著位点摘要
    pd.DataFrame({
        'Type1_Specific': list(type1_unique),
        'Type2_Specific': list(type2_unique),
        'Common': list(common_positions)
    }).to_excel(writer, sheet_name='Significant_Sites')

# 7. 可视化
fig, ax = plt.subplots(figsize=(10, 6))
x = range(1, len(type1_results) + 1)
ax.scatter(x, type1_results['Score'], label='Type-I', alpha=0.6, s=20)
ax.scatter(x, type2_results['Score'], label='Type-II', alpha=0.6, s=20)
ax.set_xlabel('序列位置')
ax.set_ylabel('功能分化得分')
ax.set_title('Type-I vs Type-II 功能分化得分比较')
ax.legend()
plt.savefig('functional_divergence_comparison.png', dpi=300)
plt.show()
```

### 示例2：计算有效位点数量并应用 FDR 校正

```python
from diverge import Effective, Fdr

# 1. 计算有效位点数量
eff = Effective(
    "protein_family.aln",
    "clade1.tree", 
    "clade2.tree"
)

# 2. 应用 FDR 校正
fdr = Fdr(
    "protein_family.aln",
    "clade1.tree", 
    "clade2.tree"
)

# 3. 获取结果
type1_effective = eff.type1_results()
type2_effective = eff.type2_results()
type1_fdr = fdr.type1_results()
type2_fdr = fdr.type2_results()

# 4. 根据 q-value 筛选可靠位点
q_threshold = 0.05
reliable_type1 = type1_fdr[type1_fdr['q-value'] < q_threshold]
reliable_type2 = type2_fdr[type2_fdr['q-value'] < q_threshold]

print(f"效FDR校正后可靠的Type-I位点: {len(reliable_type1)}")
print(f"FDR校正后可靠的Type-II位点: {len(reliable_type2)}")

# 5. 位点重叠分析
type1_positions = set(reliable_type1.index)
type2_positions = set(reliable_type2.index)
common = type1_positions.intersection(type2_positions)

print(f"同时显示Type-I和Type-II功能分化的位点数量: {len(common)}")
if common:
    print(f"重叠位点: {common}")
```

## 常见错误和解决方法

| 错误消息 | 可能原因 | 解决方法 |
|----------|----------|----------|
| `Tree file 'xxx' is not in newick or nexus format.` | 树文件格式不正确 | 确保树文件是有效的Newick或Nexus格式 |
| `Tree depth is less than 3, please check your tree` | 树的深度不足 | 需要使用深度大于等于3的有效系统发育树 |
| `FileNotFoundError: xxx does not exist` | 指定的文件不存在 | 检查文件路径是否正确 |
| `TypeError: cluster_name must be a list` | 集群名称参数类型错误 | 确保集群名称是列表类型 |
| `Need to create at least two clusters first.` | 缺少足够的树文件 | 添加更多的树文件作为参数 |
| `Gu2001 method failed. The input trees do not contain branch length!` | 树缺少分支长度 | 使用包含分支长度的树文件 |

## 最佳实践

1. **数据准备**
   - 确保序列对齐质量高，无大量空位
   - 树文件应包含分支长度
   - 树的拓扑结构应足够深（至少3层）

2. **分析流程**
   - 通常先进行 Type-I 和 Type-II 分析比较
   - 对显著位点应用 FDR 校正减少假阳性
   - 使用 Effective 分析评估实际的功能分化位点数量

3. **结果解读**
   - Type-I：涉及进化速率变化的功能分化
   - Type-II：涉及物理化学性质改变的功能分化
   - P值小于0.05的位点通常被认为显著
   - theta值反映功能分化的程度

4. **性能考虑**
   - 大型数据集可能需要较长计算时间
   - 对于超大序列比对，考虑先提取感兴趣区域
   - 确保有足够内存处理大型序列比对

## 参考信息

- DIVERGE 版本: 4.0.0+
- 输入格式:
  - 序列比对: Clustal (.aln)
  - 树文件: Newick (.nwk, .tree) 或 Nexus (.nex)
- 依赖: Python 3.6+, Biopython, NumPy, Pandas, Matplotlib
