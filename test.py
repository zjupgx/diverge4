# from diverge import _gu99cpp


# calculator = _gu99cpp.create_calculator(["./test_data/web_test/gu99/CASP.aln","./test_data/web_test/inner_clade_11.tree","./test_data/web_test/inner_clade_5.tree","./test_data/web_test/inner_clade_9.tree"])
# calculator.calculate()



# from Bio import Phylo
# tree = Phylo.read("./test_data/web_test/inner_clade_9.tree",'newick')
# Phylo.draw(tree)

# tree2 = Phylo.read("./test_data/cl2.tree",'newick')
# Phylo.draw(tree2)

# tree3 = Phylo.read('./test_data/CASP.tree','newick')
# Phylo.draw(tree3)
# tree_ex = Phylo.read("./test_data/CASP.tree",'newick')


# Phylo.write()


# tree = Phylo.read("./test_data/web_test/NTRK_1.tree",'newick')
# Phylo.draw(tree)

# tree2 = Phylo.read("./test_data/web_test/NTRK_2.tree",'newick')
# Phylo.draw(tree2)


# from diverge.super_cluster import *
# from Bio import Phylo
# trees,_,_ = super_cluster("./test_data/CASP.aln","./test_data/cl1.tree","./test_data/cl2.tree","./test_data/cl3.tree")



import pandas as pd
import matplotlib.pyplot as plt
# from altair_saver import save
# from PIL import Image

# Load the data
species_count = pd.read_csv("./web/statics/species_counts.csv")
species_count.columns = ['species', 'num']


# 创建图表
fig, ax = plt.subplots(figsize=(10, 6))  # 设置图表大小

bars = ax.bar(species_count['species'], species_count['num'], color='#97c4ed')

# 添加文本标签
for bar in bars:
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width() / 2, height, f'{height}', 
            ha='center', va='bottom', fontsize=8)

# 设置轴标签和标题样式
ax.set_xlabel('Species', fontsize=13, fontweight='bold', color='black')
ax.set_ylabel('Sequences number', fontsize=13, fontweight='bold', color='black')
ax.tick_params(axis='x', labelrotation=-45, labelsize=10, labelcolor='black')

# 去除图表边框
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_color('black')
ax.spines['bottom'].set_color('black')

# 保存为300 DPI的PNG文件
plt.tight_layout()
plt.savefig('species_count_chart.svg', dpi=300)

# 显示图表
plt.show()
