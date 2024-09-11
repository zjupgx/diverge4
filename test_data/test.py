
from Bio import Phylo,AlignIO,SeqIO
from diverge import *

tree1 = Phylo.read("./cl1.tree",'newick')
tree2 = Phylo.read("./cl2.tree",'newick')
tree3 = Phylo.read("./cl3.tree",'newick')

# 测试alignment中序列数量对结果的影响
aln_path = "./CASP.aln"
aln = AlignIO.read(aln_path, "clustal")

full_sequence_rvs = Rvs(aln_path,trees=[tree1])

print(Rvs(aln_path,trees=[tree1]).results())
print(Rvs(aln_path,trees=[tree1]).summary())

terminals_records = []
for tree in [tree1,tree2]:
  for clade in tree.find_clades(terminal=True):
    terminals_records.append(clade.name)

seq_records = []
for seq in aln:
  if seq.id in terminals_records:
    seq_records.append(seq)
print("1:",len(seq_records))

SeqIO.write(seq_records, "test.aln", "clustal")


terminals_records = []
for tree in [tree1,tree2,tree3]:
  for clade in tree.find_clades(terminal=True):
    terminals_records.append(clade.name)

seq_records = []
for seq in aln:
  if seq.id in terminals_records:
    seq_records.append(seq)

print("2:",len(seq_records))
SeqIO.write(seq_records, "test2.aln", "clustal")

print(Rvs("test.aln",trees=[tree1]).results())
print(Rvs("test.aln",trees=[tree1]).summary())


print(Rvs("test2.aln",trees=[tree1]).results())
print(Rvs("test2.aln",trees=[tree1]).summary())
# import _gu99cpp
sequence_13_rvs = Rvs("test.aln",trees=[tree1])
sequence_18_rvs = Rvs("test2.aln",trees=[tree1])
sequence_13_aln = AlignIO.read("test.aln", "clustal")
sequence_18_aln = AlignIO.read("test2.aln", "clustal")

# calculator = _gu99cpp.create_calculator(["./CASP.aln",tree1,tree2,tree3],['tree1','tree2','tree3'])
# calculator.calculate()

T = Gu99("../test_data/CASP.aln",trees=[tree1,tree2,tree3],cluster_name=["trees1","trees2","trees3"])
summary = T.summary()
results = T.results()
print(summary)
print("------------")
print(results)
# test
T = Gu99("../test_data/CASP.aln", "../test_data/cl1.tree", "../test_data/cl2.tree","../test_data/cl3.tree",cluster_name=["tree1","tree2","tree3"])
summary = T.summary()
results = T.results()
print(summary)
print("------------")
print(results)
# test gu2001
T = Gu2001("../test_data/CASP.aln", "../test_data/treelength/cl1.tree", "../test_data/treelength/cl2.tree","../test_data/treelength/cl3.tree")
summary = T.summary()
results = T.results()
print(summary)
print("------------")
print(results)
# test type2 
T = Type2("../test_data/CASP.aln", "../test_data/cl1.tree", "../test_data/cl2.tree","../test_data/cl3.tree")
summary = T.summary()
results = T.results()
print(summary)
print("------------")
print(results)
# test asym
T = Asym("../test_data/CASP.aln", "../test_data/cl1.tree", "../test_data/cl2.tree","../test_data/cl3.tree")
results = T.results()
print(results)
# test effective
T = Effective("../test_data/CASP.aln", "../test_data/cl1.tree", "../test_data/cl2.tree")
# T._help()
results1 = T.type1_results()
results2 = T.type2_results()
print(results1)
print("------------")
print(results2)
# test fdr
T = Fdr("../test_data/CASP.aln", "../test_data/cl1.tree", "../test_data/cl2.tree","../test_data/cl3.tree")
results1 = T.type1_results()
results2 = T.type2_results()
print(results1)
print("------------")
print(results2)
# test rvs
T = Rvs("../test_data/CASP.aln", "../test_data/cl1.tree")
summary = T.summary()
results = T.results()
print(summary)
print("------------")
print(results)
# test type one analysis
T = TypeOneAnalysis("../test_data/CASP.aln", "../test_data/cl1.tree", "../test_data/cl2.tree","../test_data/cl3.tree")
summary = T.summary()
results = T.results()
print(summary)
print("------------")
print(results)