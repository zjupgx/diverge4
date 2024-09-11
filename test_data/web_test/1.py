from Bio import Phylo
tree1 = Phylo.read("./inner_clade_12.tree",'newick')
Phylo.draw(tree1)
tree2 = Phylo.read("../treelength/cl1.tree",'newick')
Phylo.draw(tree2)


tree3 = Phylo.read("../CASP.tree",'newick')

from diverge i