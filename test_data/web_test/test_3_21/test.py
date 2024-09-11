from diverge.utils import CalPipe
pipeline = CalPipe("./CASP.aln","./inner_clade_14.tree","./inner_clade_2.tree",'./inner_clade_20.tree')

from diverge import _gu99cpp
gu99 = _gu99cpp.create_calculator(["./CASP.aln","./inner_clade_14.tree","./inner_clade_2.tree",'./inner_clade_20.tree'])
gu99.calculate()
gu99._r_names()
from diverge import _type2cpp
type2 = _type2cpp.create_calculator(["./CASP.aln","./inner_clade_14.tree","./inner_clade_2.tree",'./inner_clade_20.tree'])
type2.calculate()
type2._r_names()

from diverge import _rvscpp
rvs = _rvscpp.create_calculator(["./CASP.aln","./inner_clade_14.tree","./inner_clade_2.tree",'./inner_clade_20.tree'])
rvs.calculate()
rvs._r_names()

from diverge import _fdrcpp
fdr = _fdrcpp.create_calculator(["./CASP.aln","./inner_clade_14.tree","./inner_clade_2.tree",'./inner_clade_20.tree'])
fdr.calculate()
fdr._r_names()

from diverge.binding import Gu99
gu99 = Gu99("./CASP.aln","./inner_clade_14.tree","./inner_clade_2.tree",'./inner_clade_20.tree')
gu99.results()

from diverge.binding import Gu2001
gu2001 = Gu2001("./CASP.aln","./inner_clade_14.tree","./inner_clade_2.tree",'./inner_clade_20.tree')
gu2001.results()

from diverge.binding import Type2
type2 = Type2("./CASP.aln","./inner_clade_14.tree","./inner_clade_2.tree",'./inner_clade_20.tree')
type2.results()