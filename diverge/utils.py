from .binding import Gu99,Gu2001,Type2,Fdr,Effective,Asym,TypeOneAnalysis,Rvs
import pandas as pd
import numpy as np
from Bio import Phylo
from Bio.Phylo import BaseTree
from scipy.cluster.hierarchy import linkage, cut_tree,dendrogram
from pymsaviz import MsaViz
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceCalculator
import os
import tempfile
import shutil


class CalPipe():
    """CalPipe class
    wrap some calculate pipeline of diverge
    """
    def __init__(self,aln_file,*tree_files) -> None:
        """
        __init__ 
        :param aln_file: alingnment file
        :param pipe: some calculate pipeline mode, defaults to 'mode1'
          mode1: tree without branch_length,two clusters : gu99,type2,fdr,effective
          mode2: tree without branch_length,at least three clusters : gu99,type2,fdr,functional distance
          mode3: tree with branch_length,two clusters : gu99,gu2001,type2,fdr,effective
          mode4: tree with branch_length,at least three clusters : gu99,gu2001,type2,fdr,functional distance
        :param _result_summary: result summary dict
        """
        self.func_dict = {'Gu99':True,'Gu2001':False,'Type2':True,'Fdr':True,'Effective':False,'FunDist':False,'Type1Analysis':False}
        self.tree_files = tree_files
        self.aln_file = aln_file
        self._pipe_select(*self.tree_files)
        self._pipeline()
        self.result_summary = self._result_summary() 
        self.detail = [key for key in self.func_dict.keys() if self.func_dict[key]==True]
    def _pipe_select(self,*tree_files):
        cluster_num = len(tree_files)
        has_branch_length = _has_branch_length(*tree_files)
        if has_branch_length:
            self.func_dict['Gu2001'] = True
        if cluster_num >=3:
            self.func_dict['FunDist'] = True
        if cluster_num ==2:
            self.func_dict['Effective'] = True
        if cluster_num ==3:
            self.func_dict['Asym'] = True
            self.func_dict['Type1Analysis'] = True
    def _pipeline(self):
        if self.func_dict['Gu99']:
            gu99 = Gu99(self.aln_file,*self.tree_files)
            self.gu99_results = gu99.results()
            self.gu99_summary = gu99.summary()
            if self.func_dict['FunDist']:
                self.fundist_results = gu99.fundist()
        if self.func_dict['Gu2001']:
            gu2001 = Gu2001(self.aln_file,*self.tree_files)
            self.gu2001_results = gu2001.results()
            self.gu2001_summary = gu2001.summary()
        if self.func_dict['Type1Analysis']:
            toa = TypeOneAnalysis(self.aln_file,*self.tree_files)
            self.type1analysis_summary = toa.summary()
            self.type1analysis_results = toa.results()
        if self.func_dict['Type2']:
            type2 = Type2(self.aln_file,*self.tree_files)
            self.type2_summary = type2.summary()
            self.type2_results = type2.results()
        if self.func_dict['Effective']:
            effective = Effective(self.aln_file,*self.tree_files)
            self.type1_effective = effective.type1_results()
            self.type2_effective = effective.type2_results()
        if self.func_dict['Fdr']:
            fdr = Fdr(self.aln_file,*self.tree_files)
            self.type1_fdr = fdr.type2_results()
            self.type2_fdr = fdr.type2_results()
        
        if self.func_dict['Asym']:
            asym = Asym(self.aln_file,*self.tree_files)
            self.asym_results = asym.results()
    def _result_summary(self):
        _result_summary = {}
        if self.func_dict['Gu99']:
            _result_summary['gu99_summary'] = self.gu99_summary
            _result_summary['gu99_results'] = self.gu99_results
        if self.func_dict['Gu2001']:
            _result_summary['gu2001_summary'] = self.gu2001_summary
            _result_summary['gu2001_results'] = self.gu2001_results
        if self.func_dict['Type1Analysis']:
            _result_summary['type1analysis_summary'] = self.type1analysis_summary
            _result_summary['type1analysis_results'] = self.type1analysis_results
        if self.func_dict['Type2']:
            _result_summary['type2_summary'] = self.type2_summary
            _result_summary['type2_results'] = self.type2_results
        if self.func_dict['Effective']:
            _result_summary['type1_effective'] = self.type1_effective
            _result_summary['type2_effective'] = self.type2_effective
        if self.func_dict['Fdr']:
            _result_summary['type1_fdr'] = self.type1_fdr
            _result_summary['type2_fdr'] = self.type2_fdr
        if self.func_dict['FunDist']:
            _result_summary['fundist_results'] = self.fundist_results
        if self.func_dict['Asym']:
            _result_summary['asym_results'] = self.asym_results
        return _result_summary
    def __str__(self):
        strtext = "Diverge calculation pipeline\n"
        for i,j in enumerate(self.detail):
            strtext += "step{}: {}\n".format(i+1,j)
        strtext += "#####################\n"
        strtext += "You can get the result by calling the result_summary attribute or the specific attribute as follow:\n {}".format(self.result_summary.keys())
        return strtext
    def __repr__(self) -> str:
        return self.__str__()
def _has_branch_length(*tree_files):
    k = 0
    # for tree_file in tree_files,check if check_tree(tree_file) is all true or all false
    for tree_file in tree_files:
        k += _check_tree(tree_file)
    try:
        k == len(tree_files) or k == 0
    except:
        ex = Exception('tree file error,please check your tree file')
        raise ex
    if k == len(tree_files):
        return True
    elif k == 0:
        return False

def _check_tree(tree_file):
    tree = Phylo.read(tree_file, 'newick')
    has_branch_length = False
    for clade in tree.find_clades():
        if clade.branch_length:
            has_branch_length = True
            return has_branch_length
    return has_branch_length

# read tree file PHYLIP format
def read_tree(tree_file):
    tree = Phylo.read(tree_file, 'newick')
    ## set brach length to 1 for all clade
    for clade in tree.get_nonterminals():
        clade.branch_length = 1
    return tree

def tree_construct(aln):
    calculator = DistanceCalculator('identity')
    constructor = DistanceTreeConstructor()
    dm = calculator.get_distance(aln)
    tree = constructor.nj(dm)
    for clade in tree.get_nonterminals():
        clade.name = None
    return tree

# from a Bio.Phylo.BaseTree.Tree object, return a distance matrix
# distance is defined as the mini number of nodes between two leaves  
def tree_to_dist(tree):
    clades_name = [clade.name for clade in tree.get_terminals() if clade.name != '' ]
    dist = pd.DataFrame(index=clades_name, columns=clades_name)
    for clade in tree.get_terminals():
        for other in tree.get_terminals():
            if clade.name != other.name:
                dist.loc[clade.name, other.name] = tree.distance(clade, other)
    return dist.fillna(0)



# split tree into clusters
def split_tree(tree_file,n_clusters):
    tree = read_tree(tree_file)
    # construct the distance matrix 
    tree_dist = tree_to_dist(tree)
    ## transform distance matrix into a condensed matrix
    condensed_tree_dist = tree_dist.values[np.triu_indices(len(tree_dist), k=1)]
    Z = linkage(condensed_tree_dist, method='complete')
    clusters = cut_tree(Z,n_clusters=n_clusters).reshape(-1)
    ## from clusters to split tree object and return a list of split tree
    sub_trees = []
    root_distance = []
    for i in range(n_clusters):
        root = tree.common_ancestor(tree_dist.index[np.where(clusters == i)].values)
        # function = lambda x:x-1 if x>i else x
        # index = function(np.where(clusters == i)[0])
        root_distance.append(max([root.distance(i) for i in root.get_terminals()]))
        sub_trees.append(root)

    return sub_trees,root_distance

def write_subtree(tree_file,n_clusters,path):
    sub_trees,_ = split_tree(tree_file,n_clusters)
    for i in range(1,n_clusters+1):
        tree = BaseTree.Tree.from_clade(sub_trees[i-1])
        Phylo.write(sub_trees[i],os.path.join(path,f'subtree{i}'),'newick')

def create_temp_folder():
    temp_folder = tempfile.TemporaryDirectory()
    return temp_folder

def write_temp_tree(tree,temp_folder,tree_name='temp_tree',subdir=None):
    if subdir is not None:
        temp_tree_path = os.path.join(temp_folder.name,subdir,tree_name+'.tree')
        if not os.path.exists(os.path.dirname(temp_tree_path)):
            os.mkdir(os.path.join(temp_folder.name,subdir))
    else:
        temp_tree_path = os.path.join(temp_folder.name,tree_name+'.tree')
    Phylo.write(tree,temp_tree_path,'newick')
    return temp_tree_path

def plot_msa(aln,aln_format='clustal',color_scheme="None",marker_list=None,wrap_length=60,show_grid=True,show_plot=True,save_plot=False,path="./msa.png",marker='v',marker_color="red"):
    mv=MsaViz(aln,color_scheme=color_scheme,wrap_length=wrap_length,show_grid=show_grid, format =aln_format)
    if marker_list is not None:
        mv.add_markers(marker_list,marker=marker,color=marker_color)
    if show_plot:
        mv.plotfig()
    if save_plot:
        mv.savefig(path,pad_inches=0.5)

def view_cutoff_msa(aln,results,colname,color_scheme="None",cutoff=0.5,show_grid=True,show_plot=True,save_plot=False,path="./msa.png",marker='v',marker_color="red",aln_format='clustal'):
    """view msa with cutoff value
    color_scheme : str | None, optional
            Color scheme. If None, `Zappo`(AA) or `Nucleotide`(NT) is set.
            [`Clustal`|`Zappo`|`Taylor`|`Flower`|`Blossom`|`Sunset`|`Ocean`|
            `Hydrophobicity`|`HelixPropensity`|`StrandPropensity`|`TurnPropensity`|
            `BuriedIndex`|`Nucleotide`|`Purine/Pyrimidine`|`None`]
    marker : str | Marker type of matplotlib.
            See https://matplotlib.org/stable/api/markers_api.html for details.
    """
    try:
        results.loc[:,colname]
    except(KeyError):
        raise Exception(f'column {colname} not found in results')
    marker_list = results[results.loc[:,colname] > cutoff].index.tolist()
    plot_msa(aln,aln_format=aln_format,marker_list=marker_list,show_grid=show_grid,show_plot=show_plot,save_plot=save_plot,path=path,marker=marker,color_scheme=color_scheme,marker_color=marker_color)

def pre_check_tree(*tree_files):
    for tree_file in tree_files:
        tree = Phylo.read(tree_file, 'newick')
        for clade in tree.find_clades():
            clade.branch_length = 1
        # get tree depth
        tree_depth = max(tree.depths().values())
        if tree_depth < 3:
            raise Exception(f'tree depth is less than 3, please check your tree file{tree_file}')
            return False    
    return True
