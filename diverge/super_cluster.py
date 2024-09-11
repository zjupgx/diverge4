from typing import Any
from Bio import Phylo,AlignIO
import pandas as pd
import numpy as np
import os
import itertools
from tqdm import tqdm
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
# from .utils import create_temp_folder,write_temp_tree,remove_temp_folder
# from .binding import Gu99,Type2
from diverge import Gu99,Type2
from diverge.utils import printv,tqdm_joblib
from joblib import Parallel, delayed



# temp_folder = create_temp_folder()

def aln_read(aln_file):
  try:
    aln = AlignIO.read(aln_file, 'clustal')
  except ValueError:
    try:
      aln = AlignIO.read(aln_file, 'fasta')
    except ValueError:
      print("==============\n",aln_file,"==============\n")
      raise Exception("The alignment file is not in fasta or clustal format")
  return aln

def tree_construct(aln,dist_calc='identity'):
  calculator = DistanceCalculator(dist_calc)
  constructor = DistanceTreeConstructor()
  dm = calculator.get_distance(aln)
  tree = constructor.nj(dm)
  return dm,tree

def sub_dm(dm,c_list):
  clades_name = pd.Index(dm.names)[np.where(c_list)].to_list()
  sub_dm = DistanceMatrix(clades_name)
  for clade_name1,clade_name2 in itertools.combinations(clades_name,2):
    sub_dm[clade_name1,clade_name2] = dm[clade_name1,clade_name2]
  return sub_dm

def get_group_list(group_num):
  n = group_num
  nums = list(range(1, n + 1))
  num_range = range(1, n//2+1)
  group_list = []
  for num in num_range:
    combos = itertools.combinations(nums, num)
    for combo in combos:
        group1 = list(combo)
        group2 = list(set(nums) - set(combo))
        if([group2,group1] in group_list):
          continue
        else:
          group_list.extend([[group1, group2]])
  return group_list

def sep_cluster(tree_cluster,cluster_num):
  group_list = get_group_list(cluster_num)
  cluster_list = []
  for group in group_list:
    cluster = np.zeros(len(tree_cluster))
    group1,gourp2 = group[0],group[1]
    cluster[np.isin(tree_cluster,group1)] = 1
    cluster[np.isin(tree_cluster,gourp2)] = 2
    cluster_list.extend([cluster])
  return cluster_list,group_list

def tree_reconstruct(dm,cluster):
  constructor = DistanceTreeConstructor()
  cluster1_list = np.isin(cluster,1)
  cluster2_list = np.isin(cluster,2)
  sub_dm1 = sub_dm(dm,cluster1_list)
  sub_dm2 = sub_dm(dm,cluster2_list)
  tree1 = constructor.nj(sub_dm1)
  tree2 = constructor.nj(sub_dm2)
  return tree1,tree2

def get_cluster(aln,*tree_files,trees:list=[]):
  dm,tree = tree_construct(aln)
  tree_cluster = [0]*len(dm.names)
  i = 1
  tree_dict = {}
  if trees!=[]:
    for tree in trees:
      tree_id = f"subtree_{i}"
      tree_dict[i] = tree_id
      tree_terminal = [i.name for i in tree.get_terminals()]
      t_list = [dm.names.index(j) for j in tree_terminal]
      for k in t_list: 
        tree_cluster[k] = i
      i +=1
  else:
    for tree_file in tree_files:
      tree = Phylo.read(tree_file,'newick')
      tree_id = os.path.basename(tree_file).split('.')[0]
      tree_dict[i] = tree_id
      tree_terminal = [i.name for i in tree.get_terminals()]
      t_list = [dm.names.index(j) for j in tree_terminal]
      for k in t_list: 
        tree_cluster[k] = i
      i +=1
  return tree_cluster,i-1,tree_dict

def re_clean_tree(tree):
  for clade in tree.get_nonterminals():
    clade.name = None
  # return tree.root_at_midpoint()
  return tree

def get_super_cluster_list(aln,*tree_files,trees:list=[]):
  tree_cluster,cluster_num,tree_dict = get_cluster(aln,*tree_files,trees=trees)
  cluster_list,group_list = sep_cluster(tree_cluster,cluster_num)
  dm,_ = tree_construct(aln)
  super_cluster_list = []
  for cluster in cluster_list:
    tree1,tree2 = tree_reconstruct(dm,cluster)
    tree1 = re_clean_tree(tree1)
    tree2 = re_clean_tree(tree2)
    super_cluster_list.extend([[tree1,tree2]])
  return super_cluster_list,group_list,tree_dict

def process_tree(aln_file, super_cluster, sp_type):
    if sp_type == 1:
        calculator = Gu99(aln_file, trees=super_cluster)
    else:
        calculator = Type2(aln_file, trees=super_cluster)
    summary = calculator.summary()
    position = calculator.results().index.values.tolist()
    results = calculator.results().values.tolist()
    return results, position, summary
  
def get_super_cluster_pp(aln_file, *tree_files, sp_type: int = 1, trees: list = [],verbose=True):
    aln = aln_read(aln_file)
    super_cluster_list, group_list, tree_dict = get_super_cluster_list(aln, *tree_files, trees=trees)
    results_list = []
    position_list = []
    summary_list = []
    n_jobs = 8
    if verbose:
      printv(f"Running type{sp_type} functional divergence super cluster analysis...")
      with tqdm_joblib(desc="Processing super cluster groups", total=len(super_cluster_list)) as progress_bar:
        results = Parallel(n_jobs=n_jobs)(delayed(process_tree)(aln_file, super_cluster, sp_type) for super_cluster in super_cluster_list)
      printv('Finish super cluster analysis.')
    else:
      results = Parallel(n_jobs=n_jobs)(delayed(process_tree)(aln_file, super_cluster, sp_type) for super_cluster in super_cluster_list)
    for result in results:
        results_list.append(result[0])
        summary_list.append(result[2])
        if not position_list:
            position_list = result[1]
    
    results_array = np.reshape(np.array(results_list), (len(super_cluster_list), -1))
    
    return results_array, super_cluster_list, group_list, tree_dict, position_list, summary_list

class SuperCluster():
    def __init__(self, aln_file ,*tree_files,trees:list=[],sp_type:int=1,verbose=True):
        """
        Attributes:
          sp_type (int): Using SuperCluster to calculate sp_type I or sp_type II functional divergence.
          aln_file (str): The path to the alignment file.
          tree_files (str): A list of paths to the phylogenetic tree files.
          result (pandas.DataFrame): A dataframe that stores the posterior probability of every site in the alignment for each group in the supercluster.
        """
        self.pp_list, self.tree_list, self.group_list, self.tree_dict , self.position_list,self.summary_list = get_super_cluster_pp(aln_file,*tree_files,sp_type=sp_type,trees=trees,verbose=verbose)
        self.result = pd.DataFrame(self.pp_list,index=[f"{i}" for i in self.group_list],columns=self.position_list)
        self.summary = self.get_summary()
    def get_summary(self):
        col_name = [f"{i}" for i in self.group_list]
        summary_df = pd.DataFrame(index=list(self.summary_list[0].index.values)+["Qk>0.5","Qk>0.67","Qk>0.9"],columns=col_name)
        for i in range(len(self.group_list)):
          # summary_df.iloc[0,i] = f"{round(self.param_list[i][0][0],2)}±{round(self.param_list[i][1][0],2)}"
          # summary_df.iloc[1,i] = round(self.param_list[i][2][0],2)
          index_s = self.summary_list[i].shape[0]
          summary_df.iloc[:index_s,i] = self.summary_list[i].iloc[:,0]
          summary_df.iloc[index_s,i] = sum(self.pp_list[i,:]>0.5)
          summary_df.iloc[index_s+1,i] = sum(self.pp_list[i,:]>0.67)
          summary_df.iloc[index_s+2,i] = sum(self.pp_list[i,:]>0.9)
          
        return summary_df
    
    
if __name__ == "__main__":
  from Bio import AlignIO
  from diverge.utils import auto_split
  from diverge import super_cluster
  test = AlignIO.read("../CaseStudy/4.fas",'fasta')
  subtrees,subtree_clades,exclude_list= auto_split("../CaseStudy/4.fas",plot=False,exclude_level=4,tree_construct_method='upgma',cluster_method='ward')
  SuperCluster_1 = super_cluster.SuperCluster("../CaseStudy/4.fas",sp_type=1,trees=subtrees)
