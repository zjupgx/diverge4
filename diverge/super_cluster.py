from Bio import Phylo
import pandas as pd
import numpy as np
import os
import shutil
import tempfile
import dendropy
import itertools
import random
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceMatrix
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from .utils import create_temp_folder,write_temp_tree
from .binding import Gu2001

calculator = DistanceCalculator('identity')
constructor = DistanceTreeConstructor()
temp_folder = create_temp_folder()
def aln_read(aln_file,aln_type):
  aln = AlignIO.read(aln_file, aln_type)
  return aln

def tree_construct(aln):
  dm = calculator.get_distance(aln)
  tree = constructor.nj(dm)
  return dm,tree

def sub_dm(dm,c_list):
  clades_name = pd.Index(dm.names)[[c_list]].to_list()
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
  return cluster_list

def tree_reconstruct(dm,cluster):
  cluster1_list = np.isin(cluster,1)
  cluster2_list = np.isin(cluster,2)
  sub_dm1 = sub_dm(dm,cluster1_list)
  sub_dm2 = sub_dm(dm,cluster2_list)
  tree1 = constructor.nj(sub_dm1)
  tree2 = constructor.nj(sub_dm2)
  return tree1,tree2

def get_cluster(aln_file,*tree_files):
  aln = AlignIO.read(aln_file,'clustal')
  dm,tree = tree_construct(aln)
  tree_cluster = [0]*len(dm.names)
  i = 1
  for tree_file in tree_files:
    tree_terminal = [i.name for i in Phylo.read(tree_file,'newick').get_terminals()]
    list = [dm.names.index(j) for j in tree_terminal]
    for k in list: 
      tree_cluster[k] = i
    i +=1
  return tree_cluster,i-1

def re_clean_tree(tree):
  for clade in tree.get_nonterminals():
    clade.name = None
  return tree

def get_super_cluster_list(aln_file,*tree_files):
  tree_cluster,cluster_num = get_cluster(aln_file,*tree_files)
  cluster_list = sep_cluster(tree_cluster,cluster_num)
  aln = aln_read(aln_file,'clustal')
  dm,_ = tree_construct(aln)
  super_cluster_list = []
  for cluster in cluster_list:
    tree1,tree2 = tree_reconstruct(dm,cluster)
    tree1 = re_clean_tree(tree1)
    tree2 = re_clean_tree(tree2)
    super_cluster_list.extend([[tree1,tree2]])
  return(super_cluster_list)

def write_tree_file(super_cluster_list):
  tree_path_list = []
  for i in range(len(super_cluster_list)):
    tree1,tree2 = super_cluster_list[i][0],super_cluster_list[i][1]
    tree1_path = write_temp_tree(tree1,temp_folder,'cluster1',str(i))
    tree2_path = write_temp_tree(tree2,temp_folder,'cluster2',str(i))
    tree_path_list.extend([[tree1_path,tree2_path]])
  
  return tree_path_list

# def get_super_cluster_pp2(aln_file,*tree_files):
#   """
#   Make a composite Q-site plotting: each site has M number of Q values, denoted by Q1,…, QM.
#   Let Qmax=max(Q1,…, QM), which is used to represent in the composite Q-site plotting.
#   """
#   super_cluster_list = get_super_cluster_list(aln_file,*tree_files)
#   tree_path_list = write_tree_file(super_cluster_list)
#   results_list = []
#   for tree_path in tree_path_list:
#     tree1_path = tree_path[0]
#     tree2_path = tree_path[1]
#     gu2001 = Gu2001(aln_file,tree1_path,tree2_path)
#     results = gu2001.results().iloc[:,0]
#     results_list.extend([[results]])
#   results_array = np.reshape(np.array(results_list),(3, -1))
#   return results_array,super_cluster_list

import numpy as np
from concurrent.futures import ThreadPoolExecutor

def get_super_cluster_pp(aln_file, *tree_files):
    """
    Make a composite Q-site plotting: each site has M number of Q values, denoted by Q1,…, QM.
    Let Qmax=max(Q1,…, QM), which is used to represent in the composite Q-site plotting.
    """
    super_cluster_list = get_super_cluster_list(aln_file, *tree_files)
    tree_path = write_tree_file(super_cluster_list)
    results_list = []

    def process_tree_path(tree1_path, tree2_path):
        gu2001 = Gu2001(aln_file, tree1_path, tree2_path)
        results = gu2001.results().iloc[:, 0]
        return results

    with ThreadPoolExecutor() as executor:
        for result in executor.map(process_tree_path, [t[0] for t in tree_path], [t[1] for t in tree_path]):
            results_list.extend([[result]])

    results_array = np.reshape(np.array(results_list), (3, -1))
    return results_array, super_cluster_list


def super_cluster(aln_file,*tree_files):
  pp_list,tree_list = get_super_cluster_pp(aln_file,*tree_files)
  max_group = np.argmax(np.bincount(np.argmax(pp_list,axis=0)))
  trees = tree_list[max_group]
  return trees

if __name__ == "main":
  from diverge.super_cluster import *
  from Bio import Phylo
  trees = super_cluster("./test_data/CASP.aln","./test_data/cl1.tree","./test_data/cl2.tree","./test_data/cl3.tree")
  # result2,_ = get_super_cluster_pp("./test_data/CASP.aln","./test_data/cl1.tree","./test_data/cl2.tree","./test_data/cl3.tree")

  
  