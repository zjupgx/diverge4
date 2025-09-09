import math
from collections import Counter
from Bio.AlignIO import MultipleSeqAlignment
from Bio import Phylo
import numpy as np
import pandas as pd

def calculate_entropy_conservation(aas):
    valid_aas = [i for i in aas if i!="-"]
    counter = Counter(valid_aas)
    total = len(valid_aas)
    entropy = 0
    for count in counter.values():
        p = count / total
        if p > 0:
            entropy -= p * math.log2(p)
    max_entropy = math.log2(20)
    conservation = 1-(entropy/max_entropy)
    entropy_conservation_score = max_entropy - entropy
    return entropy_conservation_score,conservation

def aln2dict(aln):
  dict = {}
  for record in aln:
    dict[record.id] = record
  return dict

def get_group_alns(aln,trees):
  group_keys = [i for i in range(len(trees))]
  groups = {}
  for group_key in group_keys:
    groups[group_key] = [i.name for i in  trees[group_key].get_terminals()]
  aln_dict = aln2dict(aln)
  group_alns = {}
  for group_key in group_keys:
    group_records = [aln_dict[i] for i in groups[group_key]]
    group_alns[group_key] = MultipleSeqAlignment(group_records)
  return group_alns


def get_group_aas(group_alns,position):
  group_aas = {}
  for group_key in group_alns.keys():
    group_aas[group_key] = group_alns[group_key][:,position]
  return group_aas

def filter_results(results,aln,trees,col_threshold=0.15,group_threshold=0.3,sameaa_threshold=0.8,entropy_threshold=3):
  group_alns = get_group_alns(aln,trees)
  entropys = []
  for position in results.columns:
    group_aas = get_group_aas(group_alns,position)
    col_aas = [i for group_aa in group_aas.values() for i in group_aa]
    valid_col_aas = [i for i in col_aas if i!="-"]
    group_entropys = []
    conservations = []
    for index,group_aa in enumerate(group_aas.values()):
      if Counter(group_aa)['-']/len(group_aa) > group_threshold:
        results[position][index] = 0.1*results[position][index]
      group_entropy,conservation = calculate_entropy_conservation(group_aa)
      conservations.append(conservation)
      group_entropys.append(group_entropy)
    entropys.append(group_entropys)
    if (Counter(col_aas)['-']/len(col_aas) > col_threshold) or (np.max([i/len(valid_col_aas) for i in Counter(valid_col_aas).values()])>sameaa_threshold):
      results[position] = 0.1*results[position]
      continue
    if np.min(group_entropys)>entropy_threshold:
      results[position] = (1.2/(1+np.min(conservations)))*results[position]
      continue
  entropys = np.array(entropys)
  return(results,entropys)

