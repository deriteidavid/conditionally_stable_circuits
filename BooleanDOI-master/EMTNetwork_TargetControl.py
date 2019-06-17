import os
import itertools
import re
import random
import ast
import sys
from collections import defaultdict
import networkx as nx
import BooleanDOI_processing as BDOIp
import BooleanDOI_DOI as BDOI
import BooleanDOI_TargetControl as BDOITC


#input files
BooleanRule_filename='EMTNetwork.txt'
result_filename='EMTreduced'

#semi-input files
expanded_filename=result_filename+'_expanded_edges.txt'
node_mapping_filename=result_filename+'_node_mapping.txt'

#output files
reduced_network_filename=result_filename+'_Reduced_Rules.txt'
expanded_graphname=result_filename+'_expanded.gml'
DOI_complete_filename=result_filename+'_DOI_Complete.txt'


removelist=[]
for item in removelist:
  if os.path.isfile(item):
    os.remove(item)


max_itr=1000
repeat_times=1

#read the original graph
prefix,suffix='n',''
#write function to read the Boolean rules too
f = open(BooleanRule_filename,'r')
lines = f.readlines()
f.close()
Gread,readnodes = BDOIp.form_network(lines, sorted_nodename=False)
#readnodes should be the same each time we read the graph

#obtain reduced network
G_reduced=BDOIp.Get_reduced_network(Gread)
BDOIp.write_Boolean_rules(G_reduced,filename=reduced_network_filename)

#read or calculate the expanded network
if os.path.isfile(expanded_filename):
  g = open(expanded_filename,'r')
  lines = g.readlines()
  g.close()
  G_expanded=nx.DiGraph()
  for line in lines:
    parent,child=line.split(' ')
    G_expanded.add_edge(parent.strip(),child.strip())
else:
  G_expanded=BDOIp.Get_expanded_network(G_reduced,prefix=prefix,suffix=suffix)
  g = open(expanded_filename,'a')
  for edge in G_expanded.edges():
    g.write(str(edge[0])+' '+str(edge[1])+'\n')
  g.close()
#BDOIp.draw_gml(G_expanded,expanded_graphname)


#calculate the mapping from string nodename to index
mapping={}              #nodename to number index
inverse_mapping={}      #number index to nodename
for i,node in enumerate(readnodes):
  index=prefix+str(i)+suffix
  mapping[node]=index
  inverse_mapping[index]=node
  mapping['~'+node]='~'+index
  inverse_mapping['~'+index]='~'+node

#write the node mapping if the node mapping file does not exist
if not os.path.isfile(node_mapping_filename):
  f = open(node_mapping_filename,'a')
  for node in readnodes:
    f.write(mapping[node]+' '+node+'\n')
  f.close()

#add composite node to node mapping
for node in G_expanded.nodes():
  if node not in mapping:
    components=node.split('_')
    composite_node='_'.join([inverse_mapping[x] for x in components])
    mapping[composite_node]=node
    inverse_mapping[node]=composite_node

#define the problem
try:
  target_index=sys.argv[1]
except:
  target_index='1'

if target_index=='0':
  Target=set([mapping['EMT']])
elif target_index=='1':
  Target=set([mapping['~EMT']])
elif target_index=='2':
  Target=set([mapping[x] for x in ['~EMT','~SNAI1']])
elif target_index=='3':
  Target=set([mapping[x] for x in ['~EMT','~MEK']])
elif target_index=='4':
  Target=set([mapping[x] for x in ['~EMT','~MEK','~SNAI1']])
else:
  raise ValueError('Wrong target_index!')

try:
  forbidden_index=sys.argv[2]
except:
  forbidden_index='0'

if forbidden_index=='0':
  forbidden_nodes,avail_nodes,forbidden_node_states=[],[],[]
elif forbidden_index=='1':
  forbidden_nodes=[mapping[x] for x in ['Ecadherin','TGFB','TGFBR']]
  avail_nodes,forbidden_node_states=[],[]
elif forbidden_index=='2':
  forbidden_nodes=[mapping[x] for x in ['Ecadherin','TGFB','TGFBR','SNAI2','SNAI1','ZEB1','ZEB2','FOXC2','TWIST1','HEY1','miR200']]
  avail_nodes,forbidden_node_states=[],[]
elif forbidden_index=='3':
  forbidden_nodes=[mapping[x] for x in ['Ecadherin','TGFB','TGFBR']]
  avail_nodes=[]
  forbidden_node_states=[mapping[x] for x in readnodes if x!= 'Bcatenin_memb']
elif forbidden_index=='4':
  forbidden_nodes=[mapping[x] for x in ['Ecadherin']]
  avail_nodes,forbidden_node_states=[],[]
else:
  raise ValueError('Wrong forbidden_index!')

try:
  custom_score_index=sys.argv[3]
except:
  custom_score_index='0'


#initialization
TDOI_BFS, flagged, potential = {}, {}, {} #Domain of Influence

#first step, search through DOI of single node
BDOITC.update_single_DOI(G_expanded, result_filename, TDOI_BFS, flagged, potential)

#Target control algorithm
solutions, solution_found=BDOITC.GRASP_target_control(G_expanded, Target, max_itr, TDOI_BFS, flagged, potential,
                            custom_score_index=custom_score_index,forbidden_nodes=forbidden_nodes, avail_nodes=avail_nodes, forbidden_node_states=forbidden_node_states)

print [inverse_mapping[x] for x in Target]
#print solutions
solutions_original_name = sorted([([inverse_mapping[x] for x in solution[0]],solution[1],solutions[solution]) for solution in solutions.keys()],key = lambda x: (x[1], len(x[0]), -x[2]))
print solutions_original_name

#write down the DOI in the orignial notation
TDOI_BFS_formatted={}
normal_nodes_set=set([x for x in G_expanded.nodes() if '~' not in x and '_' not in x])
for node in readnodes:
  if mapping[node] in normal_nodes_set:
    TDOI_BFS_formatted[node]=[inverse_mapping[x] for x in TDOI_BFS[tuple([mapping[node]])]]
    TDOI_BFS_formatted[BDOITC.Negation_in_expanded(node)]=[inverse_mapping[x] for x in TDOI_BFS[tuple([mapping[BDOITC.Negation_in_expanded(node)]])]]
if not os.path.isfile(DOI_complete_filename):
  f=open(DOI_complete_filename,'a')
  f.write('Domain of influence truncated by BFS:\n')
  for node in readnodes:
    if mapping[node] in normal_nodes_set:
      f.write(str(node)+' : '+str(TDOI_BFS_formatted[node])+'\n')
      f.write(str(BDOITC.Negation_in_expanded(node))+' : '+str(TDOI_BFS_formatted[BDOITC.Negation_in_expanded(node)])+'\n')
  f.close()
