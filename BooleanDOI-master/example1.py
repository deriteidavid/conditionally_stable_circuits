import networkx as nx
import BooleanDOI_processing as BDOIp
import BooleanDOI_TargetControl as BDOItc

#input files
BooleanRule_filename = 'example1.txt'
network_name = 'example1'

#semi-input files
expanded_filename = network_name + '_expanded_edges.txt'

#Set parameters
#Note the node name for Gread is the index (integer), one can encode the nodename by adding prefix and suffix
#If the node name from the input file is not this simple, one need to create a dictionary to record the mapping
prefix,suffix='n','n'

#read the original graph
f = open(BooleanRule_filename,'r')
lines = f.readlines()
f.close()
Gread,readnodes = BDOIp.form_network(lines, sorted_nodename = False)

#form expanded network
G_expanded = BDOIp.Get_expanded_network(Gread, prefix=prefix, suffix=suffix)

#parameters for target control
max_itr = 500
custom_score_index = '3'            #can be chosen from '0'~'4'
#forbidden_nodes, avail_nodes, forbidden_node_states = [], [], []  #one can add constrain to the candidate node states by setting these three parameters
Target = set(['n2n'])
#Target = set(['n2n','n4n'])

#initialization
TDOI_BFS, flagged, potential = {}, {}, {} #Domain of Influence

#first step, search through DOI of single node
BDOItc.update_single_DOI(G_expanded, network_name, TDOI_BFS, flagged, potential)

#Target control algorithm
solutions, solution_found = BDOItc.GRASP_target_control(G_expanded, Target, max_itr, TDOI_BFS, flagged, potential, custom_score_index=custom_score_index)
print solutions
