'''
This file contains functions using GRASP algorithm to solve the target control problem in Boolean network model.
The algorithm is described in a paper titled "Target Control in Logical Models Using the Domain of Influence of Nodes"
to be publisehd in the Frontiers of Physiology.
The code is written by Gang Yang, Department of Physics, Penn State University.
'''
import os
import random
import ast
from collections import defaultdict
import networkx as nx
import BooleanDOI_DOI as BDOI



def Negation(node):
  if '.' not in node:
    return str(-int(node))
  else:
    return str(-float(node))

def Negation_in_expanded(node):
  if '~' in node:
    return node[1:]
  else:
    return '~'+node


def source_node(G):
  '''
  return source node of a graph G
  '''
  return [node for node in G.nodes() if nx.ancestors(G,node) <= set([node])]


def source_part(G):
  '''
  return the source part of Graph G, that is all nodes will be source nodes after iteratively removing source nodes
  '''
  Gprime=G.copy()
  source_part=[]
  while True:
    temp = source_node(Gprime)
    if temp!=[]:
      source_part.extend(temp)
      Gprime.remove_nodes_from(temp)
    else:
      break
  return source_part

def sink_node(G):
  '''
  return sink nodes of a graph G
  '''
  return [key for key,value in G.out_degree().iteritems() if value==0]

def sink_part(G):
  '''
  return the sink part of Graph G, that is all nodes will be sink nodes after iteratively removing source nodes
  '''
  Gprime=G.copy()
  Gprime.remove_nodes_from(source_part(G))   #to avoid repetition with source_part, kick out the source nodes first
  sink_part=[]
  while True:
    temp = sink_node(Gprime)
    if temp!=[]:
      sink_part.extend(temp)
      Gprime.remove_nodes_from(temp)
    else:
      break
  return sink_part

def custom_score1(node, measure):
  '''
  Heuristic score of a node as size of the measure, e.g. the logic Domain of Influence
  Node is a single node in the expanded graph
  The key of measure is in the form of tuple, valid measure are like TDOI_BFS, potential
  '''
  return len(measure[tuple([node])])



def custom_score2(node, measure1, measure2, ratio=1):
  '''
  heuristic score of a node as the addition of two different measure, e.g. the logic Domain of Influence
  node is a single node in the expanded graph
  The key of TDOI is in the form of tuple
  '''
  return len(measure1[tuple([node])])+ratio*len(measure2[tuple([node])])

def custom_score3(node, TDOI, Target):
  '''
  heuristic score of a node as penalized size of domain of influence.
  The penalization is done through multiply -1 for those nodes whose LDOI contains the negation of Target.
  node is a single node in the expanded graph
  The key of TDOI is in the form of tuple
  '''
  contains_negation=any([ Negation_in_expanded(x) in TDOI[tuple([node])] for x in Target])
  return (2*int(not(contains_negation))-1)*len(TDOI[tuple([node])])

def custom_score4(node, TDOI, Target, maximum_TDOI):
  '''
  heuristic score of a node as penalized size of domain of influence.
  The penalization is done through subtracting the maximum LDOI from the score for those nodes whose LDOI contains the negation of Target.
  node is a single node in the expanded graph
  The key of TDOI is in the form of tuple
  '''
  contains_negation=any([ Negation_in_expanded(x) in TDOI[tuple([node])] for x in Target])
  return len(TDOI[tuple([node])])-int(contains_negation)*maximum_TDOI


def generate_random_target(expanded_graph, composite_node_delimiter = '_', target_number = None):
  '''
  Function to generate random target for a given expanded graph.
  Notice the target set will be chosen from sink part if the total number of sink nodes of the expanded network is larger than 4 times the target set.
  Otherwise the target set will be chosen from the non-source part.
  Notice the source node is not an appropriate candidate for our control framework as there is no way to influence these source nodes.

  Parameters
  ----------
  expanded_graph: the input expanded_graph as a DiGraph object
  target_number: the size of the desired target set, if target_number is None, it will be a random number between 1 and 20% of the network size

  Returns
  -------
  Target: a set of node states
  '''
  #calculate non-composite node
  normal_nodes = [x for x in expanded_graph.nodes() if composite_node_delimiter not in x]
  #check whether the size of the target size is specified
  if not target_number > 0:
    target_number = random.randint(1,int(0.1*len(normal_nodes)))   #equivalent 20% of the normal nodes
  sinks = sink_part(expanded_graph)
  normal_sinks = [x for x in sinks if x in normal_nodes]
  #calculate candidates for the target set
  if len(normal_sinks) >= 4*target_number:
    candidates = normal_sinks[:]
  else:
    sources = set(source_part(expanded_graph))
    candidates = [x for x in normal_nodes if x not in sources]
  #make sure we have enough candidates nodes, otherwise we may generate a set that does not meet the requirement of target_number
  assert len(candidates) >= 2*target_number
  #iteratively add node to the target set
  target = set()
  for i in range(target_number):
    node = random.choice(candidates)
    target.add(node)
    candidates.remove(node)
    candidates.remove(Negation_in_expanded(node))
  assert len(target) == target_number
  return target

def update_DOI(G_expanded, nodelist, TDOI, flagged, potential):
  '''
  A helper function to calcualte the LDOI of a node set.
  If the nodelist is not stored in the measure before, do the search and store it
  '''
  if tuple(nodelist) not in TDOI:
    TDOI[tuple(nodelist)],complist,conflict,potential[tuple(nodelist)] = BDOI.truncated_node_of_influence_BFS(G_expanded,nodelist)
    flagged[tuple(nodelist)] = not(conflict==[])
    return

def solution_construction(G_expanded, Target, candidates, candidates_score, TDOI_BFS, flagged, potential):
  '''
  The construction phase of the GRASP algorithm to solve the target control problem.
  The pseudocode is available in the Algorithm Table 2 of the paper titled "Target Control in Logical Models Using the Domain of Influence of Nodes"
  to be publisehd in the Frontiers of Physiology.

  Parameters
  ----------
  G_expanded: the input expanded_graph as a DiGraph object
  Target : the target set
  candidates : the available candidate nodes for target control
  candidates_score : the chosen heuristic greedy function
  TDOI_BFS : a dictionary to store the LDOI of each node state set
  flagged : a dictionary to store whether the LDOI of each node state set contains the negation of any target node
  potential : a dictionary to store the composite nodes visited but not included during the search process for LDOI for each node state set

  Returns
  -------
  solution : a list of node states as a control solution
  fail_in_construction : a Boolean value to indicate whether we failed in the construction of a solution
  '''
  #Initialization
  alpha = random.random()
  solution = []
  fail_in_construction = False
  while True:
    #calculate greedy function for each candidate node and set the pass_score
    current_candidates_score = [candidates_score[candidate] for candidate in candidates]
    max_score = max(current_candidates_score)
    min_score = min(current_candidates_score)
    pass_score = min_score+alpha*(max_score-min_score)
    #Here we assume that restricted_candidates obtained this way will not be empty
    restricted_candidates = [x for x in candidates if candidates_score[x]>=pass_score]
    assert len(restricted_candidates) > 0
    selected_candidate = random.choice(restricted_candidates)
    solution.append(selected_candidate)
    solution.sort()
    #calculate the LDOI of the current solution
    if tuple(solution) not in TDOI_BFS:
      update_DOI(G_expanded,solution,TDOI_BFS,flagged,potential)
    #update candidate nodes
    candidates.remove(selected_candidate)
    negation_selected_candidate = Negation_in_expanded(selected_candidate)
    if negation_selected_candidate in candidates:
      candidates.remove(negation_selected_candidate)
    candidates = [x for x in candidates if x not in TDOI_BFS[tuple(solution)]]
    #check whether the target set is reached
    if Target <= TDOI_BFS[tuple(solution)]:
      break
    if len(candidates) == 0:
      fail_in_construction = True
      break
  return solution, fail_in_construction


def solution_reduction(G_expanded, Target, solution, TDOI_BFS, flagged, potential):
  '''
  The local search phase of the GRASP algorithm to solve the target control problem.
  The pseudocode is available in the Algorithm Table 2 of the appendix of the paper titled "Target Control in Logical Models Using the Domain of Influence of Nodes"
  to be publisehd in the Frontiers of Physiology.

  Parameters
  ----------
  G_expanded: the input expanded_graph as a DiGraph object
  Target : the target set
  solution : the candidate solution to be reduced during the local search phase
  TDOI_BFS : a dictionary to store the LDOI of each node state set
  flagged : a dictionary to store whether the LDOI of each node state set contains the negation of any target node
  potential : a dictionary to store the composite nodes visited but not included during the search process for LDOI for each node state set

  Returns
  -------
  reduced_solution: a list of node states as a control solution reduced from the candidate solution
  '''
  reduced_solution = solution[:]
  if len(reduced_solution) == 1:
    return reduced_solution
  #we remove node state from solution in a random order to be able to generate different final solution from the same candidate solution
  randomized_reduced_solution = solution[:]
  random.shuffle(randomized_reduced_solution)
  #iteratively remove node (state) from the candidate solution
  for node in randomized_reduced_solution:
      temp_solution = reduced_solution[:]
      temp_solution.remove(node)
      if tuple(temp_solution) not in TDOI_BFS:
        update_DOI(G_expanded, temp_solution, TDOI_BFS, flagged, potential)
      if Target <= TDOI_BFS[tuple(temp_solution)]:
        reduced_solution = temp_solution[:]
  return reduced_solution


def construct_initial_candidates(G_expanded, Target, forbidden_nodes=[], avail_nodes=[], forbidden_node_states=[]):
  '''
  The initialization function inside construction phase of the GRASP algorithm to solve the target control problem.


  Parameters
  ----------
  G_expanded: the input expanded_graph as a DiGraph object
  Target : the target set
  forbidden_nodes : the nodes that are forbiden to be used, both the positive states and the negative states will be forbidden
  avail_nodes : the nodes that are available to be used, both the positive states and the negative states will be available
  forbidden_node_states : the node states that are forbiden to be used

  notice different meaning for the default value for forbidden_nodes and avail_nodes
  when forbidden_nodes is empty, we do not forbid any node
  when avail_nodes is empyty, we allow all nodes

  Returns
  -------
  a list of candidate node states
  '''

  Target_duplicate=set()
  for node in Target:
    Target_duplicate.add(node)
    Target_duplicate.add(Negation_in_expanded(node))
  initial_candidates=set([x for x in G_expanded.nodes() if '_' not in x and x not in Target_duplicate])
  for node in forbidden_nodes:
    initial_candidates.discard(node)
    initial_candidates.discard(Negation_in_expanded(node))
  for node in forbidden_node_states:
    initial_candidates.discard(node)
  if not avail_nodes:
    return [x for x in initial_candidates]
  else:
    avail_candidates=set(avail_nodes[:])
    for node in avail_nodes:
      avail_candidates.add(Negation_in_expanded(node))
    return [x for x in initial_candidates&avail_candidates]

def write_node_property(nodes,filename,property_dict):
  f=open(filename,'a')
  for node in nodes:
    f.write(str(node)+' : '+str(property_dict[tuple([node])])+'\n')
  f.close()

def read_TDOI(TDOI_filename,TDOI_BFS):
  '''
  Read a pre-written LDOI record file with the name TDOI_filename and store into TDOI_BFS
  '''
  f=open(TDOI_filename,'r')
  lines=f.readlines()
  f.close()
  for line in lines:
    node, TDOI = line.split(' : ')
    TDOI=TDOI.strip()[5:-2]
    if TDOI=='':
      TDOI_list=[]
    elif ',' not in TDOI:
      TDOI_list=[TDOI.strip().strip("'")]
    else:
      TDOI_list=[x.strip().strip("'") for x in TDOI.split(',')]
    TDOI_BFS[tuple([node])]=set(TDOI_list)
  return

def read_flagged(flagged_filename,flagged):
  '''
  Read a pre-written flagged record file with the name flagged_filename and store into flagged
  '''
  f=open(flagged_filename,'r')
  lines=f.readlines()
  f.close()
  for line in lines:
    node, flag =line.split(' : ')
    flagged[tuple([node])]=ast.literal_eval(flag.strip())
  return

def read_potential(potential_filename,potential):
  '''
  Read a pre-written potential record file with the name potential_filename and store into potential
  '''
  f=open(potential_filename,'r')
  lines=f.readlines()
  f.close()
  for line in lines:
    node, potential_str = line.split(' : ')
    potential_str=potential_str.strip()[1:-1]
    if potential_str=='':
      potential_list=[]
    elif ',' not in potential_str:
      potential_list=[potential_str.strip().strip("'")]
    else:
      potential_list=[x.strip().strip("'") for x in potential_str.split(',')]
    potential[tuple([node])]=potential_list
  return

def update_single_DOI(G_expanded, result_filename, TDOI_BFS, flagged, potential):
  '''
  The pre-processing step of the GRASP algorithm to solve the target control problem.
  The function calculates the LDOI, whether this node is an incompatible intervention and composite nodes attached to the LDOI
  as recorded in TDOI_BFS, flagged, potential
  This function also output three files to record the above three quantites.
  The three files share the first part of the filename as result_filename.

  Parameters
  ----------
  G_expanded : the input expanded_graph as a DiGraph object
  result_filename : the prefix of the output file names
  TDOI_BFS : a dictionary to store the LDOI of each node state set
  flagged : a dictionary to store whether the LDOI of each node state set contains the negation of any target node
  potential : a dictionary to store the composite nodes visited but not included during the search process for LDOI for each node state set

  Returns
  -------
  None
  '''
  checklist = [result_filename+suffix for suffix in ['_TDOI.txt', '_flagged.txt', '_potential.txt']]
  TDOI_filename, flagged_filename, potential_filename = checklist
  result_exist= all([os.path.isfile(item) for item in checklist])
  if not result_exist:
    normal_nodes=[node for node in G_expanded.nodes() if '_' not in node]
    for node in normal_nodes:
        update_DOI(G_expanded,[node],TDOI_BFS,flagged,potential)
    write_node_property(normal_nodes,TDOI_filename,TDOI_BFS)
    write_node_property(normal_nodes,flagged_filename,flagged)
    write_node_property(normal_nodes,potential_filename,potential)
  else:
    read_TDOI(TDOI_filename,TDOI_BFS)
    read_flagged(flagged_filename,flagged)
    read_potential(potential_filename,potential)
  return

def GRASP_target_control(G_expanded, Target, max_itr, TDOI_BFS, flagged, potential, custom_score_index='0', forbidden_nodes=[], avail_nodes=[], forbidden_node_states=[], single_solution_preferred = False):
  '''
  The GRASP algorithm to solve the target control problem.
  The pseudocode is available in the paper titled "Target Control in Logical Models Using the Domain of Influence of Nodes"
  to be publisehd in the Frontiers of Physiology.

  Notice in this implementation, if the algorithm found any single node

  Parameters
  ----------
  G_expanded: the input expanded_graph as a DiGraph object
  Target : the target set
  max_itr : the number of iterations that are repeated to find the solution
  TDOI_BFS : a dictionary to store the LDOI of each node state set
  flagged : a dictionary to store whether the LDOI of each node state set contains the negation of any target node
  potential : a dictionary to store the composite nodes visited but not included during the search process for LDOI for each node state set
  candidates_score_index : the chosen heuristic greedy function,
                           '0' represents the size of LDOI, '1' represents the size of composite nodes attached to the LDOI
                           '2' represents the linear combination of the two above with the equal weight
                           '3' represents the size of LDOI with penalization, penalized by multiplied with -1 for those containing negation of any target
                           '4' represents the size of LDOI with penalization, penalized by being subtracted by maximum LDOI for those containing negation of any target
  forbidden_nodes : the nodes that are forbiden to be used, both the positive states and the negative states will be forbidden
  avail_nodes : the nodes that are available to be used, both the positive states and the negative states will be available
  forbidden_node_states : the node states that are forbiden to be used
  single_solution_preferred : A boolean value to indicate whether we only found single compatible solutions.

  notice different meaning for the default value for forbidden_nodes and avail_nodes
  when forbidden_nodes is empty, we do not forbid any node
  when avail_nodes is empyty, we allow all nodes

  Returns
  -------
  solutions: a dictionay maps a tuple to a integer value
    the first element of the tuple is the solution as a sorted tuple
    the second element of the tuple is a Boolean value to indicate whether the LDOI of the solution is incompatible or not. True means incompatible.
    The integer number is the frequency this solution is obtained during the max_itr iterations.
  A Boolean value to indicate whether we found any solution
  '''
  initial_candidates=construct_initial_candidates(G_expanded,Target, forbidden_nodes, avail_nodes, forbidden_node_states)
  solutions=defaultdict(int)
  singlenode_solution_found, compatible_singlenode_solution_found = False, False

  #if single solution is preferred, a search among single solution is performed
  if single_solution_preferred:
    for node in initial_candidates:
      if Target <= TDOI_BFS[tuple([node])]:
        solutions[(tuple([node]),flagged[tuple([node])])]+=1
        singlenode_solution_found=True
        if not flagged[tuple([node])]:
          compatible_singlenode_solution_found=True

  #second step, if single node can not solve the problem, using GRASP algorithm to find a solution
  if not compatible_singlenode_solution_found:

    candidates_score={}
    if custom_score_index=='0':
        candidates_score={candidate:custom_score1(candidate,measure=TDOI_BFS) for candidate in initial_candidates}
    elif custom_score_index=='1':
        candidates_score={candidate:custom_score1(candidate,measure=potential) for candidate in initial_candidates}
    elif custom_score_index=='2':
        candidates_score={candidate:custom_score2(candidate,measure1=TDOI_BFS,measure2=potential) for candidate in initial_candidates}
    elif custom_score_index=='3':
        candidates_score={candidate:custom_score3(candidate,TDOI_BFS,Target) for candidate in initial_candidates}
    elif custom_score_index=='4':
        maximum_TDOI=max([custom_score1(candidate,measure=TDOI_BFS) for candidate in initial_candidates])
        candidates_score={candidate:custom_score4(candidate,TDOI_BFS,Target,maximum_TDOI) for candidate in initial_candidates}
    else:
        raise ValueError('custome_score_index does not have the right value!')

    for itr in range(max_itr):

      candidates=initial_candidates[:]
      #phase 1: construct a greedy randomized solution
      solution,fail_in_construction=solution_construction(G_expanded,Target,candidates,candidates_score,TDOI_BFS,flagged,potential)

      if fail_in_construction:
        continue
      #print itr,solution

      #phase 2: local search, remove node by node that have redudant solution
      reduced_solution=solution_reduction(G_expanded,Target,solution,TDOI_BFS,flagged,potential)

      #print itr,reduced_solution

      #append solution
      solutions[(tuple(sorted(reduced_solution)),flagged[tuple(reduced_solution)])]+=1
  return solutions,len(solutions)>0



