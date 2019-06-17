'''
This file contains functions used to do pre-processing related to Boolean network models.
The code is written by Gang Yang, Department of Physics, Penn State University if not specified.
'''

from collections import deque
import networkx as nx
import qm       #Quine-McCluskey algorithm


def Getfunc(inputlist, output, onlist, prefix='n', suffix='n', equal_sign='*='):
  '''
  Return a Boolean regulatory rule in the disjuctive normal form (in the format of Booleannet)
  based on given truth table.
  For more details about Booleannet, refer to https://github.com/ialbert/booleannet.
  It is essentially a wrapper for the Quine-McCluskey algorithm to meet the format of Booleannet.


  Parameters
  ----------
  inputlist : a list of nodenames of all the parent nodes/ regulators
  output    : the nodename of the child node
  onlist    : a list of all configurations(strings) that will give a ON state of the output based on rule
              e.g. the onlist of C = A OR B will be ['11','10','01']
  prefix='n': prefix to encode the node name to avoid one node's name is a part of another node's name
  suffix='n': suffix to encode the node name to avoid one node's name is a part of another node's name
              e.g. node name '1' will become 'n1n' in the returned result
  equal_sign: the equal sign of the rule in the returned result, whose default value follows the Booleannet format

  Returns
  -------
  The Boolean rule in the disjuctive normal form.


  References
  ----------
  QM code by George Prekas.
  '''
  #pre-processing to meet the format requirement of QM and Booleannet
  inputs=[prefix+str(x)+suffix for x in inputlist]
  inputs.reverse()
  output=prefix+str(output)+suffix+equal_sign
  onindexlist=[int(x,2) for x in onlist]
  if len(inputs)>0:
    qmtemp= qm.QM(inputs)
    temprule=qmtemp.get_function(qmtemp.solve(onindexlist,[])[1])
    temprule=temprule.replace('AND','and')
    temprule=temprule.replace('OR','or')
    temprule=temprule.replace('NOT','not')
  else:
    temprule=output[:-2]
  return output+temprule+'\n'

def write_Boolean_rules(Ginput,filename):
  '''
  Generate a Boolean rule file (in the format of Booleannet) based on a Boolean network model as a DiGraph object.
  All Boolean rules are written in the disjuctive normal form.

  Parameters
  ----------
  Ginput: a DiGraph object, the Boolean rules are given in the node property of the DiGraph object.
          The Boolean network model is a DiGraph object in the output format of form_network().
          The Boolean network model can be generated through form_network function by reading a text file in the Booleannet format.
  filename: the file name that the generated Boolean rule will be written

  Returns
  -------
  This function return None. It generated a Boolean rule file for the given Boolean network model.
  '''

  G=Ginput.copy()
  f=open(filename,'a')
  f.write('#BOOLEAN RULES\n')
  for node in G.nodes():
    ON_list=[x for x in G.node[node]['update_rules'].keys() if G.node[node]['update_rules'][x]==1]
    rule=Getfunc(G.node[node]['update_nodes'],node,ON_list)
    f.write(rule)
  f.close()
  return

def draw_gml(G,filename):
  '''
  Generate a file in the gml format which can be opened in yEd for a given Boolean network model as a DiGraph object.
  '''
  Gprime=nx.DiGraph()
  Gprime.add_edges_from(G.edges())
  nx.write_gml(Gprime,filename)
  return

def Get_expanded_network(Gread,prefix='n',suffix='n',equal_sign='*='):
  '''
  Return the expanded network for a given Boolean network model.
  The Boolean network model is a DiGraph object in the output format of form_network().
  The Boolean network model can be generated through form_network function by reading a text file in the Booleannet format.
  The Boolean rules will first be converted to a disjuctive normal form before generating the expanded network.

  Parameters
  ----------
  Gread     : the given Boolean network model
  prefix='n': prefix to encode the node name to avoid one node's name is a part of another node's name
  suffix='n': suffix to encode the node name to avoid one node's name is a part of another node's name
              e.g. node name '1' will become 'n1n' in the returned result
  equal_sign: the equal sign of the rule in the returned result, whose default value follows the Booleannet format

  Returns
  -------
  The expanded network for the given Boolean network model.


  '''
  G_expand=nx.DiGraph()
  rules = []
  #first write rules for negation nodes
  negation_rules=[]
  expanded_nodes=set()
  for node in Gread.nodes():
      ON_list=[x for x in Gread.node[node]['update_rules'].keys() if Gread.node[node]['update_rules'][x]==1]
      rule=Getfunc(Gread.node[node]['update_nodes'],node,ON_list,prefix=prefix,suffix=suffix,equal_sign=equal_sign)
      OFF_list=[x for x in Gread.node[node]['update_rules'].keys() if Gread.node[node]['update_rules'][x]==0]
      negation_rule='~'+Getfunc(Gread.node[node]['update_nodes'],node,OFF_list,prefix=prefix,suffix=suffix,equal_sign=equal_sign)
      rules.append(rule)
      negation_rules.append(negation_rule)
      expanded_nodes.add(rule.split('*=')[0])
      expanded_nodes.add(negation_rule.split('*=')[0])
  #then for each line in the rules, construct Boolean network
  composite_nodes=[]
  rules.extend(negation_rules)
  for line in rules:
      child, update_rule=line.split('*=')
      update_rule=update_rule.strip()
      if update_rule[0]=='(' and update_rule[-1]==')':
        update_rule=update_rule[1:-1]
      #single parent situation
      if child[0]=='~':
        normal_child=child[1:]
      else:
        normal_child=child[:]
      normal_child=normal_child[len(prefix):len(normal_child)-len(suffix)]
      #deal with source node situation
      if not Gread.node[int(normal_child)]['update_nodes']:
        G_expand.add_node(child)     #maybe this do not need to be done
      else:
        if 'or' in update_rule:
          parents=update_rule.split(' or ')
        else:
          parents=[update_rule]
        parents.sort()
        for parent in parents:
          parent=parent.replace('not ','~').replace('(','').replace(')','')
          if 'and' in parent:
            composite_node=parent.replace(' and ','_')
            composite_nodes.append(composite_node)
            G_expand.add_edge(composite_node,child)
            for component in composite_node.split('_'):
              G_expand.add_edge(component,composite_node)
          else:
            G_expand.add_edge(parent,child)
  return G_expand.copy()

def Get_reduced_network(Gread):
  '''
  Return the reduced network for a given Boolean network model.
  The Boolean network model is a DiGraph object in the output format of form_network().
  The Boolean network model can be generated through form_network function by reading a text file in the Booleannet format.
  The algorithm iteratively remove nodes whose node state is fixed.

  Parameters
  ----------
  Gread     : the given Boolean network model


  Returns
  -------
  the reduced Boolean network for the input


  '''
  #doing network reduction here
  G_reduced=Gread.copy()
  fixed=set()
  queue=deque()
  #collecting all the source nodes
  for node in Gread.nodes():
    if not Gread.node[node]['update_nodes']:
      source_node_value_list=[int(x) for x in set(Gread.node[node]['update_rules'].values())]
      assert len(source_node_value_list)==1
      source_node_value=source_node_value_list[0]
      queue.append((node,source_node_value))


  #BFS search to reduce source node, frontier kept in queue, visited nodes kept in fixed as a set
  while len(queue)>0:
    vertex=queue.popleft()
    if vertex not in fixed:
      fixed.add(vertex)
      untested_successors=[x for x in set(G_reduced.successors(vertex[0]))-set([y[0] for y in fixed])-set([z[0] for z in queue])]
      for node in untested_successors:
        fixed_input_index=G_reduced.node[node]['update_nodes'].index(vertex[0])
        reduced_updated_rules={}
        for key,value in G_reduced.node[node]['update_rules'].iteritems():
          if key[fixed_input_index]==str(vertex[1]):
            if fixed_input_index<len(G_reduced.node[node]['update_nodes'])-1:
              new_key=key[:fixed_input_index]+key[fixed_input_index+1:]
            else:
              new_key=key[:fixed_input_index]
            reduced_updated_rules[new_key]=value
        G_reduced.node[node]['update_rules']=reduced_updated_rules
        G_reduced.node[node]['update_nodes'].remove(vertex[0])
        if len(set(reduced_updated_rules.values()))==1:
          reduced_value=reduced_updated_rules.values()[0]
          queue.append((node,reduced_value))
      #in the end remove the fixed node from the graph
      G_reduced.remove_node(vertex[0])
  return G_reduced.copy()

def read_expanded_network(filename,expanded_graphname):
  '''
  Read the expanded network generated by Jorge's Stable Motif program (written in Java).
  More details about stable motif program can be found in https://github.com/jgtz/StableMotifs
  This works only if the Boolean network does not have any source node.
  If the network has source nodes, split the expanded network generated by Jorge's program.


  Parameters
  ----------
  filename     : filename of the expanded network generated by Jorge program.


  Returns
  -------
  the expanded network as a DiGraph object,
  and a dictionary maps the node name from the notation in the expanded network to the original name


  '''
  #read the file generated by the expanded network java program
  f = open(filename,'r')   #Remember to change the name back to sample_network.txt
  lines = f.readlines()
  f.close()

  #find the position that defines the initial index and final index
  Node_ini_index=lines.index('Expanded node names\n')+1
  Edge_ini_index=lines.index('Adj list\n')+1
  Edge_end_index=lines.index('Finding stable motifs in this network...\n')-1

  #collect nodes information
  Nodes_info=lines[Node_ini_index:Edge_ini_index-2]
  Edges_info=lines[Edge_ini_index:Edge_end_index]
  Nodes_info=[x[:-1] for x in Nodes_info]
  Edges_info=[x[:-1] for x in Edges_info]
  Nodes={}
  Node_list=[]
  Edge_list=[]

  for item in Nodes_info:
    [tempindex,tempname] = item.split('\t')
    Nodes[tempindex]=tempname

  Node_list=Nodes.keys()

  for item in Edges_info:
    tempitem = item.split('\t')
    if len(tempitem)>1:
      #form the edge list
      tempedgelist=[(tempitem[0],tempitem[x]) for x in range(1,len(tempitem))]
      for tempedge in tempedgelist:
        Edge_list.append(tempedge)
    else:
      #there is no edge
      continue

  #construct the graph
  G=nx.DiGraph()
  G.add_edges_from(Edge_list)
  check=set(G.nodes())==set(Node_list)

  #save the graphical representation of the expanded network
  Gprime=nx.DiGraph()
  edgelist=[(Nodes[x[0]],Nodes[x[1]]) for x in G.edges()]
  Gprime.add_edges_from(edgelist)
  nx.write_gml(Gprime,expanded_graphname)

  return G.copy(),Nodes


def form_network(rules,sorted_nodename=True):
    '''
    Takes as input a list of rules in the format of Booleannet(e.g. sample_network.txt)

    Outputs a networkx DiGraph with node properties:
        'update_nodes': a list of regulating nodes
        'update_rules': a dictionary with binary strings as keys, corresponding
                        to the possible states of the update_nodes, and integers
                        as values, corresponding to the state of the node given
                        that input.

    Note that nodes are identified by their position in a sorted list of
    user-provided node names.

    Notice the name of one node should not be part of other node

    The code is written by Colin Campbell.
    '''
    def clean_states(x):
        #cleans binary representation of node input states
        out=x[2:]                                                               # Strip leading 0b
        return '0'*(len(inf)-len(out))+out                                      # Append leading 0's as needed

    stream = [x.rstrip('\n') for x in rules if x != '\n' and x[0]!='#']         # Remove comments and blank lines
    #I made a slight change here so that the code will be compatible with Jorge's Java code
    if sorted_nodename:
      nodes = sorted([x.split('*=',1)[0] for x in stream])                    # Generate a sorted list of node names
    else:
      nodes = [x.split('*=',1)[0] for x in stream]

    g = nx.DiGraph()
    g.graph['knockout'] = None                                                  # At creation, no node is flagged for knockout or overexpression
    g.graph['express'] = None

    for n in xrange(len(stream)):
        node = stream[n].split('*=')[0]
        rule = stream[n].split('*=')[1]
        rule = rule.replace(' AND ',' and ')                                    # Force decap of logical operators so as to work with eval()
        rule = rule.replace(' OR ',' or ')
        rule = rule.replace(' NOT ',' not ')
        if stream[n].find('True') >= 0 or stream[n].find('False') >= 0:         # For always ON or always OFF nodes
            g.add_node(nodes.index(node))                                                       # We refer to nodes by their location in a sorted list of the user-provided node names
            g.node[nodes.index(node)]['update_nodes'] = []
            g.node[nodes.index(node)]['update_rules'] = {'':str(int(eval(rule)))}
            continue

        inf = rule.split(' ')                                                   # Strip down to just a list of influencing nodes
        inf = [x.lstrip('(') for x in inf]
        inf = [x.rstrip(')') for x in inf]
        #The sort ensures that when we do text replacement (<node string>->'True' or 'False') below in this fn, we avoid problems where node 1 is a substring of node 2 (e.g. NODE1_phosph and NODE1)
        inf = sorted([x for x in inf if x not in ['','and','or','not']],key=len,reverse=True)
        inf = [x for x in set(inf)]    #04/16/2016 to allow one variable appear twice in the rule like a XOR rule

        #mod
        for i in inf: g.add_edge(nodes.index(i),nodes.index(node))                              # Add edges from all influencing nodes to target node
        g.node[nodes.index(node)]['update_nodes'] = [nodes.index(i) for i in inf]     #mod
        g.node[nodes.index(node)]['update_rules'] = {}      #mod

        bool_states = map(bin,range(2**len(inf)))
        bool_states = map(clean_states,bool_states)
        for j in bool_states:
            rule_mod = rule[:]
            for k in range(len(j)):
                if j[k]=='0':
                    rule_mod=rule_mod.replace(nodes[g.node[nodes.index(node)]['update_nodes'][k]],'False')      # Modify the rule to have every combination of True, False for input nodes
                else: rule_mod=rule_mod.replace(nodes[g.node[nodes.index(node)]['update_nodes'][k]],'True')
            g.node[nodes.index(node)]['update_rules'][j] = int(eval(rule_mod))                                  # Store outcome for every possible input
            #mod
    return g,nodes

