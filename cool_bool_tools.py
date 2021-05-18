import networkx as nx
import numpy as np
import BooleanDOI_processing as BDOIp
import BooleanDOI_TargetControl as BDOItc
import BooleanDOI_DOI as BDOI

def partial_overlap_of_states(state_big,state_part):
    '''
    Returns the overlap of two states, where state_part is a subset or equal to of state_big
    '''
    overlap=0
    for module_node in state_part.keys():
        if state_part[module_node]==state_big[module_node]:
            overlap+=1

    return float(overlap)/len(dict(state_part))

def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z

def read_rules_text(model_name):
    '''
    docstring
    '''


    rules_file='%s.txt'%model_name
    with open(rules_file,'r') as f:
        rules=f.read()
    rules=rules.replace('#BOOLEAN RULES','')
    return rules

def read_in_attractors(model_name):

    '''docstring'''

    import pandas as pd
    attractors_file_name='%s-QuasiAttractors'%model_name+'.txt'
    import pandas as pd
    df=pd.read_csv(attractors_file_name,sep='\t',index_col=-1) #index_col=-1 we cheat by making the extra tabs at the end of the rows in the files as the index column, and then resetting it to just natural incrementing integers
    df.reset_index(inplace=True,drop=True)
    a=df.to_dict('index')
    attractors={'Attractor%d'%i:a[i] for i in a.keys()}
    return attractors

def partial_overlap_of_states(state_big,state_part):
    '''
    Returns the overlap of two states, where state_part is a subset or equal to of state_big
    '''
    overlap=0
    for module_node in state_part.keys():
        if state_part[module_node]==state_big[module_node]:
            overlap+=1

    return float(overlap)/len(dict(state_part))


def plot_state_succession(states,state_labels=None,title=None, nodes=None, x_fontsize=10,state_colors=['blue','orange'],y_fontsize=14):

    '''
    docstring
    '''
    from matplotlib import pyplot as plt
    from matplotlib import colors
    import numpy as np
    if nodes==None:
    	nodes=states[0].keys()
    state_transition=[]
    for s in states:
        state_transition.append([s[k] for k in nodes])
    cmap = colors.ListedColormap(state_colors)
    steps=len(state_transition)
    plt.figure(figsize=(len(nodes),steps))
    plt.imshow(state_transition, interpolation='none',cmap=cmap)
    ax = plt.gca()
    plt.xticks(range(len(nodes)),nodes, fontsize=x_fontsize)
    if state_labels==None:
        plt.yticks(range(steps),fontsize=y_fontsize)
    else:
        plt.yticks(range(steps),state_labels,fontsize=y_fontsize)
    if title!=None:
        plt.title(title)
    ax.set_yticks(np.arange(-.5, steps,1),minor=True)
    plt.grid(which='minor', color='black', linestyle='-', linewidth=1)
    plt.show()


class TransGraph(object):
    """
    Represents a transition graph
    """
    def __init__(self, logfile, verbose=False):
        import networkx as nx
        self.graph = nx.MultiDiGraph(  )
        self.fp = open( logfile, 'wt')
        self.verbose = verbose
        self.seen = set()
        self.store = dict()
        self.colors = dict()

    def add(self, states, times=None):
        "Adds states to the transition"

        # generating the fingerprints and sto
        times = times or range(len(states))
        fprints = []
        for state in states:
            if self.verbose:
                fp = state.bin()
            else:
                fp = state.fp()
            fprints.append( fp )
            self.store[fp] = state

        self.fp.write( '*** transitions from %s ***\n' % fprints[0] )

        for head, tail, tstamp in zip(fprints, fprints[1:], times ):
            pair = (head, tail)
            self.fp.write('T=%s: %s->%s\n' %  (tstamp, head, tail) )
            if pair not in self.seen:
                self.graph.add_edge(head, tail)
                self.seen.add(pair)

    def save(self, fname, colormap={}):
        "Saves the graph as gml"
        write_gml(graph=self.graph, fname=fname, colormap=colormap)

        self.fp.write( '*** node values ***\n' )

        # writes the mapping
        first = self.store.values()[0]
        header = [ 'state' ] + first.keys()
        self.fp.write( util.join(header) )

        for fprint, state in sorted( self.store.items() ):
            line = [ fprint ]  + map(int, state.values() )
            self.fp.write( util.join(line) )

def get_async_STG_all_states(model_async,node_threshold=0):
    '''
    Given an asyncronous model the function initiates all 2**N initial states, where N is the number of nodes, and performs an update of each node then reinitializes to the previous state,
    adding each transition edge to a state transition graph (STG).It also generates a fingerprint dictionary (fp_state_dict) which associates detailed states with integer id-s that are used as
    node id-s in the STG.
    WARNING! Due to the exponential growth of the state space, with the size of the network, this can be extremely demanding for larger networks
    One can also specify a node threshold, where the function perfoms a filtering based on the number of times a state was visited.

    Inputs: model_async - boolean2 type network model initiated as asynchronous
        node_threshold - int (default=0)

    Returns: g - networx multigraph object
         edge_occurances - dict
         state_occurances - dict
         fp_state_dict - dict
    '''
    from itertools import combinations, product
    nodes=model_async.nodes.copy()
    value_combinations=list(product([False,True],repeat=len(nodes)))


    fp_state_dict={}
    fp_states=[]
    j=0
    tg=TransGraph(logfile='async_states.txt')

    for i in value_combinations:

        start_state=dict(zip(nodes,i))

        for node in nodes:

            model_async.initialize(lambda snode: start_state[snode])
            #evolution=evolve_from_state(model_async,n,initial_state=start_state,shuffler=shuffler)
            model_async.iterate(1,shuffler=lambda lines: controlled_async_pick(lines,node))

            fp_state_dict=merge_two_dicts(fp_state_dict,dict(zip(model_async.fp(),model_async.states)))
            fp_states+=model_async.fp()
            j+=1
            tg.add(model_async.states)

    G=tg.graph

    from collections import Counter
    state_occurances=Counter([fp_states[i] for i in range(1,len(fp_states),2)])
    edge_occurances=Counter([(fp_states[i], fp_states[i+1]) for i in range(0,len(fp_states)-1,2)])
    filtered_nodes=[]
    for n in G.nodes():
        if state_occurances[n]>node_threshold:
            filtered_nodes.append(n)

    g=nx.subgraph(G,filtered_nodes).copy()
    return g,edge_occurances,state_occurances,fp_state_dict #g is a nx.MultiGraph, so having the edge weights and node weights in the graph object is a bit tedious

def controlled_async_pick(lines, node):
    for l in lines:
        if l.split('*')[0].strip() == node:
            return [l]
    raise ValueError('Node not in model!')

def STG_backbone_based_on_proxy_nodes(G, PNs, cutoff=16):

    '''
    Given a state transition graph and a set of proxy nodes from the graph the algorithm determines the collapsed
    path probabilities, by summing up the probabilities of paths between all pairs of proxy nodes, such that no
    no other proxy node is included in the path.

    Inputs: G - networkx DiGraph
            PNs - list of proxy nodes (from G)
            cutoff - the meximum length of paths

    Output: G_BB - a networkx representation of the "backbone", where the nodes are the proxies and the edges the
            potential transitions between them
            G_BB_edge_weights - the proabaility of transitions between proxy nodes
    '''

    edge_weights = {}
    for e in G.edges():
        edge_weights[e]=1./G.out_degree(e[0])

    G_BB=nx.DiGraph()
    G_BB_edge_weights={}
    G_BB=nx.complete_graph(PNs,create_using=G_BB)
    PN_pair_path_edges={}
    all_path_edges=[]
    for source_PN, target_PN in G_BB.edges():
        PN_set=set(PNs)
        PN_set.remove(source_PN)
        PN_set.remove(target_PN)
        path_edges=[]
        meta_edge_probability=0
        for path in nx.all_simple_paths(G,source_PN,target_PN, cutoff=cutoff):
            if len(set(path).intersection(PN_set))==0:
                edges_on_the_path=zip(path[:-1],path[1:])
                path_probability=np.prod([edge_weights[e] for e in edges_on_the_path])
                meta_edge_probability+=path_probability

        if meta_edge_probability!=0:
            G_BB_edge_weights[(source_PN, target_PN)]=meta_edge_probability

    for source_PN, target_PN in list(G_BB.edges()):
        if (source_PN, target_PN) not in G_BB_edge_weights.keys():
            G_BB.remove_edge(source_PN, target_PN)

    return G_BB,G_BB_edge_weights


def read_and_create_extended_network(network_file, relabel=True):

    '''
    docstring
    '''

    with  open(network_file,'r') as f:
        node_list=[]
        lines = f.readlines()
        for i in lines:
            if '*=' in i:
                node_list.append(i.split('*=')[0].strip())

    Gread, readnodes = BDOIp.form_network(lines, sorted_nodename = False)

    G_expanded = BDOIp.Get_expanded_network(Gread)

    if relabel:

        base_mapping={}
        for i in range(len(node_list)):
            base_mapping['n%dn'%i]=node_list[i]

        node_label_mapping={}
        for n in G_expanded.nodes():
            n_new=n[:]
            for indexed_name, new_name in base_mapping.items():
                n_new=n_new.replace(indexed_name,new_name)
            node_label_mapping[n]=n_new
        G_expanded=nx.relabel_nodes(G_expanded,node_label_mapping)

    return G_expanded

def state_nodes_with_conflict(BFS):
    '''
    docstring
    '''
    return list(set([i for i in BFS[0] if '_' not in i]+BFS[2]))
def state_nodes_only(BFS):
    '''
    docstring
    '''
    return list(set([i for i in BFS[0] if '_' not in i]))

def _not(s):
    '''
    docstring
    '''
    if s.strip()[0]=='~':
        return s.replace('~','')
    else:
        return '~'+s

def are_subsets_consistent(a,b):

    '''
    docstring
    '''
    a_all=[]
    b_all=[]
    for i in a:
        a_all+=i.split('_')
    for i in b:
        b_all+=i.split('_')
    for i in a_all:
        if _not(i) in b_all:
            return False
    return True

def expanded_network_consistent_simple_cycles(G):
    """Find consistent simple cycles (elementary circuits) of an expanded network.


    An simple cycle, or elementary circuit, is a closed path where no
    node appears twice, except that the first and last node are the same.
    Two elementary circuits are distinct if they are not cyclic permutations
    of each other.
    A consistent cycle is one that doesn't contain any contradicting virtual or composite nodes


    This is a nonrecursive, iterator/generator version of Johnson's
    algorithm [1]_.  There may be better algorithms for some cases [2]_ [3]_.

    Parameters
    ----------
    G : NetworkX DiGraph
       A directed graph

    Returns
    -------
    cycle_generator: generator
       A generator that produces elementary cycles of the graph.  Each cycle is
       a list of nodes with the first and last nodes being the same.


    Notes
    -----
    The implementation is a modified version of [1] bades on the NetworkX implementaiton:
    https://networkx.github.io/documentation/networkx-1.9/_modules/networkx/algorithms/cycles.html#simple_cycles

    The time complexity is O((n+e)(c+1)) for n nodes, e edges and c
    consistent elementary circuits.


    References
    ----------
    .. [1] Finding all the elementary circuits of a directed graph.
       D. B. Johnson, SIAM Journal on Computing 4, no. 1, 77-84, 1975.
       http://dx.doi.org/10.1137/0204007

    """
    import cool_bool_tools as cbt
    from collections import defaultdict
    def _unblock(thisnode,blocked,B):
        stack=set([thisnode])
        while stack:
            node=stack.pop()
            if node in blocked:
                blocked.remove(node)
                stack.update(B[node])
                B[node].clear()

    # Johnson's algorithm requires some ordering of the nodes.
    # We assign the arbitrary ordering given by the strongly connected comps
    # There is no need to track the ordering as each node removed as processed.
    subG = type(G)(G.edges()) # save the actual graph so we can mutate it here
                              # We only take the edges because we do not want to
                              # copy edge and node attributes here.
    sccs = list(nx.strongly_connected_components(subG))
    while sccs:
        scc=sccs.pop()
        # order of scc determines ordering of nodes
        startnode = scc.pop()
        # Processing node runs "circuit" routine from recursive version
        path=[startnode]
        blocked = set() # vertex: blocked from search?
        closed = set() # nodes involved in a cycle
        blocked.add(startnode)
        B=defaultdict(set) # graph portions that yield no elementary circuit
        stack=[ (startnode,list(subG[startnode])) ]  # subG gives component nbrs
        while stack:
            thisnode,nbrs = stack[-1]
            if nbrs:
                nextnode = nbrs.pop()
#                    print thisnode,nbrs,":",nextnode,blocked,B,path,stack,startnode
#                    f=raw_input("pause")
                if nextnode == startnode:
                    yield path[:]
                    closed.update(path)
#                        print "Found a cycle",path,closed
                elif nextnode not in blocked and cbt.are_subsets_consistent(path,[nextnode]):
                    path.append(nextnode)
                    stack.append( (nextnode,list(subG[nextnode])) )
                    closed.discard(nextnode)
                    blocked.add(nextnode)
                    continue
            # done with nextnode... look for more neighbors
            if not nbrs:  # no more nbrs
                if thisnode in closed:
                    _unblock(thisnode,blocked,B)
                else:
                    for nbr in subG[thisnode]:
                        if thisnode not in B[nbr]:
                            B[nbr].add(thisnode)
                stack.pop()
#                assert path[-1]==thisnode
                path.pop()
        # done processing this node
        subG.remove_node(startnode)
        H=subG.subgraph(scc)  # make smaller to avoid work in SCC routine
        sccs.extend(list(nx.strongly_connected_components(H)))
