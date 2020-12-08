#from BitVector import BitVector
import numpy as np
from collections import defaultdict

class Graph:
    def __init__(self):
        self.edges = defaultdict(set)
        self.indexes = dict()
        self.nodes = dict()
        # self.bitEdges = []
        self.adj_mat = np.zeros((0, 0))
        #self.name = "graph_name"

    def add_edge(self, from_node, to_node):
        self.edges[from_node].add(to_node)
        self.edges[to_node].add(from_node)


    def build_adjmat(self, from_node, to_node):
        self.adj_mat[from_node][to_node] = 1
        self.adj_mat[to_node][from_node] = 1

    """
    def has_edge(self, from_node, to_node):
        return self.bitEdges[from_node*len(self)+to_node]

    """
    def has_edge(self, from_node, to_node):
        return self.adj_mat[from_node][to_node] == 1

    def num_edges(self):
        return sum([len(x) for x in self.edges.values()])

    """
    getNeighbors() takes a node N from a graph
    and returns a list of all neighbors of N.
    """
    def get_neighbors(self, node):
       return self.edges.get(node, [])

    def degree(self,node):
        return len(self.get_neighbors(node))

    def __repr__(self):
        return str(self.adj_mat)
    
    def __len__(self):
        return len(self.edges)
