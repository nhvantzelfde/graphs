
import sys
sys.path.insert(0, '../data_structures')
from heap import MinHeap
from hash_table import KeyValuePair

import random

class Edge(object):
    """ Class representing a single edge in a graph.

    Edges are directed edges, but graphs can add both (u,v) and (v,u) for an undirected edge.
    Vertices are assumed to be referenced by number, so u and v are integers.

    Attributes:
        orig = origin vertex number
        dest = destination vertex number
        weight = weight of the edge; default is 1 for an unweighted graph
    """

    def __init__(self, u, v, w = 1):
        """ Initializes an edge from u->v with optional weight w. """
        self.orig = u
        self.dest = v
        self.weight = w

    def __str__(self):
        """ Returns a string representation of the edge. """
        return "edge: " + str(self.orig) + "->" + str(self.dest) + ", wgt " + str(self.weight)

class Graph(object):
    """ Class representing a full graph, which can be grown one vertex or one edge at a time.

    Vertices are numbered started at 0; edges are represented by Edge objects.
    Internal attributes should not be manipulated directly; use addVertex and addEdge instead.
    As of now, edges and vertices cannot be removed once added, though this functionality could be added later.

    Attributes:
        _directed = True if the graph has directed edges and False otherwise; set upon initialization, assumed not to change
        _adj = adjacency list of all edges, indexed by vertex number
        _edge_count = number of edges added to the graph; undirected graphs add (u,v) and (v,u) for each edge, but these counts as 1
        _neg_wgt_count = the number of negative-weight edges in the graph; certain functions won't work if there are any negative-weight edges
    """
        
    def __init__(self, directed = False):
        """ Initializes an empty graph, which is by default undirected. """
        self._directed = directed
        self._adj = []
        self._edge_count = 0
        self._neg_wgt_count = 0

    def __str__(self):
        """ Returns a string representation of the graph, which is the adjacency list. """
        return str(self._adj)

    def addVertex(self):
        """ Adds a vertex to the graph, with number n+1, where n is the current number of vertices. Returns the new vertex's number. """
        self._adj.append([])
        return len(self._adj)-1

    def addEdge(self, u, v, w = 1):
        """ Adds an edge to the graph from (u,v) with optional weight w. Enforces certain rules depending on graph type. Overwrites existing edge at (u,v) if any. """
        if u >= self.vertexCount() or v >= self.vertexCount():
            return # impossible to add edge - vertices don't exist
        if self.isUndirected() and u == v:
            return # do not allow self loops in undirected graphs
        if w == 0:
            return # do not allow zero-weight edges

        # track whether there are any negative weight edges in the graph
        if w < 0: self._neg_wgt_count += 1
        
        e = self.findEdge(u,v)
        if e:
            # overwrite existing edge
            if e.weight < 0: self._neg_wgt_count -= 1
            e.weight = w
        else:
            # create new edge
            self._adj[u].append(Edge(u,v,w))
            self._edge_count += 1

        # create edge from v->u as well if the graph is undirected
        if self.isUndirected():
            e = self.findEdge(v,u)
            if e:
                # overwrite existing edge
                e.weight = w
            else:
                # create new edge
                self._adj[v].append(Edge(v,u,w))

    def findEdge(self, u, v):
        """ Returns the Edge from u->v in the graph, if it exists, and None otherwise. """
        if u >= len(self._adj) or v >= len(self._adj):
            return None # vertices don't exist

        edges = self._adj[u]
        for e in edges:
            if e.dest == v:
                return e
        return None

    def isDirected(self):
        """ Returns True if the graph is a directed graph, and False otherwise. """
        return self._directed

    def isUndirected(self):
        """ Returns True if the graph is a undirected graph, and False otherwise. """
        return not self.isDirected()

    def hasNegativeWeight(self):
        """ Returns True if the graph has any negative weight edges, and False otherwise. """
        return self._neg_wgt_count != 0
    
    def asAdjacency(self):
        """ Returns a representation of the graph in adjacency list form. """
        return self._adj

    def asMatrix(self):
        """ Returns a representation of the graph in adjacency matrix form. For sparse graphs this will use more space than an adjacency list. """
        # create matrix of zeros
        m = [[]] * self.vertexCount()
        for i in range(len(m)):
            m[i] = [0] * self.vertexCount()

        # fill in weights for existing edges
        for i, edges in enumerate(self._adj):
            for e in edges:
                m[i][e.dest] = e.weight
        return m
    
    def vertexCount(self):
        """ Returns the number of vertices in the grahh. """
        return len(self._adj)

    def edgeCount(self):
        """ Returns the number of edges in the grahh. """
        return self._edge_count
   
    def BFS(self, v):
        """ Runs breadth-first search from the source vertex v. Returns a matrix with distances from v, or None if v does not exist. """
        if v >= self.vertexCount():
            return None # can't be run; vertex doesn't exist

        # initialize the vertex queue, the "visited" set s, and the distances from v
        q = MinHeap()
        s = set()
        distance = [float("inf")] * self.vertexCount()

        # insert the source vertex into q and set its distance as 0
        q.insert(KeyValuePair(0,v))
        s.add(v)
        distance[v] = 0

        # loop through all vertices, in order of distance from v, relaxing edges out of current vertex
        kvp = q.extractMin()
        while kvp:
            v = kvp.value
            d = kvp.key

            # relax all edges out of current vertex v
            edges = self._adj[v]
            for e in edges:
                if e.dest not in s:
                    s.add(e.dest)
                    distance[e.dest] = d + 1
                    q.insert(KeyValuePair(distance[e.dest],e.dest))

            # get next-lowest distance vertex
            kvp = q.extractMin()

        return distance

    def DFS(self):
        """ Runs depth-first search from all vertices in order, return the list of parent vertices, the topological sort order, and whether there is a cycle. """
        
        # the set of seen vertices, the parent list, and fin time list, and the topological sort order
        seen = set()
        par = [None] * self.vertexCount()
        fin = [None] * self.vertexCount()
        order = []

        aux = [False,1] # aux[0] = cycle / back edge exists; aux[1] = finish count

        # start DFS for nodes in numerical order, if they have not already been visited
        for node in range(len(self._adj)):
            if not fin[node]:
                seen.add(node)
                par[node] = "S" # "S" is a starting node for a DFS
                self._DFSHelper(node,seen,par,fin,order,aux)

        # order is calculated in order of finish; topoloical order is the reverse of this
        order.reverse()
        return (par, order, aux[0])

    def _DFSHelper(self, node, seen, par, fin, order, aux):
        """ Recursively calculates DFS from a given node. Modifies data in seen, par, fin, order and aux in place. """

        # scan all edges of a given node
        edges = self._adj[node]
        for e in edges:
            # recurse if the node hasn't been seen yet
            if e.dest not in seen:
                seen.add(e.dest)
                par[e.dest] = node
                self._DFSHelper(e.dest,seen,par,fin,order,aux)
            # otherwise, determine whether there is a back edge
            else:
                if not fin[e.dest]:
                    aux[0] = True # cycle exists, destination node not finished yet

        # set finish time and incremement; add to finished order
        fin[node] = aux[1]
        aux[1] += 1
        order.append(node)
        
    def topologicalSort(self):
        """ Returns a topological sort of the vertices, if the graph is a DAG, and None if the graph is not a DAG. """
        if self.isUndirected():
            return None # undirected graph, cannot sort

        # retrieve parent list, topological order, and whether a cycle exists from DFS
        p, o, c = self.DFS()
        
        if c:
            return None # cycle exists, cannot sort
        else:
            return o # graph is acyclic; topological sort is valid

    def hasCycle(self):
        """ Returns True if the graph contains at least one cycle, and False otherwise. """
        p, o, c = self.DFS()
        return c

    def isDAG(self):
        """ Returns True if the graph is a directed, acyclic graph, and False otherwise. """
        if self.isUndirected(): return False
        p, o, c = self.DFS()
        return not c
        
    def BellmanFord(self, v):
        """ Runs Bellman-Ford shortest-path algorithm from source v, returning the list of distances from v, or None if there is a negative-weight cycle. """
        d = [float("inf")] * self.vertexCount()
        d[v] = 0

        # relax each edge V-1 times; final time determines whether there is a negative-weight cycle
        for i in range(self.vertexCount()):
            for edges in self._adj:
                for e in edges:
                    if d[e.orig] + e.weight < d[e.dest]:
                        if i < (self.vertexCount()-1):
                            d[e.dest] = d[e.orig] + e.weight
                        else:
                            return None # one more relaxing reduced a distance -> negative cycle exists
                            
        return d

    def Dijkstra(self, v):
        """ Runs Dijkstra shortest-path algorithm from source v, returning the list of distances from v, or None if any edges are negative-weight. """
        if self.hasNegativeWeight():
            return None # Dijkstra may not run properly if any weights are negative

        # initialize the queue of vertices, the set s of "seen" or complete vertices, and the d list of current distances
        q = MinHeap()
        s = set()
        d = [float("inf")] * self.vertexCount()
        d[v] = 0

        # insert source vertex in q; could insert all if using decreaseKey, but a basic Heap implementation doesn't support that
        q.insert(KeyValuePair(d[v],v))

        # loop through each edge vertex once, in order of distance from v, relaxing all edges out of that vertex
        kvp = q.extractMin()
        while kvp:
            vert = kvp.value
            dist = kvp.key

            # only relax edges from a vertex if it has not already been visited
            if vert not in s:
                s.add(vert)
                # relax edges to vertices not already completed
                for e in self._adj[vert]:
                    if e.dest not in s:
                        if d[vert] + e.weight < d[e.dest]:
                            d[e.dest] = d[vert] + e.weight
                            q.insert(KeyValuePair(d[e.dest], e.dest)) # could also decrease key, but would not work with basic Heap implementation

            # retrieve the next-lowest distance vertex
            kvp = q.extractMin()
                
        return d    

    def DAGShortestPath(self, v):
        """ Runs a shortest path algorithm on a DAG, which visits each node in topological order, relaxing its edges. Returns None if the graph is not a DAG. """
        if self.isUndirected(): return None
        p, o, c = self.DFS()
        if c: return None

        d = [float("inf")] * self.vertexCount()
        d[v] = 0
        
        for i in o:
            edges = self._adj[i]
            for e in edges:
                if d[e.orig] + e.weight < d[e.dest]:
                   d[e.dest] = d[e.orig] + e.weight

        return d
            
    def allPairsShortestPaths(self):
        pass

    def Johnson(self):
        pass

    def FloydWarshall(self):
        pass

    def minimumSpanningTree(self):
        pass

    def Prim(self):
        pass

    def Kruskal(self):
        pass

def buildCompleteGraph(vertices, directed = True, weighted = True):
    """
    Builds a complete graph (an edge from every vertex to every other vertex), with a given vertex count and other optional attributes.
    Returns the Graph object.
    """
    graph = Graph(directed)
    for i in range(vertices):
        graph.addVertex()

    for i in range(vertices):
        if directed:
            start = 0
        else:
            start = i + 1
        
        for j in range(start,vertices):
            w = 1
            if weighted:
                w = random.randint(-5,10)
                while w == 0:
                    w = random.randint(-5,10)
            graph.addEdge(i,j,w)
        
    return graph

def buildRandomGraph(vertices, edges, directed = True, weighted = True):
    """
    Builds a random graph, with given vertex and edge counts, plus other optional attributes.
    Returns the Graph object.
    """
    graph = Graph(directed)
    for i in range(vertices):
        graph.addVertex()

    i = 0
    while i < edges:
        u, v = random.randint(0,vertices-1),random.randint(0,vertices-1)
        while v == u:
            v = random.randint(0,vertices-1)
            
        if not graph.findEdge(u, v):
            w = 1
            if weighted:
                w = random.randint(-1,10)
                while w == 0:
                    w = random.randint(-1,10)
            graph.addEdge(u,v,w)
            i += 1
        
    return graph

def buildTestGraph():
    """
    Builds an example graph for use in testing.
    """
    
    graph = Graph(True)
    for i in range(6):
        graph.addVertex()

    graph.addEdge(0,4,9)
    graph.addEdge(0,1,5)
    graph.addEdge(1,2,5)
    graph.addEdge(1,1,6)
    graph.addEdge(2,0,6)
    graph.addEdge(3,1,3)
    graph.addEdge(4,3,1)
    graph.addEdge(4,1,1)
    graph.addEdge(5,5,1)
    graph.addEdge(5,1,3)

    return graph

def main():
    kvp = KeyValuePair(1,0)
    print type(kvp) == KeyValuePair

    print "Complete graph:"
    graph = buildCompleteGraph(5, False, True)
    adj = graph.asAdjacency()
    for a in adj:
        for e in a:
            print e

    print "BFS on complete graph from node 0:"
    distance = graph.BFS(0)
    print distance

    print "DFS on complete graph:"
    p, o, c = graph.DFS()
    print "p =", p, "o =", o, "cycle =",c

    print "Adjacency matrix:"
    print graph.asMatrix()
    
    print "\nRandom graph:"
    graph = buildRandomGraph(6, 7, False)
    adj = graph.asAdjacency()
    for a in adj:
        for e in a:
            print e

    print "BFS on random graph from node 0:"
    distance = graph.BFS(0)
    print distance

    print "DFS on random graph:"
    p, o, c = graph.DFS()
    print "p =", p, "o =", o, "cycle =",c

    print "Topological sort:"
    print graph.topologicalSort()
    
    print "Adjacency matrix:"
    print graph.asMatrix()

    print "Bellman-Ford from node 0"
    d = graph.BellmanFord(0)
    print d

    print "Dijsktra from node 0"
    d = graph.Dijkstra(0)
    print d

    print "DAG shortest path from node 0"
    d = graph.DAGShortestPath(0)
    print d
    
    print "\nDesigned test graph:"
    graph = buildTestGraph()
    adj = graph.asAdjacency()
    for a in adj:
        for e in a:
            print e

    print "Bellman-Ford from node 0"
    d = graph.BellmanFord(0)
    print d

    print "Dijsktra from node 0"
    d = graph.Dijkstra(0)
    print d

    print "DAG shortest path from node 0"
    d = graph.DAGShortestPath(0)
    print d

    print "Eycle exists ", graph.hasCycle()
    print "Is DAG", graph.isDAG()
    print "Has negative weight edge", graph.hasNegativeWeight()
    print "Is directed", graph.isDirected()

    
if __name__ == "__main__":
    main()
