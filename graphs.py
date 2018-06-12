
import random
import sys
sys.path.insert(0, '../sorting')
sys.path.insert(0, '../data_structures')
from heap_sort import MinHeap
from hash_table import KeyValuePair


class Edge(object):
    def __init__(self, u, v, w = 1):
        self.orig = u
        self.dest = v
        self.weight = w

    def __str__(self):
        return "edge: " + str(self.orig) + "->" + str(self.dest) + ", wgt " + str(self.weight)

class Graph(object):
    def __init__(self, directed = False):
        self.adj = []
        self.directed = directed
        self.edge_count = 0
        self.neg_weights = False

    def addVertex(self):
        self.adj.append([])
        return len(self.adj)-1

    def addEdge(self, u, v, w = 1):
        if u >= len(self.adj) or v >= len(self.adj):
            # impossible to add edge - vertices don't exist
            return
        if not self.directed and u == v:
            # do not allow self loops in undirected graphs
            return
        if w == 0:
            # do not allow zero-weight edges
            return

        # track whether there are any negative weight edges in the graph
        if w < 0: self.neg_weights = True
        
        e = self.findEdge(u,v)
        if e:
            # overwrite existing edge
            e.weight = w
        else:
            # create new edge
            self.adj[u].append(Edge(u,v,w))
        
        if not self.directed:
            # create edge from v->u as well
            e = self.findEdge(u,v)
            if e:
                # overwrite existing edge
                e.weight = w
            else:
                # create new edge
                self.adj[v].append(Edge(v,u,w))

        self.edge_count += 1

    def findEdge(self, u, v):
        if u >= len(self.adj) or v >= len(self.adj):
            # vertices don't exist
            return None

        adj_u = self.adj[u]
        for e in adj_u:
            if e.dest == v:
                return e
        return None
        
    def asAdjacency(self):
        return self.adj

    def asMatrix(self):
        # create matrix of 0s
        m = [[]] * self.vertexCount()
        for i in range(len(m)):
            m[i] = [0] * self.vertexCount()

        # fill in weights for existing edges
        for i, edges in enumerate(self.adj):
            for e in edges:
                m[i][e.dest] = e.weight
        return m
        

    def DFS(self, node):
        if node >= len(self.adj):
            # can't be run, vertex doesn't exist
            return None
        
        seen = set()
        dist = [float("inf")] * self.vertexCount()
        par = [None] * self.vertexCount()
        fin = [None] * self.vertexCount()
        seen.add(node)
        dist[node] = 0
        par[node] = "Source"

        aux = [False,0] # aux[0] = cycle / back edge exists; aux[1] = edge count
        
        self.__DFSHelper(node,seen,dist,par,fin,aux)
        
        return (dist, par, fin, aux[0])

    def __DFSHelper(self, node, seen, dist, par, fin, aux):
        edges = self.adj[node]
        for e in edges:
            
            if e.dest not in seen:
                seen.add(e.dest)
                dist[e.dest] = dist[node]+1
                par[e.dest] = node
                self.__DFSHelper(e.dest,seen,dist,par,fin,aux)
            else:
                if not fin[e.dest]:
                    # cycle exists
                    aux[0] = True
        fin[node] = aux[1]
        aux[1] += 1
        
    def BFS(self, v):
        if v >= len(self.adj):
            # can't be run, vertex doesn't exist
            return None
        
        q = MinHeap()
        s = set()
        distance = [float("inf")] * self.vertexCount()

        q.insert(KeyValuePair(0,v))
        s.add(v)
        distance[v] = 0
        
        kvp = q.extractMin()
        while kvp:
            v = kvp.value
            d = kvp.key
            
            edges = self.adj[v]
            for e in edges:
                if e.dest not in s:
                    s.add(e.dest)
                    distance[e.dest] = d+1
                    q.insert(KeyValuePair(d+1,e.dest))
            kvp = q.extractMin()

        return distance

    def BellmanFord(self, s):
        pass

    def Dijkstra(self, s):
        pass

    def allPairsShortestPaths(self):
        pass

    def Johnson(self):
        pass

    def FloydWarshall(self):
        pass

    def vertexCount(self):
        return len(self.adj)

    def edgeCount(self):
        return self.edge_count

    def __str__(self):
        return str(self.adj)

def buildCompleteGraph(vertices, directed = True, weighted = True):
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
    graph = Graph(directed)
    for i in range(vertices):
        graph.addVertex()

    i = 0
    while i < edges:
        u, v = random.randint(0,vertices-1),random.randint(0,vertices-1)
        if not directed:
            while v == u:
                v = random.randint(0,vertices-1)
        if not graph.findEdge(u, v):
            w = 1
            if weighted:
                w = random.randint(-5,10)
                while w == 0:
                    w = random.randint(-5,10)
            graph.addEdge(u,v,w)
            i += 1
        
    return graph

def main():
    kvp = KeyValuePair(1,0)
    print type(kvp) == KeyValuePair

    print "Complete graph:"
    graph = buildCompleteGraph(5, True, True)
    adj = graph.asAdjacency()
    for a in adj:
        for e in a:
            print e

    print "BFS on complete graph from node 0:"
    distance = graph.BFS(0)
    print distance

    print "DFS on complete graph from node 0:"
    d, p, f, c = graph.DFS(0)
    print "d =", d, "p =", p, "f =", f, "cycle =",c

    print "Adjacency matrix:"
    print graph.asMatrix()
    
    print "\nRandom graph:"
    graph = buildRandomGraph(6, 6, True)
    adj = graph.asAdjacency()
    for a in adj:
        for e in a:
            print e

    print "BFS on random graph from node 0:"
    distance = graph.BFS(0)
    print distance

    print "DFS on random graph from node 0:"
    d, p, f, c = graph.DFS(0)
    print "d =", d, "p =", p, "f =", f, "cycle =",c

    print "Adjacency matrix:"
    print graph.asMatrix()
    
if __name__ == "__main__":
    main()