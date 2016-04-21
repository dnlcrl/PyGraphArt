#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import numpy.random as rnd
from scipy.spatial import cKDTree


class Graph():

    """
    Graph:
        nodes: a set of nodes i
        tree_nodes: a cKDTree repretenting the point list, used to query for
            the nearest neighbours in the edges creation prcess
        edges: a list of tuples (i,j) where i,j are in nodes
        costs: a list of cost where Ci is the cost of edges[i]
        FSs: forward stars, for every node i, FSs[i] is a list of j where
            (i,j) in edges
        FS_sosts: FSs forward stars costs, for every node i, FSs[i] is a
            list of costs Cj where (i,j) in edges
        BSs: backward stars, for every node j, BSs[j] is a list of i where
            (i,j) in edges
        gamma: nodes-edges incidence matrix, a len(nodes)Xlen(edges) matrix
            where gamma[e,n] = Ce + if n is tail(e), + Ce + if n is head(e),
            0 else
        qu: adiacence matrix, a len(nodes)Xlen(nodes) matrix where qu[i,j] =
            Cij if (i,j) in nodes, -Cji if (j,i) in nodes, 0 else
        capacities: capacity vector, for each edge e in edges, capacities[e]
            represent the capacity value of e
    """
    nodes = np.array([])
    edges = np.array([])
    FSs = np.array([])
    FS_costs = np.array([])
    BSs = np.array([])
    gamma = np.array([])
    qu = np.array([])
    capacities = np.array([])
    flows = np.array([])

    def __init__(self, oriented=True, points=None, edges=None, rand=True, n=2000, max_neighbours=20, acyclic=False, negative_costs=False, maxflux=False):
        '''
        Initialize the graph's data structures:
        param:
        oriented: True if the graph is oriented
        points: a list of nodes coordinates
        edges: a list of edges i,j where i,j are in range(len(points))
        rand: if True, generate random edges searching the rand(2,
            max_neighbours) nearest neighbours of each node
        n: number of nodes, used if points is None to generate 5*n random nodes
            and select the n nodes nearest to coords 0,0
        max_neighbours: max number of neighbours to search for in the random
            edges creation process
        acyclic: True if the graph must be acyclic
        negative_costs: if acyclic is True, negative and negative_costs
            indicates if negative costs edges could be added
        maxflux: True if a flux network has to be build (using edge capacities)

        '''
        self.max_neighbours = max_neighbours if max_neighbours < n else n - 1
        self.oriented = oriented
        self.acyclic = acyclic
        self.negative_costs = negative_costs
        self.maxflux = maxflux
        if points is None:
            self.n = n  # number of nodes
            self._init_random_nodes()
        else:
            self.nodes = points
            self.n = len(points)
        self._init_nodes()
        self._init_structs()
        if not edges and not rand:
            print 'you must provide an edge list or rand=True'
            return
        if rand:
            self._init_edges()
            self._init_edges_list()
            self._init_gamma()
            self._init_qu()
        if edges:
            self.add_edges(edges)
            if self.is_cyclic(self.FSs) and not self.cyclic:
                print 'error, edges provides contains cycles and cyclic is False'
                return

        # FIXME here in some case (using italy's countries coords)
        # the edges generation process generate unreachable nodes, thus BS[i]
        # is [] and this could raise some error
        if self.acyclic:
            for i in range(self.n):
                if len(self.BSs[i]) == 0:
                    self.add_edges([[head, i] for head in self.FSs[i]])

        if maxflux:
            if self.oriented is False:
                print 'cannot apply maxflux on an unoriented graph'
            if self.negative_costs is True:
                print 'costs will still be positive'
            self._init_capacities()
            self._init_flux()

        print 'number of nodes: ', self.n
        print 'number of arches: ', len(self.edges)

    def _init_random_nodes(self):
        '''
        Initialize uniformly random nodes coordinates
        '''
        s = self.n*5  # number of random samples
        # generate random points
        self.nodes = np.array(
            1.*rnd.uniform(0, 1, size=(s, 2)), dtype=np.float32)

    def _init_nodes(self):
        '''
        Initialize nodes structures
        '''
        center = np.array([[0.5, 0.5]],  dtype=np.float32)
        nodes_tree = cKDTree(self.nodes)
        # keep only n points nearest to 0, 0
        nodes_indices = nodes_tree.query(center, k=self.n)[1][0]
        self.nodes = self.nodes[nodes_indices]
        self.nodes_tree = cKDTree(self.nodes)

    def add_nodes(self, nodes):
        for n in nodes:
            if n not in self.nodes:
                self.n += 1
                self.nodes = np.append(self.nodes, n).reshape(self.n, 2)
                self.FSs = np.append(self.FSs, np.empty(1, dtype=np.object))
                self.FS_costs = np.append(
                    self.FS_costs, np.empty(1, dtype=np.object))
                self.FSs[-1] = np.array([])
                self.FS_costs[-1] = np.array([])
                # self.BSs = np.append(self.BSs, np.array([], dtype=np.int))
        self._init_BSs()

        self._init_gamma()
        self._init_qu()

    def add_edges(self, edges):
        '''
        Adds a list of (i,j) where i,j are in range(n) and reinitialize the
        structures
        '''
        for e in edges:
            if e[1] not in self.FSs[e[0]]:
                self.edges = np.append(self.edges, None)
                self.edges[-1] = (e[0], e[1])
                self.FSs[e[0]] = np.append(self.FSs[e[0]], None)
                self.FSs[e[0]][-1] = e[1]
                self.FS_costs[e[0]] = np.append(self.FS_costs[e[0]], None)
                c = ((self.nodes[e[1]] - self.nodes[e[0]])**2)
                c = np.sqrt(c.sum(axis=-1))
                self.FS_costs[e[0]][-1] = c
        self._init_BSs()
        self._init_gamma()
        self._init_qu()

    def remove_edges(self, edges):
        '''
        Adds a list of (i,j) where i,j are in range(n) and reinitialize the
        structures
        '''
        for e in edges:
            if e[1] in self.FSs[e[0]]:
                self.edges = np.delete(self.edges, np.where(self.edges == e))
                self.FS_costs[e[0]] = np.delete(
                    self.FS_costs[e[0]], np.where(self.FSs[e[0]] == e[1]))
                self.FSs[e[0]] = np.delete(
                    self.FSs[e[0]], np.where(self.FSs[e[0]] == e[1]))
        self._init_BSs()
        self._init_gamma()
        self._init_qu()

    def _init_structs(self):
        '''
        Initialize the data structures
        '''
        self.edges = np.array([], dtype=object)
        self.FSs = np.empty((self.n), dtype=np.object)
        self.BSs = np.empty((self.n), dtype=np.object)
        self.FS_costs = np.empty((self.n), dtype=np.object)
        for i in range(self.n):
            self.FSs[i] = np.array([], dtype=np.int)
            self.BSs[i] = np.array([], dtype=np.int)
            self.FS_costs[i] = np.array([], dtype=np.float32)

    def _init_edges(self):
        '''
        Initialize random edges filling the FSs struct
        '''
        if self.acyclic:
            self._init_acyclic_edges()
            return
        for i in range(len(self.nodes)):
            # rnd.randint(, self.max_neighbours+1)
            neighbours = self.max_neighbours
            query_res = self.nodes_tree.query(
                self.nodes[i], k=neighbours+1)

            self.FSs[i] = query_res[1][1:]
            for node_index in self.FSs[i]:
                self.BSs[node_index] = np.append(self.BSs[node_index], [i])

            self.FS_costs[i] = query_res[0][1:]

    def _init_acyclic_edges(self):
        '''
        Initialize random edges filling the FSs struct
        '''
        print '''Warning: Generating acyclic edges is NP-hard, expect long times, very long for more than 500 nodes'''
        big_s = set([])
        for i in range(len(self.nodes)):
            # if len(self.BSs[i]) == 0:
            big_s.add(i)
            neighbours = rnd.randint(3, self.max_neighbours+1)
            query_res = self.nodes_tree.query(
                self.nodes[i], k=neighbours+1)
            hole = len(self.BSs[i]) != 0
            for tail, cost in zip(query_res[1][1:], query_res[0][1:]):
                if not hole:
                    self.FSs[i] = np.append(self.FSs[i], tail)
                    self.FS_costs[i] = np.append(self.FS_costs[i], cost)
                    self.BSs[tail] = np.append(self.BSs[tail], i)
                else:
                    temp_fs = self.FSs.copy()
                    temp_fs[i] = np.append(temp_fs[i], tail)
                    if not self.is_cyclic(temp_fs):
                        self.FSs[i] = temp_fs[i].copy()
                        self.FS_costs[i] = np.append(self.FS_costs[i], cost)
                        self.BSs[tail] = np.append(self.BSs[tail], i)
                    # elif self.negative_costs:
                    #     self.BSs[i] = np.append(self.BSs[i], tail)
                    #     self.FS_costs[tail] = np.append(
                    #         self.FS_costs[tail], -cost)
                    #     self.FSs[tail] = np.append(self.FSs[tail], i)

                #big_s = big_s.union(set([i for i in self.FSs[i]]))

    def _init_edges_list(self):
        '''
        Initialize the edge list and cost list
        '''
        l = sum([len(x) for x in self.FSs])
        self.edges = np.empty((l), dtype=object)
        self.costs = np.empty((l), dtype=object)
        index = 0
        for i in range(len(self.FSs)):
            for n_index in range(len(self.FSs[i])):
                self.edges[index] = [i, self.FSs[i][n_index]]
                self.costs[index] = self.FS_costs[i][n_index]
                index += 1

    def _init_BSs(self):
        '''
        Initialize BSs struct
        '''
        self.BSs = np.empty((self.n), dtype=np.object)

        for i in range(len(self.FSs)):
            for node_index in self.FSs[i]:
                if self.BSs[node_index] is None:
                    self.BSs[node_index] = np.array([], dtype=np.int)
                self.BSs[node_index] = np.append(self.BSs[node_index], [i])

    def _init_gamma(self):
        '''
        Initialize Gamma
        '''
        l = sum([len(x) for x in self.FSs])
        h = self.n
        self.gamma = np.zeros((l, h), dtype=np.float32)
        gamma_index = 0
        for x in range(self.n):
            for fs_index in range(len(self.FSs[x])):
                if self._init_gamma_edge(gamma_index, x, fs_index):
                    gamma_index += 1
        self.gamma = self.gamma[~np.all(self.gamma == 0, axis=1)]

    def _init_gamma_edge(self, gamma_index, x, y):
        '''
        Try to add a column to Gamma, return True if success
        '''
        column = np.zeros(self.n, dtype=np.float32)
        direction = 1
        if self.oriented:
            direction = -1
        column[x] = direction * self.FS_costs[x][y]
        column[self.FSs[x][y]] = self.FS_costs[x][y]
        self.gamma[gamma_index] = column
        return True

    def _init_qu(self):
        '''
        Initialize Qu
        '''
        l = h = len(self.nodes)
        self.qu = np.zeros((l, h), dtype=np.float32)
        direction = 1
        if self.oriented:
            direction = -1
        for x in self.gamma:
            nonzero = x.nonzero()[0]
            if x[nonzero[0]] > 0:
                invert = 1
            else:
                invert = -1
            self.qu[nonzero[0], nonzero[1]] = invert*x[nonzero[0]]
            self.qu[nonzero[1], nonzero[0]] = invert*direction*x[nonzero[0]]

    def _init_capacities(self):
        '''
        Inizialize capacity vector
        '''
        self.capacities = np.random.randint(1, 11, len(self.edges))

    def _init_flux(self):
        '''
        Inizialize flux vector
        '''
        self.flows = np.zeros((len(self.edges), 0), dtype=np.int)

    def edges_to_fs(self, edges):
        """
        Return the FS list for a given list of edges
        """
        return [edges[edges[:, 0] == i][:, 1] for i in range(1 + np.max(edges.flatten()))]

    def is_cyclic(self, fs):
        """
        Return True if the UNdirected graph fs has a cycle.
        fs must be a list of FSs
        """
        path = set()
        visited = set()

        def visit(vertex):
            if vertex in visited:
                return False
            visited.add(vertex)
            path.add(vertex)
            for neighbour in fs[vertex]:
                if neighbour in path or visit(neighbour):
                    return True
            path.remove(vertex)
            return False
        return any(visit(v) for v in range(len(fs)))

    def delta_plus(self, nodes):
        '''
        Returns the list of edges forwarding from a set of nodes
        '''
        bool_indices_head = np.array([x[0] in nodes for x in self.edges])
        bool_indices_tail = np.array([x[1] not in nodes for x in self.edges])
        bool_indices_edges = np.bitwise_and(
            bool_indices_head, bool_indices_tail)
        return self.edges[bool_indices_edges]

    def delta_minus(self, nodes):
        '''
        Returns the list of edges backwarding from a set of nodes
        '''
        bool_indices_head = np.array([x[0] not in nodes for x in self.edges])
        bool_indices_tail = np.array([x[1] in nodes for x in self.edges])
        bool_indices_edges = np.bitwise_and(
            bool_indices_head, bool_indices_tail)
        return self.edges[bool_indices_edges]

    def coords_edge(self, edge):
        '''
        Returns coordinates head and tail points for an edge
        '''
        return (self.nodes[edge[0]], self.nodes[edge[1]])

    def coords_edges(self, edges):
        '''
        Returns a list of coordinates head and tail points for all edge in edges
        '''
        res = np.empty((len(edges)), dtype=object)
        for r, e in zip(range(len(edges)), edges):
            if e[0] is None:
                e[0] = 0
            res[r] = self.coords_edge(e)
            if len(res[r][0]) != 2:
                print 'there is an error with the edges'
                import pdb
                pdb.set_trace()

        # v = np.vectorize(self.coords_edge, otypes=[np.object])
        # res = v(edges)
        return res

    def BFS(self, start, fs=None):
        '''
        Returns the BFS tree for the graph starting from start
        '''
        to_be_processed = np.array([start], dtype=np.int)
        known = np.array([], dtype=np.int)
        tree = np.array([], dtype=object)
        if fs is None:
            fs = self.FSs
        while len(to_be_processed) > 0:
            # pop
            current_node = to_be_processed[-1]
            to_be_processed = np.delete(to_be_processed, -1)

            for node in fs[current_node]:
                if node not in known:
                    known = np.append(known, node)
                    tree = np.append(tree, None)
                    tree[-1] = (current_node, node)
                    # push
                    to_be_processed = np.insert(to_be_processed, 0, node)

        return tree

    def BFS_path(self, start, end, fs=None):
        '''
        Returns the BFS part for the graph from start to end if one
        '''
        to_be_processed = np.array([start], dtype=np.int)
        known = np.array([], dtype=np.int)
        prev = np.array([None for i in range(self.n)])
        if fs is None:
            fs = self.FSs
        while len(to_be_processed) > 0:
            # pop
            current_node = to_be_processed[-1]
            to_be_processed = np.delete(to_be_processed, -1)

            for node in fs[current_node]:
                if node not in known:
                    known = np.append(known, node)
                    prev[node] = current_node
                    # push
                    to_be_processed = np.insert(to_be_processed, 0, node)
                    if node == end:
                        to_be_processed = []
                        break
                        # return prev
        edges = np.array([], dtype=object)
        e = end
        while True:
            if e == start:
                return edges
            if prev[e] is None:
                return np.array([])
            edges = np.append(None, edges)
            edges[0] = [prev[e], e]
            e = prev[e]

    def DFS(self, start, fs=None):
        '''
        Returns the DFS tree for the graph starting from start
        '''
        to_be_processed = np.array([start], dtype=np.int)
        known = np.array([], dtype=np.int)
        tree = np.array([], dtype=object)
        if fs is None:
            fs = self.FSs
        while len(to_be_processed) > 0:
            # pop
            current_node = to_be_processed[0]
            to_be_processed = np.delete(to_be_processed, 0)

            for node in fs[current_node]:
                if node not in known:
                    known = np.append(known, node)
                    tree = np.append(tree, None)
                    tree[-1] = (current_node, node)
                    # push
                    to_be_processed = np.insert(to_be_processed, 0, node)

        return tree

    def prim(self):
        '''
        Returns Prim's minimum spanninng tree
        '''
        big_f = set([])
        costs = np.empty((self.n), dtype=object)
        costs[:] = np.max(self.costs) + 1
        big_e = np.empty((self.n), dtype=object)
        big_q = set(range(self.n))
        tree_edges = np.array([], dtype=object)
        while len(big_q) > 0:
            v = np.argmin(costs)
            big_q.remove(v)
            costs[v] = np.Infinity
            big_f.add(v)
            if big_e[v] is not None:
                tree_edges = np.append(tree_edges, None)
                tree_edges[-1] = (big_e[v], v)

            for i, w in zip(range(len(self.FSs[v])), self.FSs[v]):
                if w in big_q and self.FS_costs[v][i] < costs[w]:
                    costs[w] = self.FS_costs[v][i]
                    big_e[w] = v
        return tree_edges

    def kruskal(self, maximum=False):
        '''
        Returns Prim's minimum spanninng tree
        '''
        sets = set([frozenset([x]) for x in range(self.n)])
        tree_edges = np.array([], dtype=object)
        sorted_index = np.argsort(self.costs)
        if maximum:
            sorted_index = sorted_index[::-1]
        sorted_edges = self.edges[sorted_index]
        while len(tree_edges) < self.n - 1 and len(sorted_edges) > 0:
            edge = sorted_edges[0]
            sorted_edges = np.delete(sorted_edges, 0)
            u = frozenset([edge[0]])
            v = frozenset([edge[1]])

            find_set_u = [s for s in sets if not u.isdisjoint(s)][0]
            find_set_v = [s for s in sets if not v.isdisjoint(s)][0]
            if not find_set_u == find_set_v:
                sets.remove(find_set_u)
                sets.remove(find_set_v)
                sets.add(find_set_u.union(find_set_v))
                tree_edges = np.append(tree_edges, None)
                tree_edges[-1] = edge

        return tree_edges, sets

    def connect_graphs(self, sets_orig, edges_orig):
        '''
        Returns the edges needed to connect unconnected graphs (sets of nodes)
        given a set of sets of nodes, select the master_graph (the biggest) one
        and search the shortest edges to connect the other sets of nodes
        '''
        master_graph = max(sets_orig, key=len)
        sets = sets_orig.copy()
        edges = np.array([], dtype=object)
        sets.remove(master_graph)
        master_tree = cKDTree(self.nodes[list(master_graph)])
        for s in sets:
            x = np.array(list(s))
            nearests = np.array([master_tree.query(self.nodes[v]) for v in x])
            tails = nearests[
                nearests[:, 0].argsort()][:, 1][:self.max_neighbours]
            heads = x[nearests[:, 0].argsort()][:self.max_neighbours]
            for head, tail in zip(heads, tails):
                edges = np.append(edges, None)
                edges[-1] = (head, tail)
                edges = np.append(edges, None)
                edges[-1] = (tail, head)
        return edges

    def dijkstra(self, source):
        '''
        Returns Dijkstra's shortest paths from source to all other nodes,
        if the (directed) graph doesn't contains negative costs
        '''
        if self.oriented is False:
            print 'cannot apply dijstra, graph is not oriented'
            return
        if len(self.costs[self.costs < 0]) > 0:
            print 'cannot apply dijstra, graph contains negative cost edges'
            return

        big_q = np.arange(self.n)
        dist = np.empty((self.n), dtype=np.float32)
        dist[:] = np.finfo(np.float32()).max
        prev = np.empty((self.n), dtype=np.object)
        prev[:] = None
        prev[source] = source
        dist[source] = 0
        while len(big_q) > 0:
            u_i = np.argmin(dist[big_q])
            u = big_q[u_i]
            big_q = np.delete(big_q, u_i)
            for v in [v for v in self.FSs[u] if v in big_q]:
                alt = dist[u] + \
                    self.FS_costs[u][np.where(self.FSs[u] == v)[0][0]]
                if alt < dist[v]:
                    dist[v] = alt
                    prev[v] = u
        edges, prev = self.fix_unheaded_edges(prev)
        return edges

    def fix_unheaded_edges(self, prev):
        '''
        Add prev to all nodes n where prev[n] is None, after dijstra some of
        the nodes could still be unreachable, thus edge[0] is None, to prevent
        errors this ugly function makes sure that all edges are reachable from
        source after dijstra, hopefully at the minimum cost.
        '''
        edges = np.array([], dtype=object)
        prev_presence = np.array([False if x is None else True for x in prev])
        while False in prev_presence:
            for v in range(len(prev)):

                if prev[v] is None:
                    for fs in self.FSs[v][np.argsort(self.FS_costs[v])]:
                        if prev_presence[fs]:
                            prev[v] = fs
                            break
                if prev[v] is None:
                    for bs in self.BSs[v]:
                        if prev_presence[bs]:
                            prev[v] = bs
                            break
            prev_presence = np.array(
                [False if x is None else True for x in prev])

        for v in range(len(prev)):
            edges = np.append(edges, None)
            edges[-1] = [prev[v], v]
        return edges, prev

    def topological_sort(self):
        '''
        Returns a list topological sorted nodes
        '''
        if self.is_cyclic(self.FSs):
            print 'cannot apply labels, graph contains cycles'
            return
        big_l = []  # Empty list that will contain the sorted elements
        # Set of all nodes with no incoming edges
        big_s = set([0])
        bs_copy = self.BSs.copy()
        while len(big_s) > 0:
            n = big_s.pop()
            big_l.append(n)
            for m in self.FSs[n]:
                bs_copy[m] = np.delete(bs_copy[m], np.where(bs_copy[m] == n))
                # bs_copy[m].remove(n)
                if len(bs_copy[m]) == 0:
                    big_s.add(int(m))
        return big_l

    def labels(self, source):
        '''
        Returns Labels-algorithm's shortest paths from source to all other
        nodes, if the (directed) graph doesn't contains cycles
        '''
        if self.oriented is False:
            print 'cannot apply labels, graph is not oriented'
            return
        arg_sort = np.array(self.topological_sort())
        sorted_BSs = self.BSs[arg_sort]

        for i in range(len(sorted_BSs)):
            bs_i = []
            for e in sorted_BSs[i]:
                try:
                    bs_i.append(np.where(arg_sort == e)[0][0])
                except:
                    print 'error'
            sorted_BSs[i] = np.array(bs_i)

        pred = np.empty((self.n), dtype=np.int)
        pred[source] = 0
        edges = np.array([], dtype=object)

        f = []
        f.append(source)
        for t in range(1, len(arg_sort)):
            prev_costs = []
            for p in sorted_BSs[t]:
                prev_costs.append(f[p] + self.qu[arg_sort[p]][arg_sort[t]])

            f.append(min(prev_costs))
            pred[t] = sorted_BSs[t][np.argmin(prev_costs)]
            edges = np.append(edges, None)
            edges[-1] = [arg_sort[pred[t]], arg_sort[t]]
        return edges

    def bellman_ford(self, source):
        '''
        Returns Labels-algorithm's shortest paths from source to all other
        nodes, if the (directed) graph doesn't contains cycles
        '''
        if self.oriented is False:
            print 'cannot apply bellman_ford, graph is not oriented'
            return
        dist = np.array([np.Infinity for x in range(self.n)], dtype=np.float32)
        pred = np.empty((self.n), dtype=np.int)
        pred[source] = source
        dist[source] = 0

        for i in np.arange(1, self.n):
            for e in range(len(self.edges)):
                if dist[self.edges[e][0]] + self.costs[e] < dist[self.edges[e][1]]:
                    dist[self.edges[e][1]] = dist[
                        self.edges[e][0]] + self.costs[e]
                    pred[self.edges[e][1]] = self.edges[e][0]

        for e in range(len(self.edges)):
            if dist[self.edges[e][1]] > dist[self.edges[e][0]] + self.costs[e]:
                print 'Error, Graph contains a negative-weight cycle'
                break

        edges = np.array([], dtype=object)
        for v in range(len(pred)):
            edges = np.append(edges, None)
            edges[-1] = [pred[v], v]

        return edges  # , prev, dist

    def floyd_warshall(self, source):
        '''
        Returns floyd_warshall's shortest paths from source to all other
        nodes, if the (directed) graph doesn't contains negative cycles
        '''
        print '''warning! apply this algorithm only if constricted, it takes\\
        O(n^3)!'''
        print 'O(n^3) = O(', self.n**3, ')'
        dist = np.empty((self.n, self.n), dtype=np.float32)
        pred = np.zeros((self.n, self.n), dtype=np.int)
        dist.fill(np.Infinity)
        for v in range(self.n):
            dist[v][v] = .0
        for e in range(len(self.edges)):
            u = self.edges[e][0]
            v = self.edges[e][1]
            dist[u][v] = self.costs[e]
            pred[u][v] = v
        for h in range(1, self.n):
            for i in range(1, self.n):
                for j in range(self.n):
                    if dist[i][h] + dist[h][j] < dist[i][j]:
                        dist[i][j] = dist[i][h] + dist[h][j]
                        pred[i][j] = pred[h][j]
            for i in range(1, self.n):
                if dist[i][i] < 0:
                    print 'Error! found negative cycle, thus the problem is inferiorly unlinmited'
                    return
        edges = np.array([], dtype=object)
        for v in range(len(pred)):
            edges = np.append(edges, None)
            edges[-1] = [pred[source][v], v]

        return edges  # , prev, dist

    def ford_fulkerson(self):
        '''
        Returns ford_fulkerson's max flux in a flux graph
        '''
        if not self.maxflux:
            print 'cannot apply for fulkerson, graph is not a flow graph'
        print '''warning! this algorithm takes O(m^2n)!'''
        print 'O(m^2n) = O(', self.n*len(self.edges)**2, ')'
        # fixgraph self
        self.remove_edges([[0, i]for i in self.FSs[0]])
        self.remove_edges([[i, self.n-1] for i in self.FSs[-1]])

        s_node = np.array([-.85, -.85])
        s = self.n-2
        t_node = np.array([0.85, 0.85])
        t = self.n-1

        self.add_nodes([s_node, t_node])

        for v in range(len(self.FSs)-2):
            if self.FSs[v] == []:
                self.add_edges([v, t])
                self.capacities = np.append(
                    self.capacities, sum(self.capacities[np.where(self.edges[:][1] == v)]))
                self.flux = np.append(self.flux, 0)

#                c_s.append(self.capacities[np.where(self.edges == [v])])
        for v in range(len(self.BSs)-2):
            if self.BSs[v] == []:
                self.add_edges([s, v])
                self.capacities = np.append(
                    self.capacities, sum(self.capacities[np.where(self.edges[:][0] == v)]))
                self.flux = np.append(self.flux, 0)
        max_flux = np.zeros(len(self.edges), dtype=np.int)
        v = 0
        optimal = False
        directed_edges = self.edges.copy()
        directed_capacities = self.capacities.copy()
        inverted_edges = self.edges.copy()
        inverted_edges.fill(None)
        inverted_capacities = self.capacities.copy()
        inverted_capacities.fill(0)
        fs = self.FSs.copy()
        while not optimal:
            for e in range(len(self.edges)):
                if max_flux[e] != self.capacities[e]:

                    fs[self.edges[e][0]] = np.append(
                        fs[self.edges[e][0]], self.edges[e][1])
                    directed_edges[e] = self.edges[e]
                    directed_capacities[e] = self.capacities[e] - max_flux[e]
                else:
                    directed_edges[e] = None
                    directed_capacities[e] = 0
                if max_flux[e] != 0:
                    fs[self.edges[e][1]] = np.append(
                        fs[self.edges[e][1]], self.edges[e][0])
                    inverted_edges[e] = self.edges[e]
                    inverted_capacities[e] = max_flux[e]
                else:
                    inverted_edges[e] = None
                    inverted_capacities[e] = 0
            big_p = self.BFS_path(s, t, fs)
            if len(big_p) == 0:
                optimal = True
                break

            u_list = []

            for e in range(len(directed_edges)):
                if directed_edges[e]:
                    for p in big_p:
                        if directed_edges[e] == p:
                            u_list.append(directed_capacities[e])
                if inverted_edges[e]:
                    for p in big_p:
                        if inverted_edges[e] == p:
                            u_list.append(inverted_capacities[e])
            big_fi = min(u_list)
            v = v + big_fi

            for e in range(len(directed_edges)):
                if len(np.where(directed_edges[e] == big_p)[0]) > 1:
                    max_flux[e] = max_flux[e] + big_fi
            for e in range(len(inverted_edges)):
                if len(np.where(directed_edges[e] == big_p)[0]) > 1:
                    max_flux[e] = max_flux[e] - big_fi
        return self.edges, max_flux

    def tsp_hamilton_cycle(self, start):
        sorted_edges = np.argsort(self.costs)
        constraints = self.costs.copy()
        constraints.fill(None)
        self.upper_bound = 0
        self.final_x = None
        edges = np.array([i for i in self.edges])
        opt_c = np.Infinity
        if self.oriented:
            print 'tsp only on unoriented graph'
            return
        problems = np.array([], dtype=np.object)
        problems = np.append(problems, None)
        problems[0] = constraints
        while len(problems) != 0:
            _ = problems[0]
            np.delete(problems, 0)
            one_tree = np.array([], dtype=np.int)
            chosen = np.empty(len(edges), dtype=np.bool)
            chosen.fill(False)
            for i in sorted_edges:
                if constraints[i] is not None:
                    if constraints[i] == 1:
                        chosen[i] = True
                        if len(one_tree) > 0:
                            one_tree = np.vstack((one_tree, [edges[i]]))
                        else:
                            one_tree = [edges[i]]

                        # visited_nodes = np.unique(
                        #   np.append(visited_nodes, edges[i]))

            for i in sorted_edges:
                todo1 = todo2 = False
                if len(one_tree) > 0:
                    bincount = np.bincount(np.hstack(one_tree))
                    if len(bincount) > edges[i][0]:
                        if bincount[edges[i][0]] <= 1:
                            todo1 = True
                    else:
                        todo1 = True
                    if len(bincount) > edges[i][1]:
                        if np.bincount(np.hstack(one_tree))[edges[i][1]] <= 1:
                            todo2 = True
                    else:
                        todo2 = True
                else:
                    todo1 = todo2 = True

                cycle = False
                if len(one_tree) > 0:
                    cycle = self.is_cyclic(
                        self.edges_to_fs(np.vstack((one_tree, [edges[i]]))))

                if (constraints[i] is None) and todo1 and todo2 and not cycle:
                    chosen[i] = True
                    if len(one_tree) > 0:
                        one_tree = np.vstack((one_tree, [edges[i]]))
                    else:
                        one_tree = [edges[i]]
                    # visited_nodes = np.unique(
                    #    np.append(visited_nodes, edges[i]))
                    if len(edges[chosen]) == self.n-1:
                        break
            for i in sorted_edges[~chosen]:
                if not np.any(one_tree == i):
                    if len(one_tree) > 0:
                        one_tree = np.vstack((one_tree, [edges[i]]))
                    else:
                        one_tree = [edges[i]]
                    chosen[i] = True
                    break
            lowerbound = np.sum(self.costs[chosen])
            if lowerbound > opt_c:
                continue
            #

            if np.all(np.bincount(np.hstack(one_tree)) == 2):
                # np.max(np.bincount(np.hstack(one_tree))) == 2:

                if lowerbound < opt_c:
                    opt_c = lowerbound
                    continue

            node = np.argmax(np.bincount(np.hstack(one_tree)))
            eds = self.edges[np.bitwise_or(constraints == 1, constraints == 0)]
            if len(eds) > 0:
                if np.sum(np.hstack(eds) == node) == 2:
                    continue
            min_cost = np.Infinity
            for i in sorted_edges:
                if node in edges[i] and self.costs[i] < min_cost:
                    branch_i = i

            constraints[branch_i] = 0
            cycle = self.is_cyclic(self.edges_to_fs(edges[constraints == 1]))
            if not cycle:
                problems = np.append(problems, None)
                problems[-1] = constraints

            constraints[branch_i] = 1
            cycle = self.is_cyclic(self.edges_to_fs(edges[constraints == 1]))
            if not cycle:
                problems = np.append(problems, None)
                problems[-1] = constraints

        return edges[constraints == 1]

    def tsp(self, start):
        res = []
        nodes = np.array([i for i in range(self.n)])
        visited = np.empty(self.n, dtype=np.bool)
        visited.fill(False)
        visited[start] = True
        node = start

        while False in visited:
            tree = cKDTree(self.nodes[~visited])
            nearest = tree.query(self.nodes[node], k=1)[1]
            t = nodes[~visited][nearest]
            res.append([node, t])
            visited[t] = True
            node = t
        return res + [node, start]
