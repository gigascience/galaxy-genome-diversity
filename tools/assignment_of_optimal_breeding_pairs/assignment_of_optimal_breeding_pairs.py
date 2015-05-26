#!/usr/bin/env python2.6

import sys
import munkres
import random

class Vertex(object):
    def __init__(self, name):
        self.name = name
        self.neighbors = {}
        self.color = 0
        self.explored = False

    def add_neighbor(self, neighbor, weight=0.0):
        if neighbor in self.neighbors:
            if self.neighbors[neighbor] != weight:
                die('multiple edges not supported')
        else:
            self.neighbors[neighbor] = weight

class Graph(object):
    def __init__(self):
        self.vertex_list = {}
        self.vertices = 0
        self.max_weight = 0.0

    def add_vertex(self, name):
        if name not in self.vertex_list:
            self.vertex_list[name] = Vertex(name)
            self.vertices += 1
        return self.vertex_list[name]

    def add_edge(self, name1, name2, weight):
        vertex1 = self.add_vertex(name1)
        vertex2 = self.add_vertex(name2)
        vertex1.add_neighbor(vertex2, weight)
        vertex2.add_neighbor(vertex1, weight)
        self.max_weight = max(self.max_weight, weight)

    def from_edge_file(self, filename):
        fh = try_open(filename)
        line_number = 0
        for line in fh:
            line_number += 1
            line = line.rstrip('\r\n')
            elems = line.split()
            if len(elems) < 3:
                die('too few columns on line {0} of {1}:\n{2}'.format(line_number, filename, line))
            name1 = elems[0]
            name2 = elems[1]
            weight = float_value(elems[2])
            if weight is None:
                die('invalid weight on line {0} of {1}:\n{2}'.format(line_number, filename, line))
            self.add_edge(name1, name2, weight)
        fh.close()

    def bipartite_partition(self):
        vertices_left = self.vertex_list.values()

        while vertices_left:
            fifo = [vertices_left[0]]
            while fifo:
                vertex = fifo.pop(0)
                if not vertex.explored:
                    vertex.explored = True
                    vertices_left.remove(vertex)

                    if vertex.color == 0:
                        vertex.color = 1
                        neighbor_color = 2
                    elif vertex.color == 1:
                        neighbor_color = 2
                    elif vertex.color == 2:
                        neighbor_color = 1

                    for neighbor in vertex.neighbors:
                        if neighbor.color == 0:
                            neighbor.color = neighbor_color
                        elif neighbor.color != neighbor_color:
                            return None, None
                        fifo.append(neighbor)

        c1 = []
        c2 = []

        for vertex in self.vertex_list.values():
            if vertex.color == 1:
                c1.append(vertex)
            elif vertex.color == 2:
                c2.append(vertex)

        return c1, c2

def try_open(*args):
    try:
        return open(*args)
    except IOError:
        die('Failed opening file: {0}'.format(args[0]))

def float_value(token):
    try:
        return float(token)
    except ValueError:
        return None

def die(message):
    print >> sys.stderr, message
    sys.exit(1)

def main(input, randomizations, output):
    graph = Graph()
    graph.from_edge_file(input)
    c1, c2 = graph.bipartite_partition()

    if c1 is None:
        die('Graph is not bipartite')

    if len(c1) + len(c2) != graph.vertices:
        die('Bipartite partition failed: {0} + {1} != {2}'.format(len(c1), len(c2), graph.vertices))

    with open(output, 'w') as ofh:
        a1 = optimal_assignment(c1, c2, graph.max_weight)
        optimal_total_weight = 0.0
        for a in a1:
            optimal_total_weight += a[0].neighbors[a[1]]

        print >> ofh, 'optimal average {0:.3f}'.format(optimal_total_weight / len(a1))

        if randomizations > 0:
            random_total_count = 0
            random_total_weight = 0.0
            for i in range(randomizations):
                a2 = random_assignment(c1, c2)
                random_total_count += len(a2)
                for a in a2:
                    random_total_weight += a[0].neighbors[a[1]]
            print >> ofh, 'random average {0:.3f}'.format(random_total_weight / random_total_count)


        for a in a1:
            print >> ofh, '\t'.join([a[0].name, a[1].name])

def optimal_assignment(c1, c2, max_weight):
    matrix = []
    assignment = []

    for v1 in c1:
        row = []
        for v2 in c2:
            row.append(max_weight + 1.0 - v1.neighbors[v2])
        matrix.append(row)

    m = munkres.Munkres()
    indexes = m.compute(matrix)
    for row, column in indexes:
        assignment.append([c1[row], c2[column]])

    return assignment

def random_assignment(c1, c2):
    assignment = []

    ## note, this assumes that graph is complete bipartite
    ## this needs to be fixed
    c1_len = len(c1)
    c2_len = len(c2)
    idx_list = list(range(max(c1_len, c2_len)))
    random.shuffle(idx_list)

    if c1_len <= c2_len:
        for i, v1 in enumerate(c1):
            assignment.append([v1, c2[idx_list[i]]])
    else:
        for i, v1 in enumerate(c2):
            assignment.append([v1, c1[idx_list[i]]])

    return assignment

################################################################################

if len(sys.argv) != 4:
    die('Usage')

input, randomizations, output = sys.argv[1:]
main(input, int(randomizations), output)
