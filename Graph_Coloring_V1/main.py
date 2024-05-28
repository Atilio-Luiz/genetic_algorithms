# client program
import networkx as nx
from ga_gc_v1 import ga_graph_coloring

##
# Receives a string with the file name and a boolean value
# that is True iff the graph is 0-indexed
##
def create_graph_from_file(filename, graph_is_0_indexed):
    #G = dict()
    G = nx.Graph()
    # Open the file using the open() function
    with open(filename, 'r') as file:
        # Use a for loop to iterate over each line in the file
        for line in file:
            if line[0] == 'c' or line[0] == 'p':
                continue
            if line[0] == 'e':
                line = line.rstrip()
                e = line.split(' ')
                u = int(e[1])
                v = int(e[2])
                if graph_is_0_indexed == False:
                    u = u - 1
                    v = v - 1
                G.add_edge(u,v)
    return G

##
# Receives a graph and runs the genetic algorithm on it
##
def generate_coloring(graph):
    gagc = ga_graph_coloring(graph, 100, 0.5, 0.1, 0.2, 0.1)
    print(gagc.run(200))

# main function
def main():
    filename = 'dsjc250.5.col'
    G = create_graph_from_file(filename, False)
    print("-----------------")
    print(G)
    max_degree = max(dict(G.degree()).values())
    print("maximum degree of the graph:",max_degree)
    print("algorithm is running .....")
    generate_coloring(G)
    print("-----------------")


if __name__ == "__main__":
    main()