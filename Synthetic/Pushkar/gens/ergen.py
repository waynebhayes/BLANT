import clean
import os
import sys

#script_path = "/home/pushkarg/networkprops/codebase/"
script_path = os.path.dirname(os.path.realpath(__file__))

sys.path.insert(0, script_path + "/networkx/networkx-2.2/networkx/")
import networkx

def create(input_file_name):
    output_file_name = "/tmp/" + os.path.basename(input_file_name) + str(os.getpid())

    # node-map
    nodemap = []  # list of lists
    clean.clean(input_file_name, output_file_name, nodemap)  # names inserted into nodemap
    networkx_graph = networkx.read_edgelist(output_file_name)
    
    # delete the temporary file
    os.remove(output_file_name)

    nodes = len(nodemap)
    edges = len(networkx_graph.edges)
    
    ergraph = networkx.generators.random_graphs.erdos_renyi_graph(nodes, edges/(nodes * (nodes-1)/2 ))
    for line in networkx.generate_edgelist(ergraph, data=False):
        print(line)
    
if __name__ == "__main__":    
    error = False
    if len(sys.argv) == 2:
        filename = sys.argv[1].strip()
        create(filename)
    else:
        error = True

    if error: 
        print("Input error! Usage: python3 ergen.py inputFile")












