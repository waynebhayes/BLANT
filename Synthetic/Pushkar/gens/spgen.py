import sys

def main(nodes, edges):
    graph = [set() for i in range(nodes)]
    e = 0
    i = 0
    while e < edges:
        i += 1
        for n in range(0, nodes):
            graph[n].add((n+i)%nodes)
            e += 1
            if e == edges: break
    
    for n in range(0, nodes):
        for x in graph[n]:
            edge = str(n) + " " + str(x)
            print(edge)

if __name__ == "__main__":     
    if len(sys.argv) != 2:
        print("USAGE: python3 spgen.py inputGraph")
        sys.exit()
    
    graphfile = sys.argv[1]
    graphnodes = {}
    graphedges = {}

    for line in open(graphfile, 'r'):
        line = line.strip()
        a,b = tuple(map(int, line.split()))

        graphnodes[a] = True
        graphnodes[b] = True
        edge = str(min(a,b)) + " " + str(max(a,b))
        graphedges[edge] = True 


    nodes = len(graphnodes.keys())
    edges = len(graphedges.keys())
    main(nodes, edges)



