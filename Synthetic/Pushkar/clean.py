def clean(input_file_name, output_file_name, id_name_store):
    file_in = open(input_file_name, 'r')
    file_ot = open(output_file_name, 'w')

    node_id = 0
    node_map = {}  # name:id

    graph = []  # a list of sets (this set contains edges)
    #print("cleaning!")
    for line in file_in.readlines():
        line = line.strip()
        if (len(line) == 0) or (line[0] == '#') or ((line.count(' ') + line.count('\t')) != 1):
            print("line skipped")
            continue

        node1, node2 = tuple(map(str, line.split()))

        if node1 not in node_map:
            node_map[node1] = node_id  # saving name-id, inside local ds
            id_name_store.append([node1])  # saving id-name, inside the passed argument
            node1 = node_id
            graph.append(set([]))
            node_id += 1
        else:
            node1 = node_map[node1]

        if node2 not in node_map:
            node_map[node2] = node_id  # saving name-id, inside local ds
            id_name_store.append([node2])  # saving id-name, inside the passed argument
            node2 = node_id
            graph.append(set([]))
            node_id += 1
        else:
            node2 = node_map[node2]

        if node1 < node2:
            graph[node1].add(node2)
        elif node2 <= node1:     # self-loops are ok (node1 == node2)
            graph[node2].add(node1)

    # now write to the output file
    for i, edgeset in enumerate(graph):
        for neighbor in edgeset:
            line = str(i) + " " + str(neighbor) + "\n"
            file_ot.write(line)

    #print("cleaned..", "elts in id_name_store:", len(id_name_store), ". elts in internal structs: ",
          #len(node_map), len(graph))


