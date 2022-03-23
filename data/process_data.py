import os
import argparse
from platform import node

SPLITOR = '      '

def process_graph(data_graph, output, skip):
    filepath = os.getcwd()

    if not os.path.exists(filepath):
        os.mkdir(filepath)

    filename = os.path.join(filepath, data_graph)
    output = os.path.join(filepath, output)

    edges = []
    with open(filename, 'r') as lines:
        index = skip
        for line in lines:
            if index > 0:
                index -= 1
                continue
            
            line = line.replace('\n', '')
            nodes = line.split(SPLITOR)

            edges.append((nodes[0], nodes[1], nodes[3].replace(' ', '')))

    edges = sorted(edges, key=lambda tup: tup[2])


    u_node = {}
    u_id = 0

    v_node = {}
    v_id = 0

    for (u, v, t) in edges:
        if u not in u_node:
            u_node[u] = u_id
            u_id += 1

        if v not in v_node:
            v_node[v] = v_id
            v_id += 1


    edges_count = 0
    nodes = ""
    for (u, v, t) in edges:
        nodes += "{} {} {}\n".format(u_node[u], v_node[v], t)
        edges_count += 1

    temp = "{} {} {}\n".format(len(u_node), len(v_node), edges_count)
    nodes = temp + nodes

    f = open(output, 'w')
    f.write(nodes)
    f.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Processing graph data 🙂')
    parser.add_argument('--data', '-d', required=True, help='Please add data graph')
    parser.add_argument('--skip', '-s', required=True, help='How many lines to skip')
    parser.add_argument('--output', '-o', required=True, help='Please add the output file')

    args = parser.parse_args()
    process_graph(args.data, args.output, int(args.skip))