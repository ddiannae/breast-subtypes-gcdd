import argparse
import sys
import networkx as nx
from random import sample
import matplotlib.pyplot as plt
import pandas as pd

def getLCsizes(G) :
    lcsizes = [len(c) for c in sorted(nx.connected_components(G), key=len, reverse=True)]
    return(lcsizes)

def getPercolationGraph(path, sample_size):

    colnames = ["label", "nodes",  "lc"]
    df = pd.DataFrame(columns=colnames)
    conds = {"basal": "Basal", "luma": "LumA", "lumb": "LumB", "her2" : "Her2+", "healthy": "Healthy"}

    for cond in conds.keys():

        print("Working with condition: ", cond)

        graph = nx.Graph()
        edges = nx.read_edgelist(path + cond + ".sif", delimiter = "\t",  data=(('MI',float), ))
        graph.add_edges_from(edges.edges())

        maxlcsizes = []
        nnodes = [graph.number_of_nodes()]
        lcsizes = getLCsizes(graph)
        if lcsizes:
            maxlcsizes.append(lcsizes[0])

        while True:
            out_nodes = sample(graph.nodes(), min(sample_size, graph.number_of_nodes()))
            graph.remove_nodes_from(out_nodes)
            lcsizes = getLCsizes(graph)
            if lcsizes:
                nnodes.append(graph.number_of_nodes())
                maxlcsizes.append(lcsizes[0])
            else:
                break
        df_cond = pd.DataFrame({colnames[0]: [cond] * len(maxlcsizes),
                                colnames[1]: nnodes,
                                colnames[2]: maxlcsizes})
        df = df.append(df_cond, sort=False)

    print("Saving plot")

    by_label = df.groupby(colnames[0])
    fig = plt.figure()
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))

    for name, group in by_label:
        plt.plot(group[colnames[1]], group[colnames[2]], label=conds[name])

    plt.legend()
    plt.gca().invert_xaxis()
    plt.grid(True)
    plt.ylabel('Nodes in Largest Component')
    plt.xlabel('Nodes in Graph')
    fig.savefig('Percolation.png', dpi=fig.dpi)

def main(args):
    parser = argparse.ArgumentParser(description="Do something.")
    parser.add_argument("-p", "--path", type=str, default= "./", required=False)
    parser.add_argument("-s", "--size", type=int, default= 50, required=False)
    args = parser.parse_args(args)
    print("Searching for files in ", args.path)
    print("Removing ", args.size, " nodes per iteration")
    getPercolationGraph(args.path, args.size)

if __name__ == "__main__":
    main(sys.argv[1:])
