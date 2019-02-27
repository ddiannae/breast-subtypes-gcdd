import argparse
import sys
import networkx as nx
from random import sample
import matplotlib.pyplot as plt
import pandas as pd

def getComponentSizes(G) :
    return([len(c) for c in sorted(nx.connected_components(G), key=len, reverse=True)])

def chunk_file(f, chunksize=4096):
    return iter(lambda: f.read(chunksize), b'')


def getPercolationGraph(path, sample_size):

    colnames = ["label", "nodes",  "lc"]
    df = pd.DataFrame(columns=colnames)
    conds = {"basal": "Basal", "luma": "LumA", "lumb": "LumB", "her2": "Her2+", "healthy": "Healthy"}

    for cond in conds.keys():

        print("Working with condition: ", cond)

        graph = nx.Graph()

        for chunk in pd.read_csv(path + cond + '.sif', delimiter="\t",
                                 chunksize=5000, names = ["source", "MI", "target"]):
            Gchunk = nx.from_pandas_edgelist(chunk, "source", "target", "MI")
            print("New graph chunk with ", str(Gchunk.number_of_nodes()), " nodes")
            graph.add_edges_from(Gchunk.edges(data = True))
            print("Entire graph has ", str(graph.number_of_nodes()),
                  " nodes and ",  str(graph.number_of_edges()), " edges")

        sortedEdges = sorted(graph.edges(data=True), key=lambda x: x[2]["MI"])
        LCSizes = []
        nnodes = [graph.number_of_nodes()]
        cSizes = getComponentSizes(graph)
        if cSizes:
            LCSizes.append(cSizes[0])

        while True:
            nToRemove = min(sample_size, graph.number_of_edges())
            edgesToRemove = sortedEdges[0:nToRemove]
            graph.remove_edges_from(edgesToRemove)
            graph.remove_nodes_from(list(nx.isolates(graph)))
            sortedEdges = [e for e in sortedEdges if e not in edgesToRemove]
            cSizes = getComponentSizes(graph)
            if cSizes:
                if cSizes[0] > graph.number_of_nodes()/2:
                    nnodes.append(graph.number_of_nodes())
                    LCSizes.append(cSizes[0])
                else:
                    nnodes.append(graph.number_of_nodes())
                    LCSizes.append(0)
            else:
                break
        df_cond = pd.DataFrame({colnames[0]: [cond] * len(LCSizes),
                                colnames[1]: nnodes,
                                colnames[2]: LCSizes})
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
    plt.ylabel('Nodes in Giant Component')
    plt.xlabel('Nodes in Graph')
    fig.set_size_inches(18.5, 10.5)
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
