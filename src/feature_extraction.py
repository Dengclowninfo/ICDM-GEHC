#!/usr/bin/env python
from utils import *
from DIANA import *
import numpy as np
import argparse
import pandas as pd
import networkx as nx
from node2vec import Node2Vec


def computeW():
    gene_to_id = load_gene_to_id()
    weight_out_file = "../out/edge_weights/" + 'weights.txt'
    N = len(gene_to_id)
    w = np.zeros((N, N))
    with open(weight_out_file) as f:
        header = f.readline()
        lines = f.readlines()
        for line in lines:
            line = line.split()
            w[int(gene_to_id[line[0]])][int(gene_to_id[line[1]])] = line[2]
            w[int(gene_to_id[line[1]])][int(gene_to_id[line[0]])] = line[2]
    return w

def load_sample_gene():
    filename = "../data/pan12gene2freq.txt"
    genes = []
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            genes.append(line[0])
    return genes

def load_gene_list():
    filename = "../data/hint_index_file.txt"
    genes = []
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            genes.append(line[1])
    return genes

if __name__ == "__main__":
    network_file = 'hint'

    # parse arguments
    description = "Feature representation learning and clustering on edge-weighted graphs"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-d', '--dimensions', type=int, required=False, default=32, help="Node2vec_dimensions")
    parser.add_argument('-p', '--return_p', type=float, required=False, default=4, help="Node2vec_return_parameter")
    parser.add_argument('-q', '--leave_q', type=float, required=False, default=1, help="Node2vec_leave_parameter")
    parser.add_argument('-wl', '--walk_length', type=int, required=False, default=60, help="Node2vec_walk_length_parameter")
    parser.add_argument('-nw', '--num_walks', type=int, required=False, default=300, help="Node2vec_num_walks_parameter")
    parser.add_argument('-k', '--clusters_number', type=int, required=False, default=127,
                        help="DIANA_Number_of_clusters")




    args = parser.parse_args()

    d = args.dimensions
    return_p = args.return_p
    leave_q = args.leave_q
    walk_length = args.walk_length
    num_walks = args.num_walks
    k = args.clusters_number

    weight_out_file = "../out/edge_weights/" + 'weights.txt'

    W = computeW()
    id_to_gene = load_id_to_gene()
    N = len(id_to_gene)
    G = nx.Graph()

    for i in range(N):
        G.add_node(i)
    nodes = list(G.nodes)
    for i in nodes:
        for j in nodes:
            if W[i][j] > 0:
                G.add_edge(i, j,weight=W[i][j])

    node2vec = Node2Vec(
        G,
        dimensions = d,
        p = return_p,
        q = leave_q,
        walk_length = walk_length,
        num_walks = num_walks,
        workers = 1,
        weight_key = 'weight',
        seed = 8
    )
    #
    model = node2vec.fit(window=3,min_count=1,batch_words=4)
    X = model.wv.vectors
    # print(X.shape)

    pd.DataFrame(X).to_csv('../out/X.csv')

    # sp_e = sp.csc_matrix(X)
    # e_path = "../out/sparse_matrix_X.npz"
    # sp.save_npz(e_path, sp_e)

    model.wv.save_word2vec_format('../out/EMBEDDING_FILENAME.vector')

    # k
    # DIANA
    print('DIANA...')
    clustAssing = DIANA(X, k)
    cluster_labels_2 = [n for n in range(0, X.shape[0])]
    for i in range(0, k):
        for j in clustAssing[i]:
            cluster_labels_2[j] = i
    colors = []
    # Traverse each node in networkx order
    for node in nodes:
        # Get the index number of this node in embedding
        idx = model.wv.key_to_index[str(node)]
        # Get the clustering result of this node
        colors.append(cluster_labels_2[idx])

    n_clu = len(set(colors))
    print(n_clu)
    id_to_gene = load_id_to_gene()
    fp2 = open('../out/cu.txt', 'w', encoding='utf-8')
    model_genes2 = set()
    for flag in range(n_clu):
        for i in range(len(nodes)):
            if colors[i] == flag:
                fp2.write(id_to_gene[nodes[i]] + '\t')
                model_genes2.add(id_to_gene[nodes[i]])
        fp2.write('\n')