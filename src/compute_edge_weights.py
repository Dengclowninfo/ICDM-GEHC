#!/usr/bin/env python

from utils import *
import networkx as nx
import os
import argparse

#compute mutex values
def compute_mutex_nsep():
    G = nx.Graph()
    for e in edge_list:
         G.add_edge(e[0], e[1])
    G.to_undirected()
    # G.remove_edges_from(G.selfloop_edges())
    G.remove_edges_from(nx.selfloop_edges(G))
    N = len(genes)
    num_samples = len(load_patients_to_indices())  # number of patients
    weight_out_file = mex_nsep_out_file
    fhout = open(weight_out_file, 'w+')
    count_1s = 0
    count_0s = 0
    for i in range(len(edge_list)):
        e = edge_list[i]
        gene1_index = e[0]
        gene2_index = e[1]
        gene1 = id_to_gene[gene1_index]
        gene2 = id_to_gene[gene2_index]

        gene1_neighbours = G[gene1_index]
        gene2_neighbours = G[gene2_index]

        union_set1 = set()
        union_set2 = set()
        count1 = 0
        count2 = 0

        if gene1 in data:
            union_set1 = set(data[gene1])
            count1 = len(data[gene1])
        if gene2 in data:
            union_set2 = set(data[gene2])
            count2 = len(data[gene2])

        for gene in gene1_neighbours:
            if id_to_gene[gene] in data and gene != gene2_index:  # exclude gene2 ?
                union_set1 = union_set1.union(data[id_to_gene[gene]])
                count1 += len(data[id_to_gene[gene]])
        for gene in gene2_neighbours:
            if id_to_gene[gene] in data and gene != gene1_index:
                union_set2 = union_set2.union(data[id_to_gene[gene]])
                count2 += len(data[id_to_gene[gene]])



        union_set1 = len(union_set1)
        union_set2 = len(union_set2)

        m1 = 0
        m2 = 0
        if count1 != 0:
            m1 = float(union_set1) / count1
        if count2 != 0:
            m2 = float(union_set2) / count2

        mutex = (m1+m2)/2

        if mutex == 1:
           count_1s += 1
        print(id_to_gene[e[0]] + "\t" + id_to_gene[e[1]] + "\t" + str(mutex),file=fhout)


    fhout.close()


#compute coverage values
def compute_cov():
    G = nx.Graph()
    for e in edge_list:
         G.add_edge(e[0], e[1])

    G.to_undirected()
    # G.remove_edges_from(G.selfloop_edges())
    G.remove_edges_from(nx.selfloop_edges(G))
    N = len(genes)

    weight_out_file = cov_out_file
    fhout = open(weight_out_file, 'w+')

    for i in range(len(edge_list)):
        e = edge_list[i]
        gene1_index = e[0]
        gene2_index = e[1]
        gene1 = id_to_gene[gene1_index]
        gene2 = id_to_gene[gene2_index]

        gene1_patients = []
        gene2_patients = []

        if gene1 in data:
            gene1_patients = data[gene1]
        if gene2 in data:
            gene2_patients = data[gene2]

        gene1_count = len(gene1_patients)
        gene2_count = len(gene2_patients)

        if gene1_count + gene2_count != 0:
            gene1_cover = float(gene1_count) / num_samples
            gene2_cover = float(gene2_count) / num_samples

            cov = gene1_cover * gene2_cover

            print(id_to_gene[e[0]]+ "\t" + id_to_gene[e[1]] + "\t" + str(cov),file = fhout)
    fhout.close()

# Calculate the confidence score with set threshold
def miRNA_SSC(g1,g2):
    g1_mir_list = []
    g2_mir_list = []
    g1_mir_set = set()
    g2_mir_set = set()
    g1_g2_mir_set = set()
    for mir in mirDIP:
        if g1 == mir[0]:
            g1_mir_list.append(mir)
            g1_mir_set.add(mir[1])
        if g2 == mir[0]:
            g2_mir_list.append(mir)
            g2_mir_set.add(mir[1])

    g1_g2_mir_set = g1_mir_set & g2_mir_set
    sum_sc = 0.0
    for g1_mir in g1_mir_list:
        if g1_mir[1] in g1_g2_mir_set:
            sum_sc += g1_mir[2]
    for g2_mir in g2_mir_list:
        if g2_mir[1] in g1_g2_mir_set:
            sum_sc += g2_mir[2]
    g1_g2_edge = 2 * len(g1_g2_mir_set)
    if len(g1_g2_mir_set) != 0:
        ans = float(sum_sc / g1_g2_edge)
    else:
        ans = ts_threshold
    return ans


#compute edge weights as product of mutex ,confidence and coverage values with set threshold
def compute_edge_weights():
    N = len(genes)
    mutex_nsep_scores = {}
    with open(path_pre + 'mutex_nsep.txt', "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip().split()
            mutex_nsep_scores[line[0]+" "+line[1]] = float(line[2])

    cov_scores = {}
    with open(path_pre + 'cov.txt', "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip().split()
            cov_scores[line[0]+" "+line[1]] = float(line[2])

    weight_out_file = edge_out_file
    fhout = open(weight_out_file, 'w+')

    for i in range(len(edge_list)):
        e = edge_list[i]
        gene1_index = e[0]
        gene2_index = e[1]
        gene1 = id_to_gene[gene1_index]
        gene2 = id_to_gene[gene2_index]
        gene1_patients = []
        gene2_patients = []

        if gene1 in data:
            gene1_patients = data[gene1]
        if gene2 in data:
            gene2_patients = data[gene2]

        gene1_count = len(gene1_patients)
        gene2_count = len(gene2_patients)

        if gene1_count + gene2_count != 0:


            mutex = mutex_nsep_scores[gene1 + ' ' + gene2]
            cov = cov_scores[gene1 + ' ' + gene2]
            if mutex < mex_threshold:
                mutex = 0
            res = mutex * cov
            isc = miRNA_SSC(id_to_gene[e[0]],id_to_gene[e[1]])
            res *= isc
            print(id_to_gene[e[0]] + "\t" + id_to_gene[e[1]] + "\t" + str(res),file=fhout)

    fhout.close()

def load_sample_gene():
    filename = "../data/pan12gene2freq.txt"
    genes = []
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            genes.append(line[0])
    return genes


if __name__ == '__main__':

    # parse arguments
    description = "Calculate the edge weight of two genes based on coverage, mutual exclusion, and confidence"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-t', '--mex_threshold', type=float, required=False, default=0.7, help="mutex_nsep_threshold")
    parser.add_argument('-ts', '--ts_threshold', type=float, required=False, default=0.2, help="ts_threshold")
    args = parser.parse_args()

    mex_threshold = args.mex_threshold

    res = 0
    n = 0
    with open('../data/mirDIP_Results_229135 _gene_6930.txt', 'r') as f:
        line = f.readline()
        while line:
            line = line.replace('\n', '').split('\t')
            res += float(line[4])
            n += 1
            line = f.readline()
    ts_threshold = round(args.ts_threshold * res / n,1)
    print("ts_threshold:",ts_threshold)
    # creating path
    path_pre = '../out/edge_weights/pre/'
    main_path = '../out/edge_weights/'

    if not os.path.exists(path_pre):
        os.makedirs(path_pre)

    # input files
    mex_nsep_out_file = path_pre + 'mutex_nsep.txt'
    cov_out_file = path_pre + 'cov.txt'
    edge_out_file = main_path + 'weights.txt'

    # load data
    num_samples = len(load_patients_to_indices())  # number of patients
    data = load_gene_vs_patient_data()
    genes = load_unique_genes()
    id_to_gene = load_id_to_gene()
    gene_to_id = load_gene_to_id()  # gene string to indices
    edge_list = load_edge_list()
    mirDIP = read_mirDIP()

    geness = load_gene_list()
    sample_gene = load_sample_gene()
    union = set(geness) & set(sample_gene)

    # print(id_to_gene)

    #computing edge weights
    compute_mutex_nsep()
    compute_cov()
    compute_edge_weights()
