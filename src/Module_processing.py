import copy
from utils import *
import scipy.sparse as sp
import numpy as np
import os
import argparse
import pandas as pd
import networkx as nx
from tqdm import tqdm

def computeW():
    gene_to_id = load_gene_to_id()
    weight_out_file = "../out/edge_weights/" + 'weights.txt'
    N = len(gene_to_id)
    w = np.zeros((N, N))

    with open(weight_out_file) as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            w[int(gene_to_id[line[0]])][int(gene_to_id[line[1]])] = line[2]
            w[int(gene_to_id[line[1]])][int(gene_to_id[line[0]])] = line[2]
    return w

# Node coverage
def CO_NODE(v_index):
    # print(node_index)
    # number of patients
    num_samples = len(load_patients_to_indices())
    gene_patients = data[id_to_gene[v_index]]
    co = float(len(gene_patients)) / num_samples
    return co

def load_sample_gene():
    filename = "../data/pan12gene2freq.txt"
    genes = []
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            genes.append(line[0])
    return genes

def cosmic(v):
    fhinput = open('../data/Census_allTue_May_23_12-08-15_2017.tsv')
    cosmic_genes = []
    for line in fhinput:
        cosmic_genes.append(line.split()[0])

    genes = load_gene_list()
    sample_gene = load_sample_gene()
    # 436
    uni = set(genes) & set(sample_gene) & set(cosmic_genes)
    l = len(set(v) & uni)
    return l

def NCG(v):
    NCG_genes = []
    with open('../data/NCG_591.txt','r')as f:
        for line in f.readlines():
            line = line.replace('\n','').split()
            NCG_genes.append(line[0])
    genes = load_gene_list()
    sample_gene = load_sample_gene()
    # 410
    uni = set(genes) & set(sample_gene) & set(NCG_genes)
    l = len(set(v) & uni)
    return l

def Node_Strength(i):
    neighbor = list(G.neighbors(i))
    node_force = 0.0
    for j in range(len(neighbor)):
        node_force += G[i][neighbor[j]]['weight']
    if len(neighbor) == 0:
        return 0
    else:
        return node_force / len(neighbor)
def k_shell(index_m,k_s=1):
    tmp_list11 = []
    subG = nx.Graph()
    for i in index_m:
        G.add_node(i)

    for i in index_m:
        for j in index_m:
            if W[i][j] > 0:
                subG.add_edge(i, j, weight=W[i][j])
    # 遍历每一个连通子图
    for sub_g in nx.connected_components(subG):
        sub_g = subG.subgraph(sub_g)
        sub = nx.k_core(sub_g, k=k_s)
        tmp11 = list(sub.nodes)
        if len(tmp11) >= min_module_size:
            tmp_list11.append(tmp11)
    return tmp_list11[:]


if __name__=='__main__':

    # parse arguments
    description = "Process cluster"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-dn', '--del_num', type=int, required=False, default=1, help="delete numbers")
    parser.add_argument('-k', "--k_threshold", type=int, required=False, default=3,help="min length of module")
    parser.add_argument('-ns', '--numstart', type=int, required=False, default=2500, help="max number of genes")
    parser.add_argument('-ne', '--numend', type=int, required=False, default=100, help="min number of genes")
    parser.add_argument('-step', '--stepsize', type=int, required=False, default=100, help="step size of decrement from num start to num end")
    args = parser.parse_args()


    num_start = args.numstart
    num_end = args.numend
    step =  args.stepsize
    min_module_size = args.k_threshold
    del_num = args.del_num

    # genetic sample data
    data = load_gene_vs_patient_data()
    gene_to_id = load_gene_to_id()
    id_to_gene = load_id_to_gene()

    if not os.path.exists("../out/txt/"):
        os.makedirs("../out/txt/")

    # initialization
    # module all gene index
    modules_gene = []
    # All modules index
    modules_list = []
    modules_end = []

    # build a graph
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
                G.add_edge(i, j, weight=W[i][j])



    cosmic_list = []
    NCG_list = []
    # 初始化
    with open('../out/cu.txt', 'r') as f:
        for line in f.readlines():
            line = line.replace('\n', '')
            tmp_list = line.split('\t')
            del tmp_list[-1]
            if len(tmp_list) >= min_module_size:
                f_tmp = []
                for j in tmp_list:
                    modules_gene.append(gene_to_id[j])
                    f_tmp.append(gene_to_id[j])
                modules_list.append(f_tmp)
        f.close()
    print('初始化:',len(modules_list),modules_list)
    print('初始化:',len(modules_gene))

    for i in nodes:
        if i not in modules_gene:
            G.remove_node(i)


    print('初始化G:', len(G.nodes))
    print("Number of modules:", len(modules_list), "Number of genes:", len(modules_gene), ", G.nodes:",G.number_of_nodes(), " G.edges:", G.number_of_edges())
    while num_start >= num_end:
        # 找到ni 评分最小的点
        ni = []
        for i in modules_gene:
            va = Node_Strength(i) * CO_NODE(i)
            ni.append(va)
        arr = np.array(ni)
        # 找到第 k 小的元素的索引
        partition_index = np.argpartition(arr, del_num)[:del_num]
        # 对前 k 小的元素进行排序，得到它们在原数组中的索引
        indices = partition_index[np.argsort(arr[partition_index])]
        # 输出后 k 小的元素在原数组中的索引
        #print(indices)
        for d in indices:
            for i in modules_list:
                if modules_gene[d] in i:
                    i.remove(modules_gene[d])
        for d in indices:
            # 从 G 中删除该点
            G.remove_node(modules_gene[d])

        # 取子连通图
        tmp_list = []
        for m in modules_list:
            # 遍历每一个连通子图
            subG = G.subgraph(m)
            for sub_g in nx.connected_components(subG):
                tmp11 = list(sub_g)
                if len(tmp11) >= min_module_size:
                    tmp_list.append(tmp11)

        modules_list.clear()
        modules_list = tmp_list[:]

        modules_gene.clear()
        for i in modules_list:
            for j in i:
                modules_gene.append(j)
        nodes.clear()
        nodes = list(G.nodes)
        for i in nodes:
            if i not in modules_gene:
                G.remove_node(i)


        print("Number of modules:",len(modules_list),"Number of genes:",len(modules_gene),", G.nodes:",G.number_of_nodes()," G.edges:",G.number_of_edges())


        if len(modules_gene) <= num_start and len(modules_gene) > num_start-step:
            print('################',num_start,'################')
            cosmic_num = 0
            NCG_num = 0
            for i in modules_list:
                mo = []
                for j in i:
                    mo.append(id_to_gene[j])
                cosmic_num += cosmic(mo)
                NCG_num += NCG(mo)
            cosmic_list.append(cosmic_num)
            NCG_list.append(NCG_num)
            # write to file
            path_txt = "../out/txt/" + str(num_start) + "_" + str(len(modules_list)) + "_" + str(len(modules_gene)) + "_" + str(cosmic_num) +'_'+str(NCG_num) +".txt"
            txtfile = open(path_txt,'w',encoding='utf-8')
            m = 0
            for i in modules_list:
                mm = []
                for j in i:
                    mm.append(id_to_gene[j])
                m1 = "".join(str(i) + ' ' for i in mm)
                txtfile.write(m1)
                m += 1
                txtfile.write('\n')
            txtfile.close()
            num_start -= step

    # Reference Gene save File
    cosmic_note_path = "../out/cosmic.txt"
    cosmic_note = open(cosmic_note_path,'w',encoding='utf-8')
    NCG_note_path = "../out/NCG.txt"
    NCG_note = open(NCG_note_path, 'w', encoding='utf-8')
    for i in reversed(cosmic_list):
        cosmic_note.write(str(i) + '\n')
    for i in reversed(NCG_list):
        NCG_note.write(str(i) + '\n')
    cosmic_note.close()
    NCG_note.close()





