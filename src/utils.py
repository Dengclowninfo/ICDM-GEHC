############################################################
#                                                          #
#                                                          #
####  util.py includes common imports for all scripts   ####
#                                                          #
#                                                          #
############################################################

import numpy as np

GENE_ID_OFFSET = 1  # the genes id starts with 1

network_name = "hint"     # ["hint", "irefindex", "multinet"]
network_beta = 0.4        # [0.4, 0.45, 0.5]

#input data
gene_vs_patiensts_file = 'genes_vs_patients_indices_gen_paran.txt'
gene_score_file = 'pan12gene2freq.txt'
patient_data_file = 'patients_to_indices_gen.txt'
# gene_index_file = 'hint_index_file.txt'
# gene_edge_list_file = 'hint_edge_file.txt'
gene_index_file = 'intersection_hint_index.txt'
gene_edge_list_file = 'intersection_hint_edge.txt'

miRNA_file = 'mirDIP_Results_229135 _gene_6930.txt'

#analysis files
cosmic_file = 'Census_allTue_May_23_12-08-15_2017.tsv'

def load_gene_vs_patient_data():
    data = {}
    with open("../data/" + gene_vs_patiensts_file) as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            data[line[0]] = line[1:]
    return data
    
def load_gene_vs_patient_data_ready():
    data = {}
    with open("../data/"+ gene_score_file) as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            data[line[0]] = float(line[1])
    return data

def load_patients_to_indices():
    data = {}
    with open("../data/" + patient_data_file) as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            data[line[0]] = line[1]
    return data

def load_id_to_gene():
    filename = ""
    if network_name == "hint":
        filename = "../data/intersection_hint_index.txt"
    else:
        filename = "../data/" + gene_index_file

    id_to_gene = {}
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            id_to_gene[int(line[0]) - GENE_ID_OFFSET] = line[1]
    return id_to_gene

def load_gene_list():
    filename = ""
    if network_name == "hint":
        filename = "../data/intersection_hint_index.txt"
    else:
        filename = "../data/" + gene_index_file

    genes = []
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            genes.append(line[1])
    return genes

def load_gene_to_id():
    filename = ""
    if network_name == "hint":
        filename = "../data/intersection_hint_index.txt"
    else:
        filename = "../data/" + gene_index_file

    id_to_gene = {}
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            id_to_gene[line[1]] = int(line[0]) - GENE_ID_OFFSET
    return id_to_gene

def load_unique_genes():
    filename = ""
    if network_name == "hint":
        filename = "../data/intersection_hint_index.txt"
    else:
        filename = "../data/" + gene_index_file

    genes = []
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            genes.append(line[1])
    return genes

def load_edge_list():
    filename = ""
    if network_name == "hint":
        filename = "../data/intersection_hint_edge.txt"
    else:
        filename = "../data/" + gene_edge_list_file

    edge_list = []
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            edge_list.append([int(line[0]) - GENE_ID_OFFSET, int(line[1]) - GENE_ID_OFFSET])
    return edge_list

def write_matrix_file(m, filename):
    with open(filename, "w+") as f:
        for i in range(len(m)):
            for j in range(len(m[i])):
                f.write(str(m[i][j])+" ")
            f.write("\n")
    return

def load_matrix_file(filename):
    n = len(load_gene_list()) # number of genes
    m = np.zeros((n, n))
    with open(filename) as f:
        lines = f.readlines()
        i = 0
        for line in lines:
            line = line.split()
            for j in range(n):
                m[i][j] = line[j]
            i = i + 1
    return m

def read_cosmic_genes():
    fhinput = open("../data/" + cosmic_file)
    cosmic_genes = []
    line = fhinput.readline()
    for line in fhinput:
        cosmic_genes.append(line.split()[0])
    return cosmic_genes

def read_mirDIP():
    mirDIP_list = []
    fhinput = open("../data/" + miRNA_file)
    for line in fhinput.readlines():
        tmp = line.strip().split()
        tmp_re = []
        tmp_re.append(tmp[0])
        tmp_re.append(tmp[2])
        tmp_re.append(float(tmp[3]))
        mirDIP_list.append(tmp_re)
    return mirDIP_list[:]