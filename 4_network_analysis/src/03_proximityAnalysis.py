#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''

@author: adapted by Gehad Youssef  based on Han et al., 2021, original code created by woochanghwang 

Running time: depends network size an number of permutations

# main
    # sub: create_network_from_sif_file - network from entrez gene interactome database
    # sub: get_drug_target_drugbank - drug_to_genes_dict[drug] = {genes}
    # sub: get_diseasome_genes_from_selectedGenes - disease_to_gene[disease]={genes}, JUST ONE DISEASE: HNSCC
    # sub: calculate_proximity_multiple_multiprocessing_return
        # sub: get_degree_bimning - bims = [(1, 1, ['154754',....]), (2, 2, ['geneID', 'geneID, ...]), ... ]
        # sub: calculate_proximity
            # sub: calculate_closest_distance - shortest path length between nodes_from to nodes_to
                # sub: check_has_path
                # sub: get_shortest_path_length_between
            # sub: get_random_nodes
                # sub: pick_random_nodes_matching_selected
                    # sub: get_degree_equivalents
            # sub: get_degree_bimning
            

'''
#%%
import time
import pandas as pd
import numpy as np
import networkx, random
import mygene
import multiprocessing as mp
#%%
""" CHECK LINE 413
    if you find that mygene doesnt find all genes: 
    abort run, and enter entrez gene ID (human) manually in line 413
    run again
"""

""" variables """
### ENTER NUMBER OF PROCESSES FOR POOL PROCESSING (KEEP BELOW 3 ON A NORMAL LAPTOP) ###
num_process = 8
### ENTER NUMBER OF PERMUTATIONS (SHOULD BE 1,000 BUT CAN DO 100 FOR QUICKER RUNS) ###
perm_num = 1000
### ENTER IDENTIFIER FOR THE DATASET, IT DOES NOT INFLUENCE THE ANALYSIS IN ANY WAY ###
disease_name = "ccRCC"

""" input addresses """
### ENTER PATH/NAME INPUT: HUMAN PROTEIN INTERACTOME DATABASE ###
address_network_file = "../databases/human_protein_interactome.txt"
### ENTER PATH/NAME INPUT: LIST KEY PROTEINS ###
address_gene_list_file = "../data/key_proteins_muanually_filtered.txt"
### ENTER PATH/NAME INPUT: DRUG BANK ###
address_drug_target_file = "/home/sr933/rcc/4_network_analysis/databases/9606.stitch_full_drug2target.tsv"
#address_drug_target_file = "../databases/.v5.1.5_interactions_900_th_400.onlyTarget.Approved.tsv"
""" output addresses """
### ENTER PATH/NAME OUTPUT: PROXIMITY ANALYSIS RESULTS AS .csv ###
address_output_file = "../data/drug_network_proximity_results_manually_filtered_all.csv" 
### ENTER PATH/NAME OUTPUT: LIST OF SIGNIFICANTLY PROXIMAL DRUGS .txt ###
address_top_drugs = '../data/significantly_proximal_drugs.txt'
### ENTER PATH/NAME OUTPUT: SIGNIFICANTLY PROXIMAL DRUGS AND THEIR TARGETS AS GENE SYMBOL .csv ###
address_top_drug_to_target = '../data/top_drugs_to_target.csv'
#%%

def check_has_path(G,node_from,node_to):
    return(networkx.has_path(G,node_from,node_to))

def get_shortest_path_length_between(G, source_id, target_id):
    return networkx.shortest_path_length(G, source_id, target_id)

def calculate_closest_distance(network, nodes_from, nodes_to):
    values_outer = []
    d_from_to_list = []
    for node_from in nodes_from:
        values = []
        for node_to in nodes_to:
            # print("from - to", node_from, node_to)
            if not check_has_path(network,node_from,node_to): continue
            val = get_shortest_path_length_between(network, node_from, node_to)
            d_from_to_list.append([val, node_from, node_to])
            values.append(val)
        if len(values) == 0:    
            continue
        d = min(values)
        values_outer.append(d)
    d_from_to_df = pd.DataFrame(d_from_to_list)
    d_from_to_df_shortest_0 = d_from_to_df.loc[d_from_to_df[0] == 0]
    d_from_to_df_shortest_1 = d_from_to_df.loc[d_from_to_df[0] == 1]
    d_from_to_df_shortest_2 = d_from_to_df.loc[d_from_to_df[0] == 2]
    # make list of genes where path length is certain length, removing duplicates from list
    targeted_key_proteins_list_0 = list(dict.fromkeys(list(d_from_to_df_shortest_0[2])))
    targeted_key_proteins_list_1 = list(dict.fromkeys(list(d_from_to_df_shortest_1[2])))
    targeted_key_proteins_list_2 = list(dict.fromkeys(list(d_from_to_df_shortest_2[2])))
    d = np.mean(values_outer)
    return d, targeted_key_proteins_list_0, targeted_key_proteins_list_1, targeted_key_proteins_list_2

def calculate_closest_distance_random(network, nodes_from, nodes_to):
    values_outer = []
    for node_from in nodes_from:
        values = []
        for node_to in nodes_to:
            # print("from - to", node_from, node_to)
            if not check_has_path(network,node_from,node_to): continue
            val = get_shortest_path_length_between(network, node_from, node_to)
            values.append(val)
        if len(values) == 0:    
            continue
        d = min(values)
        values_outer.append(d)
    d = np.mean(values_outer)
    return d

def get_degree_bimning(g, bim_size, lengths=None):
    '''

    It tried to make number of gene list( with same degree) to bim size.
    If number of gene list with some degree , it combime with ohter genes with another degree to meet bim size.
    '''
    degree_to_nodes = {}
    for node, degree in g.degree():  # .items(): # iterator in networkx 2.0
        if lengths is not None and node not in lengths:
            continue
        degree_to_nodes.setdefault(degree, []).append(node)
    values = degree_to_nodes.keys()
    # values.sort() #pytyon 2.x
    values = sorted(values)
    bims = []
    i = 0
    while i < len(values):
        low = values[i]
        val = degree_to_nodes[values[i]]
        while len(val) < bim_size:
            i += 1
            if i == len(values):
                break
            val.extend(degree_to_nodes[values[i]])
        if i == len(values):
            i -= 1
        high = values[i]
        i += 1
        # print i, low, high, len(val)
        if len(val) < bim_size:
            low_, high_, val_ = bims[-1]
            bims[-1] = (low_, high, val_ + val)
        else:
            bims.append((low, high, val))
    # bims = [(1, 1, ['154754',....]), (2, 2, ['geneID', 'geneID, ...]), ... ]
    return bims

def get_degree_equivalents(seeds, bims, g):
    seed_to_nodes = {}
    for seed in seeds:
        d = g.degree(seed)
        for l, h, nodes in bims:
            if l <= d and h >= d:
                mod_nodes = list(nodes)
                mod_nodes.remove(seed)
                seed_to_nodes[seed] = mod_nodes
                break
    return seed_to_nodes

def pick_random_nodes_matching_selected(network, bims, nodes_selected, n_random, degree_aware=True, connected=False,
                                        seed=None):
    """
    Use get_degree_bimning to get bims
    """
    if seed is not None:
        random.seed(seed)
    values = []
    nodes = network.nodes()
    for i in range(n_random):
        if degree_aware:
            if connected:
                raise ValueError("Not implemented!")
            nodes_random = set()
            node_to_equivalent_nodes = get_degree_equivalents(nodes_selected, bims, network)
            for node, equivalent_nodes in node_to_equivalent_nodes.items():
                # nodes_random.append(random.choice(equivalent_nodes))
                chosen = random.choice(equivalent_nodes)
                for k in range(20):  # Try to find a distinct node (at most 20 times)
                    if chosen in nodes_random:
                        chosen = random.choice(equivalent_nodes)
                nodes_random.add(chosen)
            nodes_random = list(nodes_random)
        else:
            if connected:
                nodes_random = [random.choice(nodes)]
                k = 1
                while True:
                    if k == len(nodes_selected):
                        break
                    node_random = random.choice(nodes_random)
                    node_selected = random.choice(network.neighbors(node_random))
                    if node_selected in nodes_random:
                        continue
                    nodes_random.append(node_selected)
                    k += 1
            else:
                nodes_random = random.sample(nodes, len(nodes_selected))
        values.append(nodes_random)
    print('..picked random nodes')
    return values

def get_random_nodes(nodes, network, bims=None, n_random=perm_num, min_bim_size=100, degree_aware=True, seed=None):
    if bims is None:
        # Get degree bims of the network
        bims = get_degree_bimning(network, min_bim_size)
    nodes_random = pick_random_nodes_matching_selected(network, bims, nodes, n_random, degree_aware,
                                                                         seed=seed)
    return nodes_random

def calculate_proximity(network, nodes_from, nodes_to, nodes_from_random=None, nodes_to_random=None, bims=None,
                        n_random=perm_num, min_bim_size=100, seed=452456, lengths=None):
    """
    Calculate proximity from nodes_from to nodes_to
    If degree bimning or random nodes are not given, they are generated
    lengths: precalculated shortest path length dictionary
    """

    nodes_network = set(network.nodes())
    nodes_from = set(nodes_from) & nodes_network
    nodes_to = set(nodes_to) & nodes_network
    
    print('..nodes_from: ',  nodes_from)
    print('..nodes_to: ',  nodes_to)
    
    if len(nodes_from) == 0 or len(nodes_to) == 0:
        return None  # At least one of the node group not in network
    d, targeted_key_proteins_list_0, targeted_key_proteins_list_1, targeted_key_proteins_list_2 = calculate_closest_distance(network, nodes_from, nodes_to)
    if bims is None and (nodes_from_random is None or nodes_to_random is None):
        bims = get_degree_bimning(network, min_bim_size,
                                                    lengths)  # if lengths is given, it will only use those nodes
    
    if nodes_from_random is None:
        nodes_from_random = get_random_nodes(nodes_from, network, bims=bims, n_random=n_random,
                                             min_bim_size=min_bim_size, seed=seed)
    if nodes_to_random is None:
        nodes_to_random = get_random_nodes(nodes_to, network, bims=bims, n_random=n_random, min_bim_size=min_bim_size,
                                           seed=seed)
        
    random_values_list = zip(nodes_from_random, nodes_to_random)
    values = np.empty(len(nodes_from_random))  # n_random
    for i, values_random in enumerate(random_values_list):
        nodes_from, nodes_to = values_random
        # values[i] = network_utilities.get_separation(network, lengths, nodes_from, nodes_to, distance, parameters = {})
        values[i] = calculate_closest_distance_random(network, nodes_from, nodes_to)
    # pval = float(sum(values <= d)) / len(values) # needs high number of n_random
    m, s = np.mean(values), np.std(values)
    if s == 0:
        z = 0.0
    else:
        z = (d - m) / s
    return d, z, (m, s), targeted_key_proteins_list_0, targeted_key_proteins_list_1, targeted_key_proteins_list_2 # (z, pval)

def calculate_proximity_multiple_multiprocessing_return(input_tuple):
    """
    Run proximity on each entries of from and to files in a pairwise manner
    output is saved in out_file (e.g., output.txt)
    disease_mode : whole - from disease_gene.tsv
    disease_mode : disease_name - from our analysis
    """
    print('..')
    
    import os
    print('process id = {}'.format(os.getpid()))
    
    network, disease_to_genes, drug_to_targets = input_tuple
    drug = drug_to_targets[0]
    nodes_from = drug_to_targets[1]

    start_time = time.time()

    n_random = perm_num
    min_bim_size = 100
    seed = 452456
    lengths = None

    # Get degree bimning
    print('..get degree bimning.')
    bims = get_degree_bimning(network, min_bim_size)

    print('..start calculate_proximity.')
    proximity_result = []
    try:
        for disease, nodes_to in disease_to_genes.items():
            # print(drug, disease)
            d, z, (m, s), targeted_key_proteins_list_0, targeted_key_proteins_list_1, targeted_key_proteins_list_2 = calculate_proximity(network, nodes_from, nodes_to, nodes_from_random=None,
                                               nodes_to_random=None, bims=bims, n_random=n_random,
                                               min_bim_size=min_bim_size, seed=seed, lengths=lengths)
            proximity_result.append((drug, disease, len(nodes_from), len(nodes_to), d, z, targeted_key_proteins_list_0, targeted_key_proteins_list_1, targeted_key_proteins_list_2))
            #print(proximity_result)
            print("-> {} times training Runtime: {:2f} Minutes".format(drug, ((time.time() - start_time) / 60)))
            return proximity_result[0]
    except:
        print('Oops! Error getting proximity results!')

def create_graph(directed=False):
    """
        Creates & returns a graph
    """
    if directed:
        g = networkx.DiGraph()
    else:
        g = networkx.Graph()
    return g

def create_network_from_sif_file(network_file_in_sif, delim=','):
    
    setEdge = set()
    f = open(network_file_in_sif)
    for line in f:
        if delim is None:
            words = line.rstrip("\n").split()
        else:
            words = line.rstrip("\n").split(delim)
        id1 = words[0]
        id2 = words[1]
        setEdge.add((id1, id2))
    f.close()
    
    if len(setEdge) == 0:
        setEdge = None
    
    g = create_graph()
    g.add_edges_from(setEdge)
    return g

def get_connected_components(G, return_as_graph_list=True):
    """
        Finds (strongly in the case of directed network) connected components of graph
        returnAsGraphList: returns list of graph objects corresponding to connected components (from larger to smaller)
        otherwise returns list of node list corresponding nodes in connected components
    """
    result_list = []

    if return_as_graph_list:
        result_list = networkx.connected_component_subgraphs(G)
    else:
        result_list = [c for c in sorted(networkx.connected_components(G), key=len, reverse=True)]

    return result_list

def get_drug_target_drugbank(address_drug_target_file):
    """
    drugbank file = [[drugid,entrez],[]...]
    group by drugid
    return = drug = gene_set
    """

    drug_target_df = pd.read_table(address_drug_target_file,sep='\t')
    drug_target_df = drug_target_df.applymap(str)
    drug_targets_df = drug_target_df.groupby('DrugID').agg({'Drug_Target':list})
    #print('\n \nDrug target df:\n ', drug_targets_df.head())
    drug_targets_dict = drug_targets_df.to_dict()['Drug_Target']

    drug_to_genes_dict = {}
    ## Check whether targets are in network

    for drug, genes in drug_targets_dict.items():
        genes = set(genes)
        drug_to_genes_dict[drug] = genes
        
    return drug_to_genes_dict

def mapping_genes_id_to(gene_list,id_from, id_to, species = 9606, as_dataframe = False):
    mg = mygene.MyGeneInfo()
    return mg.querymany(gene_list, scopes=id_from ,fields= id_to ,species=species,as_dataframe=as_dataframe)

def get_diseasome_genes_from_selectedGenes(address_gene_list_file, disease_name):
    """
    key gene names are converted to entrez gene ID using mygene
    disease dictionary is returned (just one disease)
    :param address_gene_list_file: result from my anlysis (permutation tests)
    :param disease_name: pre defined
    :return: dict(disease_to_gene[disease]={genes})
    """
    # make list of key genes found in permutation test
    gene_df = pd.read_csv(address_gene_list_file, sep='\t')
    gene_list = list(gene_df['Gene'])
    #print("gene list: ", gene_list)
    
    # translate to entrez gene
    id_from='symbol'
    id_to = 'entrezgene'
    geneID_df = mapping_genes_id_to(gene_list,id_from=id_from, id_to=id_to, species = 9606, as_dataframe = True)    
    geneID_list = list(geneID_df[id_to])
    #print("gene ID list (translated from gene name): \n", geneID_list)
    
    # dictionary: dict(disease_to_gene[disease]={genes})
    disease_to_genes = dict()
    genes = set(geneID_list)
    disease_to_genes[disease_name] = genes

    return disease_to_genes

#%%
def main():
        
    print("Start network-based proximity analysis multiprocessing. Program will tell you when finished.")
    
    #make network from whole human interactome
    #print('create_network_from_sif_file')
    network = create_network_from_sif_file(address_network_file, delim=',')
    
    # fetch drug_to_genes_dict dictionary using address_drug_target_file
    drug_to_targets = get_drug_target_drugbank(address_drug_target_file)
    
    # fetch dict(disease_to_gene[disease]={genes}) with entrez gene ID
    disease_to_genes = get_diseasome_genes_from_selectedGenes(address_gene_list_file,disease_name)
    
    ### CAN ADD ENTREZ IDs MANUALLY HERE ###
    #disease_to_genes[disease_name].add('51373')
    
    # convert dictionary to list
    drug_to_targets_list = []
    for drug, targets in drug_to_targets.items():
        drug_to_targets_list.append([drug,targets])
    print('Length list -drug to targets-:', len(drug_to_targets_list))
    print('Length dict -drug to targets-:', len(drug_to_targets))
    
    # duplicate networks according to number of drugs#
    network_list = [network] * len(drug_to_targets_list)
    
    # duplicate disease to gene dict according to number of drugs#
    disease_to_genes_list = [disease_to_genes] * len(drug_to_targets_list)
    
    # start proximity analysis - pool processing
    print('Start proximity analysis..')
    pool = mp.Pool(processes=num_process)
    start_time = time.time()
    results = pool.imap(calculate_proximity_multiple_multiprocessing_return,
                        zip(network_list, disease_to_genes_list, drug_to_targets_list))
    pool.close()
    
    results = list(results)
    results = filter(lambda x: x is not None, results)
    results_df = pd.DataFrame(columns=["Drug", "Disease", "n.source", "n.target", "d", "z","0degree.target.keyprotein", "1degree.target.keyprotein", "2degree.target.keyprotein"], data=results)
    
    # calculate coverage and specificity scores
    try:
        results_df["n.0degree"] = results_df['0degree.target.keyprotein'].map(len)
        results_df["n.1degree"] = results_df['1degree.target.keyprotein'].map(len)
        results_df["n.2degree"] = results_df['2degree.target.keyprotein'].map(len)
        results_df["coverage.score"] = ( (results_df["n.0degree"]*3) + (results_df["n.0degree"]*2) + (results_df["n.0degree"]) ) / (results_df["n.target"]*6)
        results_df["specificity.score"] = results_df["coverage.score"] / results_df["n.source"]
    except TypeError:
        print("Oops! Error while scoring")
    
    
              
    # make dictionary for target proteins
    id_from='entrezgene'
    id_to = 'symbol'
    genesymbol_df = mapping_genes_id_to(list(disease_to_genes[disease_name]),id_from=id_from, id_to=id_to, species = 9606, as_dataframe = True)
    genesymbol_dict = genesymbol_df.loc[:,'symbol'].to_dict()

    # translate target from entrez to gene name using dictionary
    for count, row in enumerate(results_df["0degree.target.keyprotein"]):
        if results_df["0degree.target.keyprotein"].iloc[count] != []:
            genesymbols = []
            try:
                if list(row) != None:
                    for element in list(row):
                        genesymbols += [genesymbol_dict[element]]
                    results_df["0degree.target.keyprotein"].iloc[count] = [genesymbols]
            except TypeError:
                print("Oops! Error while translating")
    for count, row in enumerate(results_df["1degree.target.keyprotein"]):
        if results_df["1degree.target.keyprotein"].iloc[count] != []:
            genesymbols = []
            try:
                if list(row) != None:
                    for element in list(row):
                        genesymbols += [genesymbol_dict[element]]
                    results_df["1degree.target.keyprotein"].iloc[count] = [genesymbols]
            except TypeError:
                print("Oops! Error while translating")
    for count, row in enumerate(results_df["2degree.target.keyprotein"]):
        if results_df["2degree.target.keyprotein"].iloc[count] != []:
            genesymbols = []
            try:
                if list(row) != None:
                    for element in list(row):
                        genesymbols += [genesymbol_dict[element]]
                    results_df["2degree.target.keyprotein"].iloc[count] = [genesymbols]
            except TypeError:
                print("Oops! Error while translating")

    # save all results
    results_df.to_csv(address_output_file, index=False)
     
    print('Results:')
    print(results_df)
    print("{} times training Runtime: {:2f} Minutes".format(disease_name, ((time.time() - start_time) / 60)))
    print('Program finished.')

if __name__ == '__main__':
    main()
