#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Gehad Youssef  based on Han et al., 2021

Running time: depends on input, m * n = x combinations (start proteins: m, end proteins: n), 
    ca. 4 min for x = 1,800 on MacbookPro
    i.e. roughly x/500 = y minutes 

"""

import networkx as nx
import pickle
import pandas as pd
import itertools
import matplotlib.pyplot as plt
import numpy as np
import mygene

#%%
### make sure the working folder is the src folder (in Spyder top right) ###

""" input addresses """
### ENTER PATH/NAME OF INPUT: TEXT FILE OF START PROTEIN LIST AND END PROTEIN LIST FOR NETWORK RECONSTRUCTION (ensembl protein ID) ###
start_addr = "../data/input_proteins.txt"
end_addr = "../data/input_proteins.txt"
### INSERT PATH/NAME OF INPUT: EDGELIST FROM STRING DATABASE ###
string_addr = "../databases/9606.protein.links.v11.0.noscores.threshold400.txt"
### ENTER PATH/NAME OF INPUT: LIST OF PROTEIN ID TO GENE, DICTIONARY ###
# based on Jensen Lab dictionary which is based on alias file from STRING database
address_dictionary_protein_to_gene= '../databases/list_protein_to_gene.txt'
### ENTER PATH/NAME OF INPUT: COMPARTMENT DATABASES, LIST OF GENES EXPRESSING PROTEINS THAT ARE IN MEMBRANE/NUCLEUS ###
address_compartment_membrane_genes = '../databases/COMPARTMENT_membrane_gene_names.txt'
address_compartment_nucleus_genes = '../databases/COMPARTMENT_nucleus_gene_names.txt'

""" output addresses """
### ENTER PATH/NAME OF OUTPUT: NEW NETWORK WITH PROTEIN IDs AS NODES###
address_network_genes = "../data/network_genes" #dont add '.pkl'
### ENTER PATH/NAME OF OUTPUT: NETWORK WITH GENE NAMES AS NODES (AS OPPOSED TO PROTEIN IDs) ###
address_network_proteins = "../data/network_proteins" #dont add '.pkl'
### ENTER PATH/NAME OF OUTPUT: LIST OF NODES AS GENE/PROTEIN NAMES ###
address_network_nodes_genes = '../data/network_nodes_genes.txt'
address_network_nodes_proteins = '../data/network_nodes_proteins.txt'
### ENTER PATH/NAME OF OUTPUT: LIST OF EDGES AS GENE/PROTEIN NAMES ###
address_network_edges_genes = '../data/network_edges_genes.txt'
address_network_edges_proteins = '../data/network_edges_proteins.txt'
### ENTER PATH/NAME FOR OUTPUT: list of network genes in membrane, list of network genes in nucleus ###
address_network_membrane_genes = '../data/network_membrane_gene_names.txt'
address_network_nucleus_genes = '../data/network_nucleus_gene_names.txt'
### ENTER PATH/NAME FOR OUTPUT: table of network parameters ###
address_centrality_results= "../data/Centrality_RWR_results.csv"
### ENTER PATH/NAME FOR OUTPUT: CIRCULAR PLOT OF NETWORK ###
# this is a very simple plot, just so the network can be looked at after running the code #
address_network_plot_genes = '../visualisation/network_visualisation_genes.tiff'
address_network_plot_proteins = '../visualisation/network_visualisation_proteins.tiff'

#%% FUNCTIONS

# def load_obj(file_addr):
#     with open(file_addr+ '.pkl', 'rb') as f:
#         return pickle.load(f)
    
def save_obj(obj, file_addr ):
    "function to save pickle file"
    with open(file_addr + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)

def make_shortest_path_sources_to_targets(pair, string_G):
    '''
    Function: Find shortest paths on STRING between all start/end protein combinations
    :input: pair of start/end proteins
    :return: edgelist of all shortest paths between start/end protein pair if existent
    ~ takes 4 min for an initial list of 43 proteins, i.e. 43x43 = ~1,850 combinations
    '''

    # find shortest paths between start and end protein
    all_path_list = []
    try:
        # exclude self-loops
        if pair[0] != pair[1]:
            for p in nx.all_shortest_paths(string_G, source=pair[0], target=pair[1]):
                pairs_from_p = pairwise(p)
                all_path_list += pairs_from_p
                print(pair, 'HIT!')
    except:
        print(pair, "No Path")

    # return new edge list
    return all_path_list

def mapping_genes_id_to(prot_list,id_from, id_to, species):
    mg = mygene.MyGeneInfo()
    return mg.querymany(prot_list, scopes=id_from, fields= id_to, species=species, as_dataframe=True, verbose = True)

def calc_network_centrality_RWR(network, start_list, end_list, result_save_dir):
    '''
    Function: calculate eigenvector centrality, degree centrality
    :input: network, membrane list as start list, nucleus list as end list, address
    for output, thresholds for centralities and RWR
    :return: csv file containing centrality and RWR-values for each node,
    txt-files containing above threshold genes
    '''
    network_nodes = network.nodes()

    # eigenvector centrality
    try:
        eigenvector_centrality_dict = nx.eigenvector_centrality(network) #Returns:Dictionary of nodes with eigenvector centrality as the value.
    except:
        eigenvector_centrality_dict = dict()
        for node in network_nodes:
            eigenvector_centrality_dict[node] = 0.0
            
    # degree centrality
    try:
        degree_centrality_dict = nx.degree_centrality(network) #Returns:Dictionary of nodes with degree centrality as the value.
    except:
        degree_centrality_dict = dict()
        for node in network_nodes:
            degree_centrality_dict[node] = 0.0
            
    # betweeness centrality
    try:
        between_centrality_dict = nx.algorithms.betweenness_centrality(network) #Returns: Dictionary of nodes with betweenness centrality as the value.

    except:
        between_centrality_dict = dict()
        for node in network_nodes:
            between_centrality_dict[node] = 0.0
            
    # betweeness centrality subset: starts at membrane bound proteins ends in nuclear proteins
    try:
        betweensub_centrality_dict = nx.algorithms.betweenness_centrality_subset(network, sources=start_list, targets=end_list) #Returns: Dictionary of nodes with betweenness centrality as the value.

    except:
        betweensub_centrality_dict = dict()
        for node in network_nodes:
            betweensub_centrality_dict[node] = 0.0

    # edge betweenness centrality
    edge_BC_no_Zero = dict()
    edge_BC = nx.edge_betweenness_centrality(network) #?

    for k,v in edge_BC.items():
        edge_BC_no_Zero[k] = v+0.1 # add 0.1 to each edge score to make non-zero
    
    # set non-zero edge betweenness centrality as edge weight
    nx.set_edge_attributes(network, edge_BC_no_Zero, 'weight') #set edge weights as non-zero edge betweenness-centrality score

    # pagerank 
    # returns dictionary of nodes with PageRank as value
    PR_score = nx.pagerank(network)
    # personalised pagerank: jumps back to membrane proteins when jumping
    # in COVID-paper: start_genes_for_PR = start_list #{'COVID19':1}, here:
    dict_membrane = dict.fromkeys(start_list,1)
    PRsub_score = nx.pagerank(network, personalization=dict_membrane) 
    
    # export results as csv
    network_property_df = pd.DataFrame(
        columns=['Eigen', 'Degree', 'Between', 'Between Sub', 'RWR', 'RWR Sub'])
    for node in network_nodes:
        network_property_df.loc[node] = [eigenvector_centrality_dict[node], degree_centrality_dict[node],
                                         between_centrality_dict[node], betweensub_centrality_dict[node], 
                                         PR_score[node], PRsub_score[node]]

    network_property_df.to_csv(result_save_dir)
    
def make_circular_plot_network(nodes, G_gene, address_network_plot, start_list):
    '''
    Function: make a circular plot of the network, save as tiff in visualisation folder
    :input: network, list of nodes
    :return: none
    '''
    color_map = []
    for node in nodes:
        if node in start_list:
            color_map.append('lightcoral')
        else:
            color_map.append('lightcyan')   
            
    plt.subplots(figsize=(30,30))
    
    n=len(nodes)
    angle = []
    angle_dict = {}
    for i, node in zip(range(n),nodes):
        theta = 2.0*np.pi*i/n
        angle.append((np.cos(theta),np.sin(theta)))
        angle_dict[node] = theta
    pos = {}
    for node_i, node in enumerate(nodes):
        pos[node] = angle[node_i]
        
    description = nx.draw_networkx_labels(G_gene, pos, font_size=20 )
    for node, t in description.items():
        t.set_rotation(angle_dict[node]*360.0/(2.0*np.pi))
                                          
    nx.draw_circular(G_gene, node_color=color_map, with_labels=False, font_size=18, font_color ='k', edge_color = 'grey')
    plt.savefig(address_network_plot)
    
#%% MAIN
# NETWORK RECONSTRUCTION

print("Start program. Program will tell you when finished.")

# open start / end protein list for network reconstruction
start_df = pd.read_csv(start_addr, sep='\t', header = None)
start_list_prot = start_df.iloc[:,0].to_list()
end_df = pd.read_csv(end_addr, sep='\t', header = None)
end_list_prot = end_df.iloc[:,0].to_list()

print("Create network:")

with open(string_addr) as network_f:
    string_network_edges = [x.strip().split(',') for x in network_f.readlines()]
print('open STRING list done')

# STRING edgelist to network
string_G = nx.Graph(string_network_edges)
print('make STRING network done')

ppi_deg_pair_list = list(itertools.product(start_list_prot, end_list_prot))
print('start/end protein combinations list done')

shortest_paths_result =[]
for i in range(len(ppi_deg_pair_list)):
    shortest_paths_result += make_shortest_path_sources_to_targets(ppi_deg_pair_list[i], string_G)

# make new graph with hidden layers
G_prot = nx.Graph()
G_prot.add_edges_from(shortest_paths_result)

# save new network with protein IDs as nodes
save_obj(G_prot, address_network_proteins)
#%% MAIN
# TRANSLATION PROTEIN TO GENE

print('\nRenaming network nodes (Protein ID to gene symbol).')

# translate  protein IDs using local dictionary from JENSEN lab
# open dictionary (ensembl protein IDs to gene symbols, JENSEN lab)
network_dict = pd.read_csv(address_dictionary_protein_to_gene,  delimiter=',', header=None)
network_dict = dict(sorted(network_dict.values.tolist()))
#rename nodes in network
G_gene = nx.relabel.relabel_nodes(G_prot, mapping=network_dict, copy=True)

if len(G_gene.nodes) == len(G_prot.nodes):
    print('Renamed all nodes successfully.')
else:
    print('Not all nodes renamed, or more than one protein corresponds to the same gene.')

# also tranlsate the original input lists of proteins for later use
start_list_gene = []
for protein in start_list_prot:
    if protein in network_dict:
        start_list_gene += [network_dict[protein]]
    else:
        start_list_gene += [protein]
end_list_gene = []
for protein in end_list_prot:
    if protein in network_dict:
        end_list_gene += [network_dict[protein]]
    else:
        end_list_gene += [protein]

# save list of nodes and network with gene names
network_nodes_genes = pd.DataFrame(G_gene.nodes)
network_nodes_genes.to_csv(address_network_nodes_genes, header=None, index=False)
network_nodes_proteins = pd.DataFrame(G_prot.nodes)
network_nodes_proteins.to_csv(address_network_nodes_proteins, header=None, index=False)
network_edges_genes = pd.DataFrame(G_gene.edges)
network_edges_genes.to_csv(address_network_edges_genes, header=None, index=False)
network_edges_proteins = pd.DataFrame(G_prot.edges)
network_edges_proteins.to_csv(address_network_edges_proteins, header=None, index=False)
save_obj(G_gene, address_network_genes)

#%% MAIN
# PROTEIN LOCALISATION

print('Finding proteins that are localised in nucleus/membrane.')
# open pre-edited lists of genes in membrane and in nucleus (based on COMPARTMENT database on 30th July 2021)
df_comp_mem = pd.read_csv(address_compartment_membrane_genes, delimiter='\t', header=None)
df_comp_nuc = pd.read_csv(address_compartment_nucleus_genes, delimiter='\t', header=None)
list_comp_mem = df_comp_mem.iloc[:,0].to_list()
list_comp_nuc = df_comp_nuc.iloc[:,0].to_list()

# convert network nodes to list
nodes = list(G_gene.nodes)

# make list of genes (=nodes) that are in membrane/nucleus
network_membrane_gene_names = []
network_nucleus_gene_names = []

for node in nodes:
    if node in list_comp_mem:
        network_membrane_gene_names += [node]
    if node in list_comp_nuc:
        network_nucleus_gene_names += [node]

# save lists as txt
network_membrane_gene_names_list = pd.DataFrame(list(network_membrane_gene_names))
network_nucleus_gene_names_list = pd.DataFrame(list(network_nucleus_gene_names))
network_membrane_gene_names_list.to_csv(address_network_membrane_genes, header=None, index=False)
network_nucleus_gene_names_list.to_csv(address_network_nucleus_genes, header=None, index=False)

membrane_list = list(network_membrane_gene_names)
nucleus_list = list(network_nucleus_gene_names)

#%% MAIN
# ANALYSE NETWORK

# give lists, network, output address to function
print("Analyzing pathway network using multiple centrality methods.")
calc_network_centrality_RWR(G_gene, membrane_list, nucleus_list, address_centrality_results)

#%% MAIN
# VISUALISE NETWORK

print("Generating circular plot of network.")
make_circular_plot_network(list(G_prot.nodes), G_prot, address_network_plot_proteins, start_list_prot)
make_circular_plot_network(list(G_gene.nodes), G_gene, address_network_plot_genes, start_list_gene)
print("Program finished.")
