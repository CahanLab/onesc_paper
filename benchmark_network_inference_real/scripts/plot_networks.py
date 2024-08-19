import graphviz
import pandas as pd
import pylab
import networkx as nx
import os
import matplotlib.pyplot as plt
plt.style.use("default")

input_path = '../output/curated_networks/'
output_path = '../output/plot_networks/'

def create_graphViz(network_df, name):
    network_df.index = list(range(0, network_df.shape[0]))
    g = graphviz.Digraph('G', filename=name+".gv", graph_attr={'fontsize': '100'})
    for temp_index in network_df.index:
        if network_df.loc[temp_index, 'Type'] == "+":
            g.edge(network_df.loc[temp_index, "TF"], network_df.loc[temp_index, "TG"])
        else:
            g.edge(network_df.loc[temp_index, "TF"], network_df.loc[temp_index, "TG"], arrowhead = 'tee')
    return g

def create_networkx(network_df):
    G=nx.DiGraph()
    network_df.index = list(range(0, network_df.shape[0]))
    for temp_index in network_df.index:
        if network_df.loc[temp_index, 'Type'] == "+":
            G.add_edge(network_df.loc[temp_index, "TF"], network_df.loc[temp_index, "TG"], label='activation', color='#f5428d')
        else:
            G.add_edge(network_df.loc[temp_index, "TF"], network_df.loc[temp_index, "TG"], label='activation', color='#4c36f5')
    return G

def plot_networkx(G, save_output_path):
    edges = G.edges()
    colors = [G[u][v]['color'] for u,v in edges]
    #nx.draw(G, edge_color=colors, with_labels=True, font_weight='bold')

    #legend stuff
    pos=nx.circular_layout(G)
    h1 = nx.draw_networkx_nodes(G, pos=pos, node_color = 'w',
                            alpha = 0.9, node_size = 50, linewidths=1)
    h2 = nx.draw_networkx_edges(G, pos=pos, width=2, edge_color=colors, arrowstyle="->")
    h3 = nx.draw_networkx_labels(G, pos=pos, font_size=12, font_color='#000000',font_weight='bold')
    plt.savefig(os.path.join(save_output_path, "network.png"), dpi=150)
    plt.clf()

for dataset in os.listdir(input_path):
    if dataset == '.DS_Store':
        continue
    for method in os.listdir(os.path.join(input_path, dataset)):
        if method == '.DS_Store':
            continue
        curated_network = pd.read_csv(os.path.join(input_path, dataset, method, 'curated_network_density.csv'), index_col = 0)
        if os.path.isdir(os.path.join(output_path, dataset, method)) == False: 
            os.makedirs(os.path.join(output_path, dataset, method))
        G = create_networkx(curated_network)
        plot_networkx(G, os.path.join(output_path, dataset, method))