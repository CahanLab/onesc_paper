import graphviz
import pandas as pd
import pylab
import networkx as nx
import os
import matplotlib.pyplot as plt
plt.style.use("default")
from matplotlib.lines import Line2D


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
            G.add_edge(network_df.loc[temp_index, "TF"], network_df.loc[temp_index, "TG"], label='activation', color='#f768a4')
        else:
            G.add_edge(network_df.loc[temp_index, "TF"], network_df.loc[temp_index, "TG"], label='repression', color='#4c36f5')
    return G

def plot_networkx(G, save_output_path):
    edges = G.edges()
    colors = [G[u][v]['color'] for u,v in edges]
    #nx.draw(G, edge_color=colors, with_labels=True, font_weight='bold')

    #legend stuff
    pos=nx.circular_layout(G)
    h1 = nx.draw_networkx_nodes(G, pos=pos, node_color = 'w',
                            alpha = 0.9, node_size = 1900, linewidths=1, edgecolors='black')
   
    h2 = nx.draw_networkx_edges(G, pos=pos, width=2, edge_color=colors, arrowstyle="->", arrowsize=15, min_source_margin=30,  
                                min_target_margin=30, connectionstyle='arc3, rad = 0.1')
    h3 = nx.draw_networkx_labels(G, pos=pos, font_size=12, font_color='#000000',font_weight='bold')

    # Create custom legend handles
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='w', markersize=5, label='Regulation Type'),  # Example node
        Line2D([0], [0], color='#f768a4', lw=2, label='Activation'),  # Example edge type 1
        Line2D([0], [0], color='#4c36f5', lw=2, label='Repression'),  # Example edge type 2
        # Add more legend elements as needed for different colors or styles
    ]

    # Adding the legend
    plt.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1, 1))
    plt.subplots_adjust(right=1.2, top = 1.2)

    plt.savefig(os.path.join(save_output_path, "network.png"), bbox_inches='tight', dpi=300)
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