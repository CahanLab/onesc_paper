import numpy as np 
import pandas as pd
import scanpy as sc
import os 
import matplotlib.pyplot as plt 
import pickle
plt.style.use('default')
sc.set_figure_params(dpi_save = 150)

##### plotting out the data for myeloids #####
if os.path.isdir(os.path.join("../output/myeloid_progenitors_full/figures")) == False:
    os.makedirs(os.path.join("../output/myeloid_progenitors_full/figures"))
sc.settings.figdir = "../output/myeloid_progenitors_full/figures/"

adata = sc.read_h5ad("../input/myeloid_progenitors/filtered_adata.h5ad")
sc.pl.umap(adata, color = 'cell_types', palette = {'CMP': '#66c2a5', 'GMP': '#8da0cb', 'Monocytes': '#e5c494', 'Granulocytes': '#e78ac3', 'MEP': '#a6d854', 'Erythrocytes': '#fc8d62', 'MK': "#ffd92f"}, save = 'ct.png')

sc.pl.umap(adata, color = 'dpt_pseudotime', save = 'pt.png')

exp_df = pd.read_csv("../output/myeloid_progenitors_full/train_exp.csv", index_col = 0)

sc.pl.dotplot(adata, var_names = exp_df.index, groupby = 'cell_types', categories_order = ['CMP', 'MEP', 'GMP', 'MK', 'Erythrocytes', 'Granulocytes', 'Monocytes'], save = '_network_genes.png')

for temp_gene in exp_df.index:
    print(temp_gene)
    sc.pl.umap(adata, color = temp_gene, save = "_" + temp_gene + '.png')
