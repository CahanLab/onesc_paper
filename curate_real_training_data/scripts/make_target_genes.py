import pandas as pd
import numpy as np 
import scanpy as sc 
import os

adata = sc.read_h5ad("../input/myeloid_progenitors/filtered_adata.h5ad")
adata.uns['log1p']["base"] = None

sc.tl.rank_genes_groups(adata, groupby="cell_types", use_raw=True)

big_df = pd.DataFrame()
for ct in np.unique(adata.obs['cell_types']):
    print(ct)
    dedf = sc.get.rank_genes_groups_df(adata, group=ct)
    dedf['cell_types'] = ct
    big_df = pd.concat([big_df, dedf], axis = 0)

big_df.to_csv("../output/myeloid_progenitors_full/DE_genes.csv")