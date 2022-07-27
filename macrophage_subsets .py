#Import data
import scvelo as scv
scv.logging.print_version()
import numpy as np
scv.set_figure_params('scvelo') # for beautified visualization
import os
adata=scv.read("loom files/120a_velocyted.loom")
adata2=scv.read("loom files/l20b_velocyted.loom")
adata3=scv.read("loom files/SHAM-CTRL-CR.loom")
adata.obs.index=list(map(lambda x: "L120a_"+x ,adata.obs.index))
adata.var_names_make_unique()
adata2.obs.index=list(map(lambda x: "L120b_"+x,adata2.obs.index))
adata2.var_names_make_unique()

adata3.obs.index=list(map(lambda x: "SHAM-CTRL_"+ x.split('-CR:')[1].split('x')[0]+'-1',adata3.obs.index))
adata3.var_names_make_unique()
adata=adata.concatenate(adata2,adata3,index_unique=None)

del adata2
del adata3
adata.var_names_make_unique()

print(adata.obs_names[-10:].tolist())
#test
a = 'SHAM-CTRL-CR:TTGGTTTCAAGCGGATx'
"SHAM-CTRL_"+x.split('-CR:')[1].split('x')[0]+'-1'

## Integrate with Seurat   -whole dataset
import pandas as pd
umap=pd.read_csv("three_data_UMAP_10X.csv",index_col=0)
adata=adata[umap.index,]

scv.pl.proportions(adata)

adata

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=1000)

scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)

adata.obsm["X_umap"]=umap.values

scv.set_figure_params('scvelo') 

scv.pl.velocity_embedding_stream(adata, basis='umap')

## Interest sub-set
import pandas as pd
umap=pd.read_csv("rmc_UMAP_10X_1.csv",index_col=0)
adata_mc=adata[umap.index,]
adata_mc.obs["leiden"]=pd.read_csv("rmc_clusterinfo_10X_1.csv")["seurat_clusters"].values.astype("str")


### Import loom file (of the subset) generated using as.loom by Seurat

import scvelo as scv
adata_mc =scv.read("rmc1.loom")
adata_mc

### Merge/Filter interest cells

sceMerge = scv.utils.merge(adata,adata_mc)

sceMerge

scv.pp.filter_and_normalize(sceMerge, min_shared_counts=20, n_top_genes=1000)

scv.pp.moments(sceMerge, n_pcs=30, n_neighbors=30)

scv.tl.recover_dynamics(sceMerge)
scv.tl.velocity(sceMerge, mode='dynamical',groupby='batch')

sceMerge.obsm["X_umap"]=umap.values

scv.set_figure_params('scvelo') 

###0:blue #145DA0, 1: orange #FFA384,2: red #EF7C8E;3: cyan #85D2D0; 4:green #81B622;5: yellow #CAA23D; 6: brown #84444C ; 
scv.pl.velocity_embedding_stream(sceMerge, basis='umap',color="seurat_clusters",figsize=(7,7),palette=["#F85C70",'#B68D40','#59981A','#21B6A8',"#04D4F0","#A16AE8","#FA26A0"],title='Macrophages')


print(sceMerge.var['velocity_genes'].sum(), sceMerge.n_vars)
top_genes = sceMerge.var_names[sceMerge.var.fit_likelihood.argsort()[::-1]]
scv.pl.scatter(sceMerge, basis=top_genes[:10], ncols=5)

scv.pl.scatter(sceMerge, color='root_cells',figsize=(7,7))
#scv.pl.scatter(adata, basis=[top_genes])


