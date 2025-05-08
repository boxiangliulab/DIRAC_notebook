import os
import random

import pandas as pd
import numpy as np
import torch
import scanpy as sc
import matplotlib.pyplot as plt
import anndata
import sklearn
import scipy
import yaml
import scvi

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

data_path = "/home/project/11003054/changxu/Projects/DIRAC/Section-4"
methods = "SCANVI"
sim_methods = "layers"
omics = "RNA"
colormaps= ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#FFFF33", "#F781BF", "#999999",
            "#E5D8BD", "#B3CDE3", "#CCEBC5", "#FED9A6", "#FBB4AE", "#8DD3C7", "#BEBADA", "#80B1D3", "#B3DE69", "#FCCDE5",
            "#BC80BD", "#FFED6F", "#8DA0CB", "#E78AC3", "#E5C494", "#CCCCCC", "#FB9A99", "#E31A1C", "#CAB2D6","#6A3D9A", 
            "#B15928"]
cell_type_col = 'cell.type'
n_per_class = 100

results_list = []
for i in range(8):
    adata = anndata.read_h5ad(os.path.join(data_path, "scMultiSim_data", f"sim-{sim_methods}", f"sim_{omics}_{i}.h5ad"))
    sc.pp.filter_genes(adata, min_cells=3)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.scale(adata)
    adata.obs['celltype_num'] = pd.Categorical(adata.obs[cell_type_col]).codes
    pairs = dict(enumerate(pd.Categorical(adata.obs[cell_type_col]).categories))
    
    for j in range(1, 6):
        save_path = os.path.join(data_path, "Results", f"merged_{omics}_{methods}_{sim_methods}_{i}_{j}")
        if not os.path.exists(save_path):
            os.makedirs(save_path)

        colormaps_clusters = dict(set(zip(adata.obs["cell.type"].unique(), colormaps)))
        
        labels = np.repeat("Unknown", adata.shape[0])
        labels = labels.astype("<U43")
        for x in np.unique(adata.obs[cell_type_col]):
            idx = np.where((adata.obs[cell_type_col] == x)& (adata.obs["batch"] != j))[0]
            sampled = np.random.choice(idx, np.min([n_per_class, len(idx)]))
            labels[sampled] = adata.obs[cell_type_col][sampled]
        adata.obs["celltype_scanvi"] = labels
        scvi.model.SCVI.setup_anndata(
            adata, 
            layer="counts",
            batch_key='batch', 
            labels_key="celltype_scanvi", 
        )
        scvi_model = scvi.model.SCVI(
            adata,
            n_latent=30,
            n_layers=2,
        )
        scvi_model.train(100)
        scanvi_model = scvi.model.SCANVI.from_scvi_model(scvi_model, "Unknown")
        scanvi_model.train(25)

        SCANVI_LATENT_KEY = "X_scANVI"
        adata.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation(adata)
        adata.obs[f"{methods}"] = scanvi_model.predict(adata)
        target_adata = adata[adata.obs["batch"] == j]
        
        metrics_all = {
                "Accuracy Score":
                    float(sklearn.metrics.accuracy_score(target_adata.obs[cell_type_col], target_adata.obs[f"{methods}"])),
                "Precision Score":
                    float(sklearn.metrics.precision_score(target_adata.obs[cell_type_col], target_adata.obs[f"{methods}"], average='weighted')),
                "Recall Score":
                    float(sklearn.metrics.recall_score(target_adata.obs[cell_type_col], target_adata.obs[f"{methods}"], average='weighted')),
                "F1 Score":
                    float(sklearn.metrics.f1_score(target_adata.obs[cell_type_col], target_adata.obs[f"{methods}"], average='weighted'))}
        print(metrics_all)
        with open(os.path.join(save_path, "annotate_mertics.yaml"), "w", encoding="utf-8") as g:
            yaml.dump(metrics_all, g, default_flow_style=False, allow_unicode=True)
        results_list.append(metrics_all)
        sc.pl.embedding(target_adata, basis='spatial', color=['cell.type', f'{methods}'], title=['Ground truth',f"{methods} ACC {metrics_all['Accuracy Score']}"], s=30, show=False, palette=colormaps_clusters)
        plt.savefig(os.path.join(save_path, f"{omics}_{methods}_{sim_methods}.pdf"), dpi=300)
        plt.show()        
        target_adata.obs.to_csv(os.path.join(save_path, f"{omics}_{methods}_{sim_methods}.csv"))
        

summary_df = pd.DataFrame(results_list)
summary_df.to_csv(os.path.join(data_path, "Results", f"sim_metrics_{omics}_{methods}_{sim_methods}.csv"))