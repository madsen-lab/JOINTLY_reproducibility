### Imports
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scanpy.external as sce
import scarches as sca
from scipy import sparse
from sklearn.preprocessing import LabelEncoder

### Preprocessing
## Import and setup data
adata_raw = ad.read_h5ad("data/Adipose/watlas_filtered.h5ad")
adata_raw.X = sparse.csr_matrix(adata_raw.raw.X.copy())

## Create batch indices
le = LabelEncoder()
adata_raw.obs['batch_indices'] = le.fit_transform(adata_raw.obs["sample"].values)
adata_raw.obs['batch_indices'] = adata_raw.obs['batch_indices'].astype("category")

## Convert labels to categorical
adata_raw.obs['labels.l2'] = adata_raw.obs['labels.l2'].astype("category")

## Copy the raw object
adata = adata_raw.copy()

## Save a copy of the raw counts
adata.layers["raw"] = sparse.csr_matrix(adata.raw.X.copy())

## Normalize data
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

## Find HVGs
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='batch_indices')

## Subset the raw object to only HVGs
adata_hvg = adata_raw[:, adata.var.highly_variable]
adata_hvg = adata_hvg.copy()

### Train models
## Train an scVI model
sca.models.SCVI.setup_anndata(adata_hvg, batch_key='batch_indices', labels_key='labels.l2')
vae = sca.models.SCVI(adata_hvg, n_layers=2, encode_covariates=True, deeply_inject_covariates=False, use_layer_norm="both", use_batch_norm="none", n_latent = 15, gene_likelihood="nb")
vae.train()

## Train an scANVI model
scanvae = sca.models.SCANVI.from_scvi_model(vae, unlabeled_category = "Unknown")
scanvae.train(max_epochs=20, n_samples_per_label=100)

## Save the latent representation to a text file
latent = scanvae.get_latent_representation()
latent_df = pd.DataFrame(latent)
latent_df.to_csv('data/Adipose/latent/latent_final.csv',header = False, index= False)

## Save each model
scanvae.save('data/Adipose/models/scANVI', overwrite=True)