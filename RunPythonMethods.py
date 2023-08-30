## Imports
def warn(*args, **kwargs):
    pass

import warnings
warnings.warn = warn
import os
import sys
sys.stderr = open('tmp.err.log', 'w')
sys.stdout = open('tmp.out.log', 'w')
import scvi
import scanpy as sc
import anndata as ad
import scanorama
from scipy import sparse
from sklearn.preprocessing import LabelEncoder
import numpy as np
import pandas as pd

## Import and setup data
adata = ad.read_h5ad("tmp.h5ad")
adata.layers["raw"] = sparse.csr_matrix(adata.raw.X.copy())
adata.X = sparse.csr_matrix(adata.raw.X.copy())

## Create batch indices
le = LabelEncoder()
adata.obs['batch_indices'] = le.fit_transform(adata.obs["batch"].values)
adata.obs['batch_indices'] = adata.obs['batch_indices'].astype("category")

## Normalize and log transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

## Import variable genes
f = open("tmp.features", "r")
hvg = f.read().splitlines()
f.close()

## Process variable genes
adata.var['highly_variable'] = [True if g in hvg else False for g in adata.var_names]

## Subset to variable genes
adata = adata[:, adata.var.highly_variable]

### Integrate using scanorama
## Split into batches
split = []
batch_categories = adata.obs["batch_indices"].cat.categories
for i in batch_categories:
	split.append(adata[adata.obs["batch_indices"] == i].copy())

## Batch-aware scaling
for i in split:
	sc.pp.scale(i)

## Integrate
for rep in range(5):
	corrected = scanorama.integrate_scanpy(split, verbose = 0)
	scan = ad.concat(split)
	arr = scan.obsm["X_scanorama"]
	df = pd.DataFrame(arr, index=scan.obs.index)
	file_out = "tmp_Scanorama_rep" + str(rep) + ".txt"
	df.to_csv(file_out)

### Integrate using scVI
## Copy the adata object
adata = adata.copy()

## Setup early stopping
## Copy from here: https://discourse.scverse.org/t/mapping-old-scanvi-code-to-0-16-2/494
early_stopping_kwargs = {
    "early_stopping": True, 
    "early_stopping_monitor": 'elbo_validation',
    "early_stopping_patience": 10,
    "early_stopping_min_delta": 0.0,
}

## Loop across replicates
for rep in range(5):
  ## Setup the model
  scvi.model.SCVI.setup_anndata(adata, layer="raw", batch_key="batch_indices")

  ## Create the model
  vae = scvi.model.SCVI(adata)

  ## Train the model
  vae.train(max_epochs = 500, **early_stopping_kwargs, check_val_every_n_epoch=1)

  ## Insert latent space
  latent = vae.get_latent_representation()
  df = pd.DataFrame(latent, index=adata.obs.index)
  
  ## Write objects
  file_out = "tmp_scVI_rep" + str(rep) + ".txt"
  df.to_csv(file_out)


