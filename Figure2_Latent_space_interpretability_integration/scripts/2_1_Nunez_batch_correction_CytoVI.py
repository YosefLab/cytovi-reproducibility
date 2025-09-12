import cytovi
import anndata as ad
import scanpy as sc

# read data
adata = ad.read_h5ad('../data/2024-01-16_model_eval_multi_batch/2024-06-14_aurora_batch_figure_scaling_variants.h5ad')

model_path = '/home/labs/amit/floriani/Lab/PROJECTS/FlowVI/models/dist_eval/2024-06-14_aurora_batch_figure_new_hyper_'

train_kwargs = {
        "n_epochs_kl_warmup": 250
        }

model_kwargs = {
        'protein_likelihood': 'normal',
        'n_latent': 20
        }


# train model with min max scaling
model_key = 'minmax'

cytovi.CytoVI.setup_anndata(
    adata,
    layer="scaled",
    batch_key='batch'
)

model = cytovi.CytoVI(adata, **model_kwargs)
model.train(**train_kwargs)
model.save(f'{model_path}{model_key}', overwrite=True)

# compute latent and umap
adata.obsm[f"X_CytoVI_{model_key}"] = model.get_latent_representation()
sc.pp.neighbors(adata, use_rep=f"X_CytoVI_{model_key}", key_added= f'neighbors_CytoVI_{model_key}')
sc.tl.umap(adata, neighbors_key = f'neighbors_CytoVI_{model_key}')
adata.obsm[f"X_umap_CytoVI_{model_key}"] = adata.obsm['X_umap'].copy()


# train model with z-score scaling
model_key = 'zscore'

cytovi.CytoVI.setup_anndata(
    adata,
    layer="scaled_z",
    batch_key='batch'
)

model = cytovi.CytoVI(adata, **model_kwargs)
model.train(**train_kwargs)
model.save(f'{model_path}{model_key}', overwrite=True)

# compute latent and umap
adata.obsm[f"X_CytoVI_{model_key}"] = model.get_latent_representation()
sc.pp.neighbors(adata, use_rep=f"X_CytoVI_{model_key}", key_added= f'neighbors_CytoVI_{model_key}')
sc.tl.umap(adata, neighbors_key = f'neighbors_CytoVI_{model_key}')
adata.obsm[f"X_umap_CytoVI_{model_key}"] = adata.obsm['X_umap'].copy()


# train model with rank scaling
model_key = 'rank'

cytovi.CytoVI.setup_anndata(
    adata,
    layer="scaled_rank",
    batch_key='batch'
)

model = cytovi.CytoVI(adata, **model_kwargs)
model.train(**train_kwargs)
model.save(f'{model_path}{model_key}', overwrite=True)

# compute latent and umap
adata.obsm[f"X_CytoVI_{model_key}"] = model.get_latent_representation()
sc.pp.neighbors(adata, use_rep=f"X_CytoVI_{model_key}", key_added= f'neighbors_CytoVI_{model_key}')
sc.tl.umap(adata, neighbors_key = f'neighbors_CytoVI_{model_key}')
adata.obsm[f"X_umap_CytoVI_{model_key}"] = adata.obsm['X_umap'].copy()


adata.write('/home/labs/amit/floriani/Lab/PROJECTS/FlowVI/data/2024-01-16_model_eval_multi_batch/2024-07-10_aurora_batch_new_hyper_latents_umap.h5ad')