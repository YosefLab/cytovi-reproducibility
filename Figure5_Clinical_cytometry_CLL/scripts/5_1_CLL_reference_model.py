import cytovi
import anndata as ad
import scanpy as sc

# read data
adata =  ad.read_h5ad('../data/raw/Flow cytometry/CLL/full_cohort/processed/2024-08-30_CLL_full_cohort_reference_model_100k.h5ad')


encoder_marker_list = [var for var in adata.var_names if var not in ['Kappa', 'Lambda']]

train_kwargs = {
        "n_epochs_kl_warmup": 30}

# train CyTOF models
cytovi.CytoVI.setup_anndata(
    adata,
    layer="scaled",
    batch_key="batch",
)

model = cytovi.CytoVI(adata, n_latent=5, prior_mixture_k=5, encoder_marker_list=encoder_marker_list)
model.train(**train_kwargs)
model.save('/home/labs/amit/floriani/Lab/PROJECTS/FlowVI/models/CLL/2024-09-02_CLL_reference_model', overwrite=True)


# get latent space
adata.obsm["X_CytoVI"] = model.get_latent_representation()

# use FlowVI latent space for UMAP generation
sc.pp.neighbors(adata, use_rep="X_CytoVI")
sc.tl.umap(adata)

adata.obsm["X_umap_CytoVI"] = adata.obsm["X_umap"]

# save anndata to disk
adata.write('../data/raw/Flow cytometry/CLL/full_cohort/processed/2024-09-02_CLL_full_cohort_reference_model_100k_latent.h5ad')