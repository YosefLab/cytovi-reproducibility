import cytovi
import anndata as ad

# read data
adata = ad.read_h5ad('../data/raw/CyTOF/twins/surface/batch1/annotated/2024-05-22_norm_samples_batch1_ann_100k.h5ad')

rm_vars = ['CD138']
adata = adata[:, [var for var in adata.var_names if var not in rm_vars]].copy()

# train CyTOF models
cytovi.CytoVI.setup_anndata(
    adata,
    layer="scaled"
)

model = cytovi.CytoVI(adata, protein_likelihood = 'normal', n_latent=10)
model.train(n_epochs_kl_warmup = 200)
model.save('/home/labs/amit/floriani/Lab/PROJECTS/FlowVI/models/dist_eval/2024-10-11_CyTOF_normal_scaled', overwrite=True)

# train CyTOF models
cytovi.CytoVI.setup_anndata(
    adata,
    layer="scaled",
)

model = cytovi.CytoVI(adata, protein_likelihood = 'beta', n_latent=10)
model.train(n_epochs_kl_warmup = 200)
model.save('/home/labs/amit/floriani/Lab/PROJECTS/FlowVI/models/dist_eval/2024-10-11_CyTOF_beta', overwrite=True)


# train CyTOF models
cytovi.CytoVI.setup_anndata(
    adata,
    layer="transformed",
)

model = cytovi.CytoVI(adata, protein_likelihood = 'normal', n_latent=10)
model.train(n_epochs_kl_warmup = 200)
model.save('/home/labs/amit/floriani/Lab/PROJECTS/FlowVI/models/dist_eval/2024-10-11_CyTOF_normal', overwrite=True)