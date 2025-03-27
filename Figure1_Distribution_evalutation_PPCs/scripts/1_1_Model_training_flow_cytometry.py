import cytovi
import anndata as ad

# read data
adata = ad.read_h5ad('../data/raw/Spectral flow/Nunez/For Chiquito/annotated/2024-01-24_norm_sample_batch1_ann.h5ad')

# train aurora models
cytovi.CytoVI.setup_anndata(
    adata,
    layer="scaled"
)

model = cytovi.CytoVI(adata, protein_likelihood = 'normal', n_latent=20)
model.train(n_epochs_kl_warmup=200)
model.save('/home/labs/amit/floriani/Lab/PROJECTS/FlowVI/models/dist_eval/2024-07-10_aurora_normal_scaled', overwrite=True)

# train aurora models
cytovi.CytoVI.setup_anndata(
    adata,
    layer="scaled",
)

model = cytovi.CytoVI(adata, protein_likelihood = 'beta', n_latent=20)
model.train(n_epochs_kl_warmup=200)
model.save('/home/labs/amit/floriani/Lab/PROJECTS/FlowVI/models/dist_eval/2024-07-10_aurora_beta', overwrite=True)


# train aurora models
cytovi.pp.arcsinh(adata, global_scaling_factor=2000, transform_scatter=True)

cytovi.CytoVI.setup_anndata(
    adata,
    layer="transformed",
    batch_key='batch'
)

model = cytovi.CytoVI(adata, protein_likelihood = 'normal', n_latent=20)
model.train(n_epochs_kl_warmup=200)
model.save('/home/labs/amit/floriani/Lab/PROJECTS/FlowVI/models/dist_eval/2024-07-10_aurora_normal', overwrite=True)