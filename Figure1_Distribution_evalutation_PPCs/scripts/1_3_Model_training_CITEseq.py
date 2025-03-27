import cytovi
import anndata as ad

# read data
adata = ad.read_h5ad('../data/raw/CITE_seq/Seurat_PBMC/2024-05-22_cite_seq_processed_P7_only.h5ad')

# train CyTOF models
cytovi.CytoVI.setup_anndata(
    adata,
    layer="scaled"
)

model = cytovi.CytoVI(adata, protein_likelihood = 'normal', n_latent=20)
model.train(n_epochs_kl_warmup = 200)
model.save('/home/labs/amit/floriani/Lab/PROJECTS/FlowVI/models/dist_eval/2024-07-10_CITE_seq_normal_scaled', overwrite=True)

# train CyTOF models
cytovi.CytoVI.setup_anndata(
    adata,
    layer="scaled",
)

model = cytovi.CytoVI(adata, protein_likelihood = 'beta', n_latent=20)
model.train(n_epochs_kl_warmup = 200)
model.save('/home/labs/amit/floriani/Lab/PROJECTS/FlowVI/models/dist_eval/2024-07-10_CITE_seq_beta', overwrite=True)


# train CyTOF models
cytovi.CytoVI.setup_anndata(
    adata,
    layer="transformed",
)

model = cytovi.CytoVI(adata, protein_likelihood = 'normal', n_latent=20)
model.train(n_epochs_kl_warmup = 200)
model.save('/home/labs/amit/floriani/Lab/PROJECTS/FlowVI/models/dist_eval/2024-07-10_CITE_seq_normal', overwrite=True)