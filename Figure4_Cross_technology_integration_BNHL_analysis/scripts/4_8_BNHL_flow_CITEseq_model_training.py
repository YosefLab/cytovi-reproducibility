import cytovi
import anndata as ad
import scanpy as sc
import numpy as np

# read data
flow_data_dir = '/home/projects/amit/floriani/Lab/PROJECTS/FlowVI/data/raw/Flow cytometry/BNHL/'
adata = ad.read_h5ad(flow_data_dir + '2025-01-28_BNHL_flow_citeseq_combined.h5ad')

# read RNA
citeseq_data_dir = '/home/projects/amit/floriani/Lab/PROJECTS/FlowVI/data/raw/CITE_seq/BNHL/'
adata_rna = ad.read_h5ad(citeseq_data_dir + '2025-01-17_BNHL_CITEseq_combined_RNA_TotalVI_imputed.h5ad')

# select markers to encode
bb_markers = list(adata.var_names[~np.any(adata.layers['_nan_mask'] == 0, axis=0)])
marker = ['CD3']
enc_marker = [var for var in bb_markers if var not in marker]

# train flow model
cytovi.CytoVI.setup_anndata(
        adata,
        layer="imputed_std",
        batch_key='technology',
        sample_key='PatientID',
        labels_key='harmonized_cell_type'
    )

model = cytovi.CytoVI(adata, encoder_marker_list = enc_marker)
model.train(n_epochs_kl_warmup = 5, max_epochs = 300, batch_size = 2000, plan_kwargs={"scale_adversarial_loss": 4})

model.save('/home/projects/amit/floriani/Lab/PROJECTS/FlowVI/models/BNHL/2025-01-28_BNHL_flow_citeseq_cross', overwrite=True)

# get latent space
adata.obsm['X_CytoVI'] = model.get_latent_representation()

# save adata
adata.write_h5ad(flow_data_dir + '2025-01-28_BNHL_flow_citeseq_combined_integrated.h5ad')

# run DA analysis
print('Running RNA imputation')
adata_imputed_rna = model.impute_rna_from_reference(reference_batch='CITE_seq', adata_rna = adata_rna, layer_key='denoised_rna', use_rep = 'X_CytoVI', return_query_only = True)

# save adata
adata_imputed_rna.write_h5ad(flow_data_dir + '2025-01-28_BNHL_flow_cross_RNA_imputed.h5ad')