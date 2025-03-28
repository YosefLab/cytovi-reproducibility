import cytovi
import anndata as ad
import scanpy as sc

# read data
flow_data_dir = '/home/labs/amit/floriani/Lab/PROJECTS/FlowVI/data/raw/Flow cytometry/BNHL/'
adata = ad.read_h5ad(flow_data_dir + '2024-11-19_BNHL_flow_both_panels_combined.h5ad')
adata = cytovi.pp.subsample(adata, n_obs_group=10000, groupby='PatientID')

# rename Entity categories
adata.obs['Entity'] = adata.obs['Entity'].replace({'DLBCL, non-GCB': 'DLBCL', 'DLBCL, GCB': 'DLBCL'}).astype('category')

# save to disk
adata.write_h5ad(flow_data_dir + '2024-11-26_BNHL_flow_both_panels_combined_10k_per_patient.h5ad')

# train flow model
cytovi.CytoVI.setup_anndata(
    adata,
    layer="scaled",
    batch_key='batch',
    sample_key='PatientID'
)

model = cytovi.CytoVI(adata)
model.train(n_epochs_kl_warmup = 50, max_epochs = 800, batch_size = 1024)
model.save('/home/labs/amit/floriani/Lab/PROJECTS/FlowVI/models/BNHL/2024-11-26_flow_imputation_10K_per_patient', overwrite=True)

# get latent space
adata.obsm['X_CytoVI'] = model.get_latent_representation()

# impute missing markers
adata.layers['imputed'] = model.get_normalized_expression(n_samples = 10)

# save adata
adata.write_h5ad(flow_data_dir + '2024-11-26_BNHL_flow_both_panels_combined_10k_per_patient_integrated_imputed.h5ad')

# run DA analysis
print('Running DA analysis')
log_probs = model.differential_abundance(batch_size = 10000)

# save DA results
log_probs.to_csv(flow_data_dir + '2024-11-26_BNHL_flow_both_panels_combined_10k_per_patient_DA_log_probs.csv')

# save adata
adata.write_h5ad(flow_data_dir + '2024-11-26_BNHL_flow_both_panels_combined_10k_per_patient_integrated_imputed_DA_res.h5ad')