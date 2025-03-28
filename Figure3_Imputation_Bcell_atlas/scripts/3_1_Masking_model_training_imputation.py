import cytovi
import anndata as ad
import numpy as np
import os

# read data
os.chdir('/home/labs/amit/floriani/Lab/PROJECTS/FlowVI/notebooks/')
adata = ad.read_h5ad('../data/2024-02-20_hyperopt_flow_batch/2024-08-20_nunez_norm_50k_synth_batch.h5ad')
model_path = '/home/labs/amit/floriani/Lab/PROJECTS/FlowVI/models/imputation_eval/2024-08-20_flow_50k_synth_batch_impute_'


rm_vars = ['FSC-H', 'SSC-H', 'SSC-B-A', 'SSC-B-H', 'CXCR3', 'PD1']

adata = adata[:, [var for var in adata.var_names if var not in rm_vars]].copy()

train_kwargs = {
        "max_epochs": 1000,
        "n_epochs_kl_warmup": 100
        }

model_kwargs = {
        'n_latent': 20,
        'prior_mixture': True,
        'prior_mixture_k': 20,
        'encode_backbone_only': True
        }


imp_arrays = []

for marker in adata.var_names:
    print(f'Processing {marker} now.')
    adata_masked = cytovi.pp.mask_markers(adata, markers = marker, batch_key = 'synth_batch', masked_batch = 1)

    cytovi.CytoVI.setup_anndata(
                adata_masked,
                layer="scaled",
                batch_key='synth_batch'
    )
    model = cytovi.CytoVI(adata_masked, **model_kwargs)
    model.train(**train_kwargs)
    model.save(f'{model_path}{marker}', overwrite=True)

    # Get the imputed values (expression values) for the specific marker
    imp_array = model.get_normalized_expression(adata_masked, n_samples=50, protein_list=[marker], return_mean=False, transform_batch=None)
    
    # Select cells of interest, filtering by the batch if necessary
    imp_array_eval = imp_array[:, adata.obs['synth_batch'] == 1, :]
    
    # Append the array for this marker
    imp_arrays.append(imp_array_eval)

# Stack the list into a numpy array of shape (n_markers, n_samples, n_cells, n_features)
imp_array_stack = np.stack(imp_arrays, axis=0)

imp_array_stack = np.squeeze(imp_array_stack)

# save
np.save(f'{model_path}imputed_array.np', imp_array_stack)