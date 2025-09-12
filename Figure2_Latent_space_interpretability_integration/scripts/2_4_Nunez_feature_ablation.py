import cytovi
import anndata as ad
import random


model_path = '/home/projects/amit/floriani/Lab/PROJECTS/FlowVI/models/marker_ablation/2025-08-15_marker_ablation'

adata = ad.read_h5ad('/home/projects/amit/floriani/Lab/PROJECTS/FlowVI/data/2024-01-16_model_eval_multi_batch/2024-01-26_nunez_norm.h5ad')
adata = cytovi.pp.subsample(adata, n_obs= 20000, groupby = 'batch')


markers = adata.var_names.tolist()
n_markers = len(markers)

train_kwargs = {
        "n_epochs_kl_warmup": 50
        }

for i in range(n_markers -1):
    masked_markers = random.sample(markers, i+1)
    print(f"Masked markers: {masked_markers}")

    adata_temp = cytovi.pp.mask_markers(adata, markers = masked_markers)

    cytovi.CytoVI.setup_anndata(
        adata_temp,
        layer="scaled",
        batch_key='batch'
    )

    model = cytovi.CytoVI(adata_temp)
    model.train(**train_kwargs)
    model.save(f'{model_path}Masking_it_{i+1}', overwrite=True)

    adata.uns[f'masking_it_{i+1}'] = masked_markers
    adata.obsm[f'X_masking_it_{i+1}'] = model.get_latent_representation()
    adata.write(f'{model_path}/2025-08-15_aurora_batch_marker_masking_latents.h5ad')