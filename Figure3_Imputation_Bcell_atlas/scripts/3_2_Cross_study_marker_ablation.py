import cytovi
import anndata as ad

adata = ad.read_h5ad('../data/2024-06-14_imputation_eval_batch/2024-08-23_kreutmair_nunez_annotated_imputed.h5ad')
adata_sub = cytovi.pp.subsample(adata, n_obs = 20000, groupby='batch')

model_path = '/home/projects/amit/floriani/Lab/PROJECTS/FlowVI/models/imputation_eval/2025-08-26_cross_study_imp_BB_ablation'
bb_markers = [*adata_sub.var_names[adata_sub.var.all(axis = 1)]]

train_kwargs = {
        'max_epochs': 400,
        'n_epochs_kl_warmup': 10}

# iteratively mask markers
for mask_marker in bb_markers:

    adata_masked = cytovi.pp.mask_markers(adata_sub, markers = mask_marker, batch_key = 'batch', masked_batch = 1)
            
    cytovi.CytoVI.setup_anndata(
        adata_masked,
        layer="scaled",
        batch_key='batch',
    )

    model = cytovi.CytoVI(adata_masked)
    model.train(**train_kwargs)
    model.save(f'{model_path}/2025-08-26_2025-08-26_cross_study_imp_BB_ablation_{mask_marker}', overwrite=True)

    adata_sub.layers[f'imputed_{mask_marker}'] = model.get_normalized_expression(adata_masked, n_samples = 10)

# add baseline model
cytovi.CytoVI.setup_anndata(
    adata_sub,
    layer="scaled",
    batch_key='batch',
)

model = cytovi.CytoVI(adata_sub)
model.train(**train_kwargs)
model.save(f'{model_path}/2025-08-26_2025-08-26_cross_study_imp_BB_ablation_baseline', overwrite=True)

adata_sub.layers['imputed_baseline'] = model.get_normalized_expression(adata_sub, n_samples = 10)

adata_sub.write(f'{model_path}/2025-08-26_cross_study_imp_BB_ablation_results_adata.h5ad')