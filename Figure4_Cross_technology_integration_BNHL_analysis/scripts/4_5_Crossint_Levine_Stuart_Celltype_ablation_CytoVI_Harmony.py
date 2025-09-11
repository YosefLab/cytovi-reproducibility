from harmony import harmonize
import cytovi
import anndata as ad

ct_model_path = '/home/projects/amit/floriani/Lab/PROJECTS/FlowVI/models/cross_tech_int/2025-08-20_Levine_Stuart/2025-08-22_Celltype_ablation'

adata = ad.read_h5ad(f'{ct_model_path}/../2025-08-21_Levine_Stuart_combined_adata_all_methods.h5ad')

train_kwargs = {
        "n_epochs_kl_warmup": 50,
        'plan_kwargs': {"scale_adversarial_loss": 1}
        }

bb_markers = adata.var_names[adata.var.all(axis=1)]

for ct in adata.obs['harmonized_labels'].cat.categories:
    if ct != 'Unknown':
        print(f"Training CytoVI/Harmony with ct: {ct}")

        adata_ct = ad.read_h5ad(f'{ct_model_path}/2025-08-22_Levine_Stuart_combined_adata_{ct}_ablation.h5ad')

        cytovi.CytoVI.setup_anndata(
                adata_ct,
                layer="std_scaled",
                batch_key='batch'
        )

        model = cytovi.CytoVI(adata_ct)
        model.train(**train_kwargs)
        model.save(f'{ct_model_path}/CytoVI_{ct}', overwrite=True)
        adata_ct.obsm['X_CytoVI'] = model.get_latent_representation()

        # harmony
        adata_ct_bb = adata_ct[:, bb_markers].copy()
        adata_ct.obsm['X_harmony'] = harmonize(adata_ct_bb.layers['std_scaled'], adata_ct_bb.obs, batch_key = 'batch')

        # save
        adata_ct.write(f'{ct_model_path}/2025-08-22_Levine_Stuart_combined_adata_{ct}_ablation_CytoVI_Harmony.h5ad')
