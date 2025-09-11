from harmony import harmonize
import cytovi
import anndata as ad

model_path = '/home/projects/amit/floriani/Lab/PROJECTS/FlowVI/models/cross_tech_int/2025-08-20_Levine_Stuart/2025-08-21_Noise_benchmarking'

adata = ad.read_h5ad(f'{model_path}/../2025-08-21_Levine_Stuart_combined_adata_all_methods_noise_added.h5ad')

sigma_range = (0.1, 0.2, 0.5, 1, 2)

train_kwargs = {
        "n_epochs_kl_warmup": 50,
        'plan_kwargs': {"scale_adversarial_loss": 1}
        }

bb_markers = adata.var_names[adata.var.all(axis=1)]
adata_bb = adata[:, bb_markers].copy()

for sigma in sigma_range:
        print(f"Training CytoVI/Harmony with noise sigma: {sigma}")

        cytovi.CytoVI.setup_anndata(
                adata,
                layer=f"std_scaled_{sigma}",
                batch_key='batch'
        )

        model = cytovi.CytoVI(adata)
        model.train(**train_kwargs)
        model.save(f'{model_path}/CytoVI_{sigma}', overwrite=True)
        adata.obsm[f'X_CytoVI_{sigma}'] = model.get_latent_representation()

        # harmony
        adata.obsm[f'X_harmony_{sigma}'] = harmonize(adata_bb.layers[f'std_scaled_{sigma}'], adata_bb.obs, batch_key = 'batch')


adata.write(f'{model_path}/2025-08-21_Levine_Stuart_combined_adata_all_methods_noise_added_CytoVI_Harmony.h5ad')