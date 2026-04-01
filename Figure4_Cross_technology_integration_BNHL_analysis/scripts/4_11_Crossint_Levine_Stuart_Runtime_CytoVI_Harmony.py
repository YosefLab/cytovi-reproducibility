# runtime benchmarking cytovi, harmony
import cytovi
import anndata as ad
import time
import pandas as pd


# read anndata
model_path = '/home/projects/amit/floriani/Lab/PROJECTS/FlowVI/models/cross_tech_int/2025-08-20_Levine_Stuart'
adata = ad.read_h5ad(f'{model_path}/2026-01-14_Levine_Stuart_combined_adata_MARIO_filtered_nobatch.h5ad')

# settings
obs_range = (500, 5000, 50000)
train_kwargs = {
        "n_epochs_kl_warmup": 50,
        'plan_kwargs': {"scale_adversarial_loss": 1}
        }

time_dict = {}

for n_obs in obs_range:
    # subsample
    adata_sub = cytovi.pp.subsample(adata, n_obs=n_obs, groupby='batch')

    # CytoVI
    t0 = time.perf_counter()

    cytovi.CytoVI.setup_anndata(
        adata_sub,
        layer="std_scaled",
        batch_key='batch'
    )

    model = cytovi.CytoVI(adata_sub)
    model.train(**train_kwargs)

    adata_sub.obsm['X_CytoVI'] = model.get_latent_representation()

    t1 = time.perf_counter()

    time_dict[f'CytoVI_{n_obs}'] = t1 - t0

    # harmony
    from harmony import harmonize
    adata_sub_bb = adata_sub[:, model.backbone_markers].copy()
    t0 = time.perf_counter()
    adata_sub_bb.obsm['X_harmony'] = harmonize(adata_sub_bb.layers['std_scaled'], adata_sub_bb.obs, batch_key = 'batch')
    t1 = time.perf_counter()

    time_dict[f'Harmony_{n_obs}'] = t1 - t0

# save results
results = pd.DataFrame.from_dict(time_dict, orient='index', columns=['time'])
results.to_csv(f'{model_path}/2026-03-31_Levine_Stuart_Runtime_CytoVI_Harmony.csv')