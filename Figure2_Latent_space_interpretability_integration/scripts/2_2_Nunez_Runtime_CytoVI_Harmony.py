import cytovi
import anndata as ad
import time
import pandas as pd


model_path = '/home/projects/amit/floriani/Lab/PROJECTS/FlowVI/models/runtime_eval/2025-08-18_Nunez'
adata = ad.read_h5ad('/home/projects/amit/floriani/Lab/PROJECTS/FlowVI/data/2024-01-16_model_eval_multi_batch/2024-01-26_nunez_norm.h5ad')

# settings
obs_range = (1000, 10000, 100000, 1000000)
train_kwargs = {
        "n_epochs_kl_warmup": 50
        }

time_dict = {}

for n_obs in obs_range:
    # subsample
    adata_sub = cytovi.pp.subsample(adata, n_obs=n_obs, groupby='batch')

    # CytoVI
    t0 = time.perf_counter()

    cytovi.CytoVI.setup_anndata(
        adata_sub,
        layer="scaled",
        batch_key='batch'
    )

    model = cytovi.CytoVI(adata_sub)
    model.train(**train_kwargs)

    adata_sub.obsm['X_CytoVI'] = model.get_latent_representation()

    t1 = time.perf_counter()

    time_dict[f'CytoVI_{n_obs}'] = t1 - t0

    # harmony
    from harmony import harmonize
    t0 = time.perf_counter()
    adata_sub.obsm['X_harmony'] = harmonize(adata_sub.layers['scaled'], adata_sub.obs, batch_key = 'batch')
    t1 = time.perf_counter()

    time_dict[f'Harmony_{n_obs}'] = t1 - t0

    adata_sub.write(f'{model_path}/2025-08-18_Nunez_Runtime_{n_obs}.h5ad')

# save results
results = pd.DataFrame.from_dict(time_dict, orient='index', columns=['time'])
results.to_csv(f'{model_path}/2025-08-18_Nunez_Runtime_CytoVI_Harmony.csv')