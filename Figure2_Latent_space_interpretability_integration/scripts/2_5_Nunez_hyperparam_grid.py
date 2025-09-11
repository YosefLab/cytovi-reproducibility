import cytovi
import anndata as ad

import itertools

model_path = '/home/projects/amit/floriani/Lab/PROJECTS/FlowVI/models/hyperopt/2025-08-15_hypergrid_nlatent_nhidden'

adata = ad.read_h5ad('/home/projects/amit/floriani/Lab/PROJECTS/FlowVI/data/2024-01-16_model_eval_multi_batch/2024-01-26_nunez_norm.h5ad')

n_hidden = (64, 128, 256)
n_latent = (10, 20, 40)

train_kwargs = {
        "n_epochs_kl_warmup": 250
        }

combinations = list(itertools.product(n_hidden, n_latent))

for nh, nl in combinations:
    print(f"n_hidden={nh}, n_latent={nl}")

    cytovi.CytoVI.setup_anndata(
        adata,
        layer="scaled",
        batch_key='batch'
    )

    model_key = f'nlatent_{nl}_nhidden_{nh}'

    model = cytovi.CytoVI(adata, n_hidden=nh, n_latent=nl)
    model.train(**train_kwargs)
    model.save(f'{model_path}{model_key}', overwrite=True)

    adata.obsm[f'X_{model_key}'] = model.get_latent_representation()

adata.write(f'{model_path}/2025-08-15_Nunez_hyperparam_grid_latents.h5ad')