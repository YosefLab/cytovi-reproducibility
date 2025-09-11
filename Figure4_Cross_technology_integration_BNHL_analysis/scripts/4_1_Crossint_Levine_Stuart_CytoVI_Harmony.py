import cytovi
import anndata as ad
import pandas as pd


model_path = '/home/projects/amit/floriani/Lab/PROJECTS/FlowVI/models/cross_tech_int/2025-08-20_Levine_Stuart'
adata = ad.read_h5ad(f'{model_path}/2025-08-20_Levine_Stuart_combined_adata.h5ad')

# settings
train_kwargs = {
        "n_epochs_kl_warmup": 50,
        'plan_kwargs': {"scale_adversarial_loss": 1}
        }


cytovi.CytoVI.setup_anndata(
        adata,
        layer="std_scaled",
        batch_key='batch'
    )

model = cytovi.CytoVI(adata)
model.train(**train_kwargs)

adata.obsm['X_CytoVI'] = model.get_latent_representation()


# harmony
from harmony import harmonize
adata_bb = adata[:, model.backbone_markers].copy()
adata.obsm['X_harmony'] = harmonize(adata_bb.layers['std_scaled'], adata_bb.obs, batch_key = 'batch')

adata.write(f'{model_path}/2025-08-20_Levine_Stuart_combined_adata_cytovi_harmony.h5ad')
