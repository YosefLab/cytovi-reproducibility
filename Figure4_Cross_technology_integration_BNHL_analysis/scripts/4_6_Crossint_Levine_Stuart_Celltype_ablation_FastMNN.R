library(batchelor)
library(dplyr)

ct_model_path = '/home/projects/amit/floriani/Lab/PROJECTS/FlowVI/models/cross_tech_int/2025-08-20_Levine_Stuart/2025-08-22_Celltype_ablation'

obs  <- read.csv(paste0(ct_model_path, '/../2025-08-20_Levine_Stuart_combined_adata_obs.csv'))[, -1]
cell_types  <- unique(obs$harmonized_labels)


for (ct in cell_types){
    if (ct == 'Unknown'){
        next
    }
    obs_ct  <- read.csv(paste0(ct_model_path, '/2025-08-22_Levine_Stuart_combined_adata_', ct, '_ablation_obs.csv'))[, -1]
    expr_ct  <- read.csv(paste0(ct_model_path, '/2025-08-22_Levine_Stuart_combined_adata_', ct, '_ablation_expr.csv'))[, -1]

    markers  <- colnames(expr_ct)

    expr_batch1  <- t(expr_ct[obs_ct$batch == 0, ])
    expr_batch2  <- t(expr_ct[obs_ct$batch == 1, ])

    # FastMNN
    out  <- fastMNN(expr_batch1, expr_batch2, d = 10)
    embd  <-  reducedDim(out, "corrected")

    write.csv(embd, paste0(ct_model_path, '/2025-08-25_Levine_Stuart_combined_adata_', ct, '_ablation_embd.csv'), row.names = FALSE)
}