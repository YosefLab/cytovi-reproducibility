library(batchelor)
library(dplyr)

noise_model_path = '/home/projects/amit/floriani/Lab/PROJECTS/FlowVI/models/cross_tech_int/2025-08-20_Levine_Stuart/2025-08-21_Noise_benchmarking'

obs  <- read.csv(paste0(noise_model_path, '/../2025-08-20_Levine_Stuart_combined_adata_obs.csv'))[, -1]

sigma_range  <-  c(0.1, 0.2, 0.5, 1, 2)

for (sigma in sigma_range){
    expr  <- read.csv(paste0(noise_model_path, '/2025-08-21_Levine_Stuart_combined_adata_all_methods_noise_added_sigma_', sigma, '_expr.csv'))[, -1]

    markers  <- colnames(expr)

    expr_batch1  <- t(expr[obs$batch == 0, ])
    expr_batch2  <- t(expr[obs$batch == 1, ])

    # FastMNN
    out  <- fastMNN(expr_batch1, expr_batch2, d = 10)
    embd  <-  reducedDim(out, "corrected")

    write.csv(embd, paste0(noise_model_path, '/2025-08-21_Levine_Stuart_combined_adata_all_methods_noise_added_sigma_', sigma, '_embd_fastmnn.csv'), row.names = FALSE)
}