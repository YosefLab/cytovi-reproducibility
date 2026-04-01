# benchmarking of Cycombine & FastMNN
library(cyCombine)
library(batchelor)
library(dplyr)
library(SingleCellExperiment)

model_path <- '/home/projects/amit/floriani/Lab/PROJECTS/FlowVI/models/cross_tech_int/2025-08-20_Levine_Stuart'
adata_name <- '2026-01-14_Levine_Stuart_combined_adata_MARIO_filtered_nobatch_bb'

expr <- read.csv(paste0(model_path, '/', adata_name, '_expr.csv'))[, -1]
obs  <- read.csv(paste0(model_path, '/', adata_name, '_obs.csv'))[, -1]

markers <- colnames(expr)

expr_batch1 <- expr[obs$batch == 0, ]
expr_batch2 <- expr[obs$batch == 1, ]

obs_range <- c(500, 5000, 50000)

time_dict <- c()

for (n_obs in obs_range) {
  
  n_per_batch <- floor(n_obs / 2)
  
  # subsample each batch directly
  idx1 <- sample(seq_len(nrow(expr_batch1)), size = min(n_per_batch, nrow(expr_batch1)))
  idx2 <- sample(seq_len(nrow(expr_batch2)), size = min(n_per_batch, nrow(expr_batch2)))
  
  df1 <- expr_batch1[idx1, , drop = FALSE]
  df2 <- expr_batch2[idx2, , drop = FALSE]
  
  # FastMNN
  t0 <- Sys.time()
  
  out <- fastMNN(t(df1), t(df2), d = 10)
  embd <- reducedDim(out, "corrected")
  
  t1 <- Sys.time()
  
  time_dict[paste0("FastMNN_", n_obs)] <- as.numeric(difftime(t1, t0, units = "secs"))
  
  # cyCombine
  expr_sub <- rbind(df1, df2)
  expr_sub$batch <- c(rep(0, nrow(df1)), rep(1, nrow(df2)))
  
  t0 <- Sys.time()
  
  labels <- expr_sub %>%
    create_som(markers = markers, rlen = 10)
  
  corrected <- expr_sub %>%
    correct_data(label = labels, markers = markers)
  
  t1 <- Sys.time()
  
  time_dict[paste0("Cycombine_", n_obs)] <- as.numeric(difftime(t1, t0, units = "secs"))
}

# save results
results <- data.frame(time = as.numeric(time_dict))
rownames(results) <- names(time_dict)

write.csv(
  results,
  paste0(model_path, '/2026-03-31_Levine_Stuart_Runtime_FastMNN_Cycombine.csv')
)