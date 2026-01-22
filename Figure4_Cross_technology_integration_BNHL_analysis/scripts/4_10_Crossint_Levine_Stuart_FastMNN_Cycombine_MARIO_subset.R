library(cyCombine)
library(batchelor)
library(dplyr)

model_path  <-  '/home/projects/amit/floriani/Lab/PROJECTS/FlowVI/models/cross_tech_int/2025-08-20_Levine_Stuart'
adata_name  <- '2026-01-14_Levine_Stuart_combined_adata_MARIO_filtered_nobatch_bb'

expr  <- read.csv(paste0(model_path, '/', adata_name, '_expr.csv'))[, -1]
obs  <- read.csv(paste0(model_path, '/', adata_name, '_obs.csv'))[, -1]
markers  <- colnames(expr)

expr_batch1  <- t(expr[obs$batch == 0, ])
expr_batch2  <- t(expr[obs$batch == 1, ])

# FastMNN
out  <- fastMNN(expr_batch1, expr_batch2, d = 10)
embd  <-  reducedDim(out, "corrected")
write.csv(embd, paste0(model_path, '/', adata_name, '_FastMNN.csv'), row.names = FALSE)

# cycombine
expr[,'batch']  <- obs$batch

labels  <- expr %>% 
    create_som(markers = markers, rlen = 10)

corrected <- expr %>%
    correct_data(label = labels, markers = markers)

write.csv(corrected, paste0(model_path, '/', adata_name, '_Cycombine.csv'), row.names = FALSE)    
