library(cyCombine)
library(batchelor)
library(dplyr)

model_path  <-  '/home/projects/amit/floriani/Lab/PROJECTS/FlowVI/models/runtime_eval/2025-08-18_Nunez'

time_df <- data.frame(
)

obs_range <- c(1000, 10000, 100000, 1000000)

for (n_obs in obs_range){
    expr  <- read.csv(paste0(model_path, '/2025-08-18_Nunez_Runtime_', format(n_obs, scientific = FALSE), '_expr.csv'))[, -1]
    obs  <- read.csv(paste0(model_path, '/2025-08-18_Nunez_Runtime_',format(n_obs, scientific = FALSE), '_obs.csv'))[, -1]
    markers  <- colnames(expr)

    expr_batch1  <- t(expr[obs$batch == 0, ])
    expr_batch2  <- t(expr[obs$batch == 1, ])

    # FastMNN
    start <- Sys.time()
    out  <- fastMNN(expr_batch1, expr_batch2, d = 20)
    embd  <-  reducedDim(out, "corrected")
    end <- Sys.time()

    time_fastmnn  <- end - start
    method  <- paste0('FastMNN_', n_obs)

    time_df <- rbind(time_df, data.frame(
    method = method,
    time = time_fastmnn
    ))

    write.csv(embd, paste0(model_path, '/2025-08-18_Nunez_Runtime_', format(n_obs, scientific = FALSE), '_embd_fastmnn.csv'), row.names = FALSE)

    # cycombine
    expr[,'batch']  <- obs$batch

    start  <- Sys.time()
    labels  <- expr %>% 
        cyCombine::normalize(markers = markers, norm_method = 'scale') %>%
        create_som(markers = markers, rlen = 10)

    corrected <- expr %>%
        correct_data(label = labels, markers = markers)
    end  <- Sys.time()

    time_cycombine  <- end - start

    method  <- paste0('Cycombine_', n_obs)
    time_df <- rbind(time_df, data.frame(
    method = method,
    time = time_cycombine
    ))

    write.csv(corrected, paste0(model_path, '/2025-08-18_Nunez_Runtime_', format(n_obs, scientific = FALSE), '_embd_cycombine.csv'), row.names = FALSE)    
}

# Save the time data frame
write.csv(time_df, paste0(model_path, '/2025-08-18_Nunez_Runtime_time_FastMNN_Cycombine.csv'), row.names = FALSE)