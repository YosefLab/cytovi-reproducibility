# run mario 
import pandas as pd
import time
from mario.match import Mario 
from mario.match import pipelined_mario

# read data for MARIO
model_path = '/home/projects/amit/floriani/Lab/PROJECTS/FlowVI/models/cross_tech_int/2025-08-20_Levine_Stuart'
expr_df = pd.read_csv(f'{model_path}/2026-01-13_Levine_Stuart_combined_adata_MARIO_expr.csv', index_col=0)
obs = pd.read_csv(f'{model_path}/2026-01-13_Levine_Stuart_combined_adata_MARIO_obs.csv', index_col=0)


expr_batch1 = expr_df.loc[(obs['batch'] == 0).values].copy()
expr_batch1 = expr_batch1.loc[:, expr_batch1.sum(axis=0) != 0]

expr_batch2 = expr_df.loc[(obs['batch'] == 1).values].copy()
expr_batch2 = expr_batch2.loc[:, expr_batch2.sum(axis=0) != 0]


# settings
obs_range = (500, 5000, 50000)

time_dict = {}

for n_obs in obs_range:
    print('Running MARIO for n_obs =', n_obs)
    n_per_batch = n_obs // 2

    # subsample each batch directly
    df1 = expr_batch1.sample(n=min(n_per_batch, len(expr_batch1)), random_state=0).copy()
    df2 = expr_batch2.sample(n=min(n_per_batch, len(expr_batch2)), random_state=0).copy()

    # MARIO
    t0 = time.perf_counter()
    final_matching_lst, embedding_lst = pipelined_mario(data_lst=[df1, df2])
    t1 = time.perf_counter()

    time_dict[f'MARIO_{n_obs}'] = t1 - t0

    # MARIO without batching
    t0 = time.perf_counter()
    final_matching_lst_no_batch, embedding_lst_no_batch = pipelined_mario(data_lst=[df1, df2], n_batches=1)
    t1 = time.perf_counter()

    time_dict[f'MARIO_no_batch_{n_obs}'] = t1 - t0

# save results
results = pd.DataFrame.from_dict(time_dict, orient='index', columns=['time'])
results.to_csv(f'{model_path}/2026-03-31_Levine_Stuart_Runtime_MARIO.csv')
