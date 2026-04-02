from src.utils.hdf_to_df import hdf_to_df
from src.utils.methods import mult
import pandas as pd
import numpy as np
from datetime import datetime

comps = ['x', 'y', 'z']

fname = r'/home/sroy/Documents/MMS_events/Tail_Reconnection/20170810_121800_20170810_121900.h5'

df_dict = hdf_to_df(fname)

B1 = df_dict['b_gse_1']
B1 = B1.drop('mag', axis=1)

start_MVA = datetime(2017, 8, 10, 12, 18)
end_MVA = datetime(2017, 8, 10, 12, 19)

term1 = mult(B1, B1)[start_MVA:end_MVA].mean(axis=0)

term2 = pd.Series()

for i in comps:
    for j in comps:
        term2[f'{i}{j}'] = np.mean(B1[i][start_MVA:end_MVA]) * np.mean(B1[j][start_MVA:end_MVA])

M = (term1 - term2).values.reshape(3, 3)

eigvals, eigvecs = np.linalg.eig(M)