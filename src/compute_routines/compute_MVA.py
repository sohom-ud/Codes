from src.utils.hdf_to_df import hdf_to_df
from src.utils.methods import mult
import pandas as pd
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt

comps = ['x', 'y', 'z']

fname = r'/home/sroy/Documents/MMS_events/MP_Reconnection/20160214_204150_20160214_204210.h5'

df_dict = hdf_to_df(fname)

B1 = df_dict['b_gse_1']
B1 = B1.drop('mag', axis=1)

start_MVA = datetime(2016, 2, 14, 20, 41, 55)
end_MVA = datetime(2016, 2, 14, 20, 41, 57, 500000)

M = np.cov(B1[start_MVA:end_MVA].values.T, bias=True)

eigvals, eigvecs = np.linalg.eigh(M)

eigvals = eigvals[::-1]
eigvecs = eigvecs[:, ::-1]

B_LMN = pd.DataFrame(B1.values @ eigvecs, index=B1.index, columns=['L', 'M', 'N'])

fig, axs = plt.subplots(2, 1, sharex=True)
axs[0].plot(B1['x'], 'r')
axs[0].plot(B1['y'], 'g')
axs[0].plot(B1['z'], 'b')

axs[1].plot(B_LMN['L'], 'r')
axs[1].plot(B_LMN['M'], 'g')
axs[1].plot(B_LMN['N'], 'b')

xlim_start = datetime(2016, 2, 14, 20, 41, 54)
xlim_end = datetime(2016, 2, 14, 20, 41, 58)

plt.xlim(xlim_start, xlim_end)

plt.show()