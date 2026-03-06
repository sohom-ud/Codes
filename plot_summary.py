import matplotlib.pyplot as plt
from src.utils.hdf_to_df import hdf_to_df
import pandas

interval = r'20170711_223330_20170711_223430'
fname = rf'/home/sroy/Documents/IWF_research/MMS_events/Tail_Reconnection/{interval}.h5'

df_dict = hdf_to_df(fname)

B = df_dict['b_gse_1']

ni = df_dict['N_ion_1']
ne = df_dict['N_elc_1']

vi = df_dict['v_spincorr_ion_1']
ve = df_dict['v_spincorr_elc_1']

fig, axs = plt.subplots(4, 1, sharex=True)

axs[0].plot(B.index, B['x'], 'r')
axs[0].plot(B.index, B['y'], 'g')
axs[0].plot(B.index, B['z'], 'b')

axs[1].plot(ni.index, ni.values, 'k')
axs[1].plot(ne.index, ne.values, 'b')

axs[2].plot(vi.index, vi['x'], 'r')
axs[2].plot(vi.index, vi['y'], 'g')
axs[2].plot(vi.index, vi['z'], 'b')

axs[3].plot(ve.index, ve['x'], 'r')
axs[3].plot(ve.index, ve['y'], 'g')
axs[3].plot(ve.index, ve['z'], 'b')

plt.show()