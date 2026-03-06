import os
import matplotlib.pyplot as plt
from src.compute_routines.compute_PiD_functions import *
from src.utils.hdf_to_df import hdf_to_df
from src.compute_routines.scale_filter_functions import compute_Pi_uu, compute_Pi_bb
import datetime

base_dir = rf'/home/sroy/Documents/IWF_research/MMS_events/Dipolarization_Fronts'

data_dir = r'Data_Zhong2019'
plot_dir = r'Figures_Zhong2019'

for interval in os.listdir(os.path.join(base_dir, data_dir)):

        fname = os.path.join(base_dir, data_dir, rf'{interval}')

        df_dict = hdf_to_df(fname)

        B = df_dict['b_gse_1']
        ni = df_dict['N_ion_1']
        ne = df_dict['N_elc_1']
        E = df_dict['edp_dce_gse_1']

        P_FAC = compute_avg_P_FAC(fname, species='ion', reselectron=True)
        gradv_FAC = compute_gradv_FAC(fname, species='ion', reselectron=True)

        PS1 = P_FAC['bb'] * gradv_FAC['bb']
        PS2 = P_FAC['kk'] * gradv_FAC['kk'] + P_FAC['nn'] * gradv_FAC['nn']
        PS3 = P_FAC['kb'] * gradv_FAC['kb'] + P_FAC['nb'] * gradv_FAC['nb']
        PS4 = P_FAC['bk'] * gradv_FAC['bk'] + P_FAC['bn'] * gradv_FAC['bn'] + \
                P_FAC['kn'] * gradv_FAC['kn'] + P_FAC['nk'] * gradv_FAC['nk']

        fig, axs = plt.subplots(5, 1, sharex=True)

        axs[0].plot(B.index, B['x'], 'r')
        axs[0].plot(B.index, B['y'], 'g')
        axs[0].plot(B.index, B['z'], 'b')

        axs[1].plot(-PS1, 'k', label='PS1')
        axs[2].plot(-PS2, 'k', label='PS2')
        axs[3].plot(-PS3, 'k', label='PS3')
        axs[4].plot(-PS4, 'k', label='PS4')

        for ax in axs:
            ax.legend()

        plt.savefig(os.path.join(base_dir, 'Figures_Zhong2019', rf'{interval}_PS1to4_FAC.png'))