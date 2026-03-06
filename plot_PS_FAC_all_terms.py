import os
import matplotlib.pyplot as plt
from src.compute_routines.compute_PiD_functions import *
from src.utils.hdf_to_df import hdf_to_df
from src.utils.methods import norm
from datetime import datetime

base_dir = rf'/home/sroy/Documents/IWF_research/MMS_events/MSH_Reconnection'

# interval = '20170617_202403_20170617_202411'
# interval = '20180415_043241_20180415_043246'
interval = '20151209_050300_20151209_050400'

fname = os.path.join(base_dir, rf'{interval}.h5')

df_dict = hdf_to_df(fname)

B = df_dict['b_gse_1']

P_FAC = compute_P_FAC(fname, species='ion', probe=1)
gradv_FAC = compute_gradv_FAC(fname, species='ion', probe=1)

T1 = P_FAC['bb'] * gradv_FAC['bb']
T2 = P_FAC['kk'] * gradv_FAC['kk']
T3 = P_FAC['nn'] * gradv_FAC['nn']

T4 = P_FAC['bk'] * gradv_FAC['kb']
T5 = P_FAC['bn'] * gradv_FAC['nb']

T6 = P_FAC['bk'] * gradv_FAC['bk']
T7 = P_FAC['bn'] * gradv_FAC['bn']

T8 = P_FAC['kn'] * gradv_FAC['kn']
T9 = P_FAC['kn'] * gradv_FAC['nk']

fig, axs = plt.subplots(3, 3, sharex=True)

axs[0, 0].plot(-T1, label='T1')
axs[0, 1].plot(-T2, label='T2')
axs[0, 2].plot(-T3, label='T3')

axs[1, 0].plot(-T4, label='T4')
axs[1, 1].plot(-T5, label='T5')
axs[1, 2].plot(-T6, label='T6')

axs[2, 0].plot(-T7, label='T7')
axs[2, 1].plot(-T8, label='T8')
axs[2, 2].plot(-T9, label='T9')

start_time = datetime(2015, 12, 9, 5, 3, 56)
end_time = datetime(2015, 12, 9, 5, 3, 58)

for i in range(3):
    for j in range(3):
        axs[i, j].legend()
        axs[i, j].set_xlim(start_time, end_time)

plt.suptitle('MMS1')

plt.show()