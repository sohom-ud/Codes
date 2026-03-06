import os
import matplotlib.pyplot as plt
from src.compute_routines.compute_PiD_functions import *
from src.utils.hdf_to_df import hdf_to_df
from src.utils.methods import norm
import datetime

# import pickle
# import inspect

# def load_vars(filename):
#     caller_vars = inspect.stack()[1].frame.f_locals
#     with open(filename, 'rb') as f:
#         saved_vars = pickle.load(f)
#     caller_vars.update(saved_vars)

# load_vars('PDU.pkl')

base_dir = rf'/home/sroy/Documents/IWF_research/MMS_events/Tail_Reconnection'

# base_dir = rf'/home/sroy/Documents/IWF_research/MMS_events/MSH_Reconnection'

interval = '20170617_202403_20170617_202411'
# interval = r'20170711_223330_20170711_223430'

# interval = r'20151209_050304_20151209_'


fname = os.path.join(base_dir, rf'{interval}.h5')

df_dict = hdf_to_df(fname)

B = df_dict['b_gse_1']

ve1 = compute_v_FAC(fname, species='elc', probe=1)

B_reselectron = resample(B, ve1)
B_reselectron = B_reselectron.drop('mag', axis=1)

# PS1_ion = compute_PS1(fname, species='ion', probe=4, reselectron=True)
# PS2_ion = compute_PS2(fname, species='ion', probe=4, reselectron=True)
# PS3_ion = compute_PS3(fname, species='ion', probe=4, reselectron=True)
# PS4_ion = compute_PS4(fname, species='ion', probe=4, reselectron=True)

PS1_elc = compute_PS1(fname, species='elc', probe=4, reselectron=True)
PS2_elc = compute_PS2(fname, species='elc', probe=4, reselectron=True)
PS3_elc = compute_PS3(fname, species='elc', probe=4, reselectron=True)
PS4_elc = compute_PS4(fname, species='elc', probe=4, reselectron=True)

# PS5_ion, PS6_ion, PS7_ion, PS8_ion = compute_geom_PS(fname, species='elc', probe=1, reselectron=True)

# PS5_elc, PS6_elc, PS7_elc, PS8_elc = compute_geom_PS(fname, species='elc', probe=1, reselectron=True)

fig, axs = plt.subplots(6, 1, sharex=True)

axs[0].plot(B.index, B['x'], 'r', label='x')
axs[0].plot(B.index, B['y'], 'g', label='y')
axs[0].plot(B.index, B['z'], 'b', label='z')

# axs[1].plot(PS1_ion.index, -PS1_ion.values, 'k', label='-PS1')
# axs[2].plot(PS2_ion.index, -PS2_ion.values, 'k', label='-PS2')
# axs[3].plot(PS3_ion.index, -PS3_ion.values, 'k', label='-PS3')
# axs[4].plot(PS4_ion.index, -PS4_ion.values, 'k', label='-PS4')

axs[1].plot(PS1_elc.index, -PS1_elc.values, 'k', label='-PS1')
axs[2].plot(PS2_elc.index, -PS2_elc.values, 'k', label='-PS2')
axs[3].plot(PS3_elc.index, -PS3_elc.values, 'k', label='-PS3')
axs[4].plot(PS4_elc.index, -PS4_elc.values, 'k', label='-PS4')

# axs[1].plot(PS5_ion.index, -PS5_ion.values, 'k', label='-PS5')
# axs[2].plot(PS6_ion.index, -PS6_ion.values, 'k', label='-PS6')
# axs[3].plot(PS7_ion.index, -PS7_ion.values, 'k', label='-PS7')
# axs[4].plot(PS8_ion.index, -PS8_ion.values, 'k', label='-PS8')

axs[5].plot(ve1.index, ve1['b'], 'r', label='b')
axs[5].plot(ve1.index, ve1['k'], 'g', label='k')
axs[5].plot(ve1.index, ve1['n'], 'b', label='n')

for ax in axs[1:]:
    ax.axhline(0, color='r')
    ax.legend()

axs[0].legend()

plt.suptitle('MMS1')

plt.show()