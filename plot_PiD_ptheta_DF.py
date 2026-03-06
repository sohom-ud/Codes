import os
import matplotlib.pyplot as plt
from src.compute_routines.compute_PiD_functions import *
from src.utils.hdf_to_df import hdf_to_df
from src.compute_routines.scale_filter_functions import compute_Pi_uu, compute_Pi_bb
import datetime

base_dir = rf'/home/sroy/Documents/IWF_research/MMS_events/Dipolarization_Fronts/Data_Zhong2019'

interval = '20170519_094558_20170519_094758'

fname = os.path.join(base_dir, rf'{interval}.h5')

df_dict = hdf_to_df(fname)

B = df_dict['b_gse_1']
ni = df_dict['N_ion_1']
ne = df_dict['N_elc_1']
E = df_dict['edp_dce_gse_1']

# PiDi = compute_PiD(fname, species='ion')
# PiDe = compute_PiD(fname, species='elc')

# pthi = compute_ptheta(fname, species='ion')
# pthe = compute_ptheta(fname, species='elc')

# Q_D_ion = compute_Q_D(fname, species='ion')
# Q_omega_ion = compute_Q_omega(fname, species='ion')

# Q_D_elc = compute_Q_D(fname, species='elc')
# Q_omega_elc = compute_Q_omega(fname, species='elc')

# Q_j = compute_Q_j(fname)

# Pi_uu_i = compute_Pi_uu(fname, t=1000, unit='ms', win_gauss=0, species='ion', reselectron=True)
# Pi_uu_e = compute_Pi_uu(fname, t=1000, unit='ms', win_gauss=0, species='elc', reselectron=True)

Pi_bb_i = compute_Pi_bb(fname, t=1000, unit='ms', win_gauss=0, probe=1, species='ion', reselectron=True)
Pi_bb_e = compute_Pi_bb(fname, t=1000, unit='ms', win_gauss=0, probe=1, species='elc', reselectron=True)

# PDU_ion = compute_PDU(fname, species='ion')
# PDU_elc = compute_PDU(fname, species='elc')

# avgP = compute_avg_P(fname, species='elc')

fig, axs = plt.subplots(10, 1, sharex=True)

axs[0].plot(B.index, B['x'], 'r')
axs[0].plot(B.index, B['y'], 'g')
axs[0].plot(B.index, B['z'], 'b')

axs[1].plot(ni.index, ni.values, 'k')
axs[1].plot(ne.index, ne.values, 'r')

axs[2].plot(E.index, E['x'], 'r')
axs[2].plot(E.index, E['y'], 'g')
axs[2].plot(E.index, E['z'], 'b')

# axs[2].plot(Q_D_ion.index, Q_D_ion.values, 'k')
# axs[3].plot(Q_omega_ion.index, Q_omega_ion.values, 'k')

# axs[4].plot(Q_D_elc.index, Q_D_elc.values, 'k')
# axs[5].plot(Q_omega_elc.index, Q_omega_elc.values, 'k')

axs[4].plot(Pi_bb_i.index, Pi_bb_i.values, 'k')
axs[5].plot(Pi_bb_e.index, Pi_bb_e.values, 'r')

# axs[5].plot(Pi_uu_i.index, Pi_uu_i.values, 'k')
# axs[5].plot(Pi_uu_e.index, Pi_uu_e.values, 'r')

# axs[6].plot(PiDi.index, -PiDi.cumsum(), 'k')
# axs[7].plot(pthi.index, -pthi.cumsum(), 'k')
# axs[8].plot(PiDe.index, -PiDe.cumsum(), 'k')
# axs[9].plot(pthe.index, -pthe.cumsum(), 'k')

# axs[4].plot(PDU_ion.index, -PDU_ion.values, 'r')
# axs[5].plot(PDU_elc.index, -PDU_elc.values, 'r')

# start_time = datetime.datetime(2017, 5, 19, 9, 58, 22)
# end_time = datetime.datetime(2017, 6, 24, 23, 58, 40)

# plt.xlim(start_time, end_time)

plt.show()