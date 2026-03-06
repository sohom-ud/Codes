import os
import matplotlib.pyplot as plt
from src.compute_routines.compute_PiD_functions import *
from src.utils.hdf_to_df import hdf_to_df
import datetime

import pickle
import inspect

def load_vars(filename):
    caller_vars = inspect.stack()[1].frame.f_locals
    with open(filename, 'rb') as f:
        saved_vars = pickle.load(f)
    caller_vars.update(saved_vars)

load_vars('PDU.pkl')

base_dir = rf'/home/sroy/Documents/IWF_research/MMS_events/Tail_Reconnection'

interval = '20170617_202403_20170617_202411'

fname = os.path.join(base_dir, rf'{interval}.h5')

df_dict = hdf_to_df(fname)

B = df_dict['b_gse_1']
ni = df_dict['N_ion_1']
ne = df_dict['N_elc_1']

# r1 = df_dict['r_gse_1_res_vi']
# r1 = r1.drop('mag', axis=1)

# plt.plot(r_gse[0].data[:, 0]/(r1['x'].values * 1000))

# fig, axs = plt.subplots(3, 1)

# axs[0].plot(r1['x'])
# axs[1].plot(r1['y'])
# axs[2].plot(r1['z'])

# PiD_shear_new = out_pdu['pid_shear'].data * 1e9
# PDU_new = out_pdu['PDU'].data * 1e9
# PiD_normal_new = out_pdu['pid_norm'].data * 1e9

# PDU_i = compute_PDU(fname, species='ion', reselectron=False)
# PDU_e = compute_PDU(fname, species='elc', reselectron=False)

# PiDi_shear = compute_PiD_shear(fname, species='ion', reselectron=False)
# PiDe_shear = compute_PiD_shear(fname, species='elc', reselectron=False)

# PiDe_normal = compute_PiD_normal(fname, species='elc', reselectron=False)

P_FAC = compute_avg_P_FAC(fname, species='elc', reselectron=True)
gradv_FAC = compute_gradv_FAC(fname, species='elc', reselectron=True)

PS1 = P_FAC['bb'] * gradv_FAC['bb']
PS2 = P_FAC['kk'] * gradv_FAC['kk'] + P_FAC['nn'] * gradv_FAC['nn']
PS3 = P_FAC['kb'] * gradv_FAC['kb'] + P_FAC['nb'] * gradv_FAC['nb']
PS4 = P_FAC['bk'] * gradv_FAC['bk'] + P_FAC['bn'] * gradv_FAC['bn'] + \
          P_FAC['kn'] * gradv_FAC['kn'] + P_FAC['nk'] * gradv_FAC['nk']


# PiDi = compute_PiD(fname, species='ion', reselectron=False)
# PiDe = compute_PiD(fname, species='elc', reselectron=False)

# pthi = compute_ptheta(fname, species='ion', reselectron=False)
# pthe = compute_ptheta(fname, species='elc', reselectron=False)

fig, axs = plt.subplots(4, 1, sharex=True)

# axs[0].plot(PiDe_shear.index, -PiDe_shear.values, 'b')
# axs[1].plot(PDU_e.index, -PDU_e.values, 'r')
# axs[2].plot(PiDe_normal.index, -PiDe_normal.values, 'g')

# axs[0].plot((-PiDe_shear.values.flatten() - PiD_shear_new)/(np.abs(PiD_shear_new) + 1e-5))
# axs[1].plot((-PDU_e.values.flatten() - PDU_new)/(np.abs(PDU_new) + 1e-5))
# axs[2].plot((-PiDe_normal.values.flatten() - PiD_normal_new)/(np.abs(PiD_normal_new) + 1e-5))

# axs[0].plot(-PiDe_shear.cumsum().values, label=r'PiDe_shear')
# axs[0].plot(np.cumsum(PiD_shear_new))

# axs[1].plot(-PDU_e.cumsum().values, label=r'PDU_e')
# axs[1].plot(np.cumsum(PDU_new))

# axs[2].plot(-PiDe_normal.cumsum().values, label=r'PiDe_normal')
# axs[2].plot(np.cumsum(PiD_normal_new))

axs[0].plot(PS1.index, -PS1.values, 'k')
axs[1].plot(PS2.index, -PS2.values, 'k')
axs[2].plot(PS3.index, -PS3.values, 'k')
axs[3].plot(PS4.index, -PS4.values, 'k')

# axs[2].plot(PiDi_shear.index, -PiDi_shear.values, 'b')

# axs[0].plot(pthi.index, -pthi.values, 'r')
# axs[1].plot(PiDi.index, -PiDi.values, 'b')

# axs[2].plot(pthe.index, -pthe.values, 'r')
# axs[3].plot(PiDe.index, -PiDe.values, 'b')

# axs[0].plot(D_ion.index, D_ion['xx'], 'k', label='xx')
# axs[1].plot(D_ion.index, D_ion['yy'], 'k', label='yy')
# axs[2].plot(D_ion.index, D_ion['zz'], 'k', label='zz')

# # axs[0].plot(D_ion.index, D_ion['xy'], 'k', label='xy')
# # axs[1].plot(D_ion.index, D_ion['xz'], 'k', label='xz')
# # axs[2].plot(D_ion.index, D_ion['yz'], 'k', label='yz')

for ax in axs:
    ax.axhline(0)
    ax.legend(fontsize=12)

# # axs[2].set_ylim(-1, 1)

# start_time = datetime.datetime(2017, 6, 17, 20, 24, 3)
# end_time = datetime.datetime(2017, 6, 17, 20, 24, 11)

# plt.xlim(start_time, end_time)

# axs[1].set_ylim(-0.4, 0.4)

# fig, axs = plt.subplots(4, 1, sharex=True)

# axs[0].plot(pthi.index, -pthi.values, 'r')

# axs[1].plot(pthe.index, -pthe.values, 'r')

# axs[2].plot(PiDi.index, -PiDi.values, 'b')
# axs[3].plot(PiDe.index, -PiDe.values, 'b')

plt.show()