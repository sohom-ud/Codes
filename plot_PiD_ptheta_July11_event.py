import os
import matplotlib.pyplot as plt
from src.compute_routines.compute_PiD_functions import *
from src.utils.hdf_to_df import hdf_to_df
import datetime
from pyspedas.mms import fpi

trange = ['2017-07-11/22:33:30', '2017-07-11/22:34:30']

probe = 3

data = fpi(
    trange=trange,
    probe=probe,
    data_rate='brst',
    datatype=f"des-moms",
    center_measurement=True,
    notplot=True
)

Te_perp = data['mms3_des_tempperp_brst']

Te_para = data['mms3_des_temppara_brst']

base_dir = rf'/home/sroy/Documents/IWF_research/MMS_events/Tail_Reconnection'

interval = '20170711_223330_20170711_223430' # July 11 2017 event
# interval = '20170617_202403_20170617_202411' # Roberts 2023 event

fname = os.path.join(base_dir, rf'{interval}.h5')

df_dict = hdf_to_df(fname)

B = df_dict['b_gse_1']

ni = df_dict['N_ion_1']
ne = df_dict['N_elc_1']

vi = df_dict['v_spincorr_ion_1']
ve = df_dict['v_spincorr_elc_1']

PiDi = compute_PiD(fname, species='ion')
PiDe = compute_PiD(fname, species='elc')

pthi = compute_ptheta(fname, species='ion')
pthe = compute_ptheta(fname, species='elc')

avgP = compute_avg_P(fname, species='elc')

fig, axs = plt.subplots(5, 1, figsize=(6, 6), sharex=True)

plt.rcParams['font.family'] = 'serif'
# plt.rcParams['text.usetex'] = True

axs[0].plot(B.index, B['x'], 'r', label=r'x', lw=1)
axs[0].plot(B.index, B['y'], 'g', label=r'y', lw=1)
axs[0].plot(B.index, B['z'], 'b', label=r'z', lw=1)

axs[1].plot(vi.index, vi['x'], 'r', label=r'x', lw=1)
axs[1].plot(vi.index, vi['y'], 'g', label=r'y', lw=1)
axs[1].plot(vi.index, vi['z'], 'b', label=r'z', lw=1)

axs[2].plot(ve.index, ve['x'], 'r', label=r'x', lw=1)
axs[2].plot(ve.index, ve['y'], 'g', label=r'y', lw=1)
axs[2].plot(ve.index, ve['z'], 'b', label=r'z', lw=1)

# axs[3].plot(pthe.index, -pthe.values, 'r')
# axs[3].plot(PiDe.index, -PiDe.values, 'b', lw=1)
axs[3].plot(PiDe.index, -PiDe.cumsum(), 'b', lw=1)

# axs[4].plot(pthi.index, -pthi.values, 'k')
# axs[5].plot(PiDi.index, -PiDi.values, 'k')

# axs[4].plot(pthe.index, -pthe.values, 'r', lw=1)
axs[4].plot(pthe.index, -pthe.cumsum(), 'r', lw=1)

# axs[4].set_yscale('log')

# axs[1].set_ylim(0, 0.15)
# axs[3].set_ylim(-0.6, 1.2)
axs[4].set_ylim(25, 70)

# start_time = datetime.datetime(2017, 7, 11, 22, 34)
# end_time = datetime.datetime(2017, 7, 11, 22, 34, 4)

start_time = datetime.datetime(2017, 7, 11, 22, 33, 58)
end_time = datetime.datetime(2017, 7, 11, 22, 34, 8)

# start_time = datetime.datetime(2017, 6, 17, 20, 24, 3)
# end_time = datetime.datetime(2017, 6, 17, 20, 24, 11)

axs[0].legend(bbox_to_anchor=(1, 1), loc='upper left')
axs[1].legend(bbox_to_anchor=(1, 1), loc='upper left')
axs[2].legend(bbox_to_anchor=(1, 1), loc='upper left')
# axs[3].legend(bbox_to_anchor=(1.1, 0.1))
# axs[4].legend(bbox_to_anchor=(1.1, 0.1))

axs[0].set_ylabel(r'B[nT]')
axs[1].set_ylabel(r'$v_i$[km/s]')
axs[2].set_ylabel(r'$v_e$[km/s]')
axs[3].set_ylabel(r'Pi-D$_{elc}$')
axs[4].set_ylabel(r'$p\theta_{elc}$')

for ax in axs:
    ax.grid(which='both')

plt.xlim(start_time, end_time)

plt.tight_layout()

plt.savefig('PiD_ptheta_July11_zoomed_out.png', dpi=600)

# axs[4].plot(avgP.index, avgP['xx'])

# plt.show()