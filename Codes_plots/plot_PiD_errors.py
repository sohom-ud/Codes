import matplotlib.pyplot as plt
from src.compute_routines.compute_transport import *
from src.compute_routines.compute_transport_errors import *
from src.compute_routines.compute_PiD_functions import *
from src.compute_routines.compute_PiD_errors import *
from src.utils.resample import resample
from datetime import datetime

# fname = r'/home/sroy/Documents/MMS_events/Tail_Reconnection/20170617_202403_20170617_202411.h5' # Lu(2020)
fname = r'/home/sroy/Documents/MMS_events/Tail_Reconnection/20170810_121800_20170810_121900.h5' # Hubbert(2022)
# fname = r'/home/sroy/Documents/MMS_events/MP_Reconnection/20160214_204150_20160214_204210.h5'
# fname = r'/home/sroy/Documents/MMS_events/MP_Reconnection/20151208_112042_20151208_112045.h5' # Burch(2016)
# fname = r'/home/sroy/Documents/MMS_events/MSH_Reconnection/20151130_002340_20151130_002415.h5' # Voros(2017)

plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15

df_dict = hdf_to_df(fname)

B1 = df_dict['b_gse_1']
ni1 = df_dict['N_ion_1']
ne1 = df_dict['N_elc_1']

species='elc'

PiD = compute_PiD(fname, species=species)
ptheta = compute_ptheta(fname, species=species)
PS = compute_PS(fname, species=species)

PiD_err = compute_PiD_err(fname, species=species)
ptheta_err = compute_ptheta_err(fname, species=species)
PS_err = compute_PS_err(fname, species=species)

fig, axs = plt.subplots(4, 1, figsize=(10, 8), sharex=True)
# start_time = datetime(2015, 11, 30, 0, 23, 50)
# end_time = datetime(2015, 11, 30, 0, 23, 58)

axs[0].plot(B1.index, B1['x'], 'r', label=r'$x$')
axs[0].plot(B1.index, B1['y'], 'g', label=r'$y$')
axs[0].plot(B1.index, B1['z'], 'b', label=r'$z$')

axs[1].plot(ni1.index, ni1, 'k', label='ion')
# axs[1].plot(ni1.index, resample(ne1, ni1), 'b', label='elc')
axs[1].plot(ne1.index, ne1, 'b', label='elc')

axs[2].plot(-PiD, 'k')
axs[2].fill_between(PiD.index, -PiD_err, PiD_err, color='gray', alpha=0.5)
axs[2].axhline(0, color='r', ls='--')

axs[3].plot(-ptheta, 'k')
axs[3].fill_between(ptheta.index, -ptheta_err, ptheta_err, color='gray', alpha=0.5)
axs[3].axhline(0, color='r', ls='--')

# plt.xlim(start_time, end_time)

plt.xlabel(r'Datetime(UTC)', fontsize=15)

axs[0].set_ylabel(r'B[nT]', fontsize=15)
axs[1].set_ylabel(r'n[cm$^-3$]', fontsize=15)
axs[2].set_ylabel(r'-$\Pi_{ij}D_{ij}$[nW/m$^3$]', fontsize=15)
axs[3].set_ylabel(r'-$p\theta$[nW/m$^3$]', fontsize=15)

axs[0].legend(fontsize=15, loc='upper left')
axs[1].legend(fontsize=15, loc='upper left')

for ax in axs:
    ax.grid(which='both')

# plt.suptitle('Ions', fontsize=25)
plt.suptitle('Electrons', fontsize=25)

plt.subplots_adjust(hspace=0.05)
plt.tight_layout()

# plt.savefig(rf'/home/sroy/Documents/Codes/Plots/PiD_errors/Voros2017_magnetosheath_{species}.png', dpi=600)
# plt.savefig(rf'/home/sroy/Documents/Codes/Plots/PiD_errors/Burch2016_magnetopause_{species}.png', dpi=600)
# plt.savefig(rf'/home/sroy/Documents/Codes/Plots/PiD_errors/Lu2020_magnetotail_{species}.png', dpi=600)
plt.savefig(rf'/home/sroy/Documents/Codes/Plots/PiD_errors/Hubbert2022_magnetotail_{species}.png', dpi=600)