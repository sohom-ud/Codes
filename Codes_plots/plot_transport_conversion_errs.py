import matplotlib.pyplot as plt
from src.compute_routines.compute_transport import *
from src.compute_routines.compute_transport_errors import *
from src.compute_routines.compute_PiD_functions import *
from src.compute_routines.compute_PiD_errors import *
from datetime import datetime

# fname = r'/home/sroy/Documents/MMS_events/MP_Reconnection/20180415_043241_20180415_043246.h5'
# fname = r'/home/sroy/Documents/MMS_events/Tail_Reconnection/20170617_202403_20170617_202411.h5'
# fname = r'/home/sroy/Documents/MMS_events/MP_Reconnection/20160214_204150_20160214_204210.h5'
# fname = r'/home/sroy/Documents/MMS_events/MP_Reconnection/20151208_112042_20151208_112045.h5'
fname = r'/home/sroy/Documents/MMS_events/MSH_Reconnection/20151130_002340_20151130_002415.h5'

plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15

df_dict = hdf_to_df(fname)

B1 = df_dict['b_gse_1']

species='ion'

# Compute kinetic energy transport and error
Ef_transport = compute_kinetic_energy_transport(fname, species=species)*1e9
Ef_transport_err = compute_Ef_transport_err(fname, species=species)*1e9

Eth_transport = compute_thermal_energy_transport(fname, species=species)*1e9
Eth_transport_err = compute_Eth_transport_err(fname, species=species)*1e9

PW_transport = compute_pressure_work_transport(fname, species=species)*1e9
PW_transport_err = compute_pressure_work_transport_err(fname, species=species)*1e9

S_transport = compute_div_Poynting_flux(fname, res='e')*1e9
S_transport_err = compute_Poynting_flux_transport_err(fname, res='e')*1e9

h_transport = compute_div_q(fname, species=species)*1e9
h_transport_err = compute_heatflux_transport_err(fname, species=species)*1e9

fig, axs = plt.subplots(6, 1, figsize=(10, 12), sharex=True)
start_time = datetime(2015, 11, 30, 0, 23, 50)
end_time = datetime(2015, 11, 30, 0, 23, 58)

axs[0].plot(B1.index, B1['x'], 'r', label=r'$B_x$(GSE)')
axs[0].plot(B1.index, B1['y'], 'g', label=r'$B_y$(GSE)')
axs[0].plot(B1.index, B1['z'], 'b', label=r'$B_z$(GSE)')

axs[1].plot(Ef_transport, 'k', label=r'$\nabla\cdot\mathbf{F}_{K,  ' + species[0] + '}$')
axs[1].fill_between(Ef_transport.index, -Ef_transport_err, Ef_transport_err, color='gray', alpha=0.5)
axs[1].axhline(0, color='r', ls='--')

axs[2].plot(Eth_transport, 'k', label=r'$\nabla\cdot\mathbf{F}_{T,  ' + species[0] + '}$')
axs[2].fill_between(Eth_transport.index, -Eth_transport_err, Eth_transport_err, color='gray', alpha=0.5)
axs[2].axhline(0, color='r', ls='--')

axs[3].plot(PW_transport, 'k', label=r'$\nabla\cdot\mathbf{F}_{P, ' + species[0] + '}$')
axs[3].fill_between(PW_transport.index, -PW_transport_err, PW_transport_err, color='gray', alpha=0.5)
axs[3].axhline(0, color='r', ls='--')

axs[4].plot(h_transport, 'k', label=r'$\nabla\cdot\mathbf{q}_' + species[0] + '$')
axs[4].fill_between(h_transport.index, -h_transport_err, h_transport_err, color='gray', alpha=0.5)
axs[5].axhline(0, color='r', ls='--')

axs[5].plot(S_transport, 'k', label=r'$\nabla\cdot\mathbf{S} $')
axs[5].fill_between(S_transport.index, -S_transport_err, S_transport_err, color='gray', alpha=0.5)
axs[5].axhline(0, color='r', ls='--')

plt.xlim(start_time, end_time)

plt.xlabel(r'Datetime(UTC)', fontsize=15)

axs[0].set_ylabel(r'[nT]', fontsize=15)

for ax in axs:
    ax.legend(fontsize=15, loc='upper left')
    ax.grid(which='both')

for ax in axs[1:]:
    ax.set_ylabel(r'[nW/m$^3$]', fontsize=15)

plt.subplots_adjust(hspace=0.05)
plt.tight_layout()

plt.savefig(rf'/home/sroy/Documents/Codes/Plots/Flux_transport_errors/Voros2017_magnetosheath_{species}.png', dpi=600)