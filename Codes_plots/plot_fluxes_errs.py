import matplotlib.pyplot as plt
from src.compute_routines.compute_transport import *
from src.compute_routines.compute_transport_errors import *

def norm(A):
    return np.sqrt(A['x']**2 + A['y']**2 + A['z']**2)

fname = r'/home/sroy/Documents/MMS_events/MP_Reconnection/20180415_043241_20180415_043246.h5'
# fname = r'/home/sroy/Documents/MMS_events/Tail_Reconnection/20170617_202403_20170617_202411.h5'
# fname = r'/home/sroy/Documents/MMS_events/MP_Reconnection/20160214_204150_20160214_204210.h5'

plt.rcParams['text.usetex'] = True

# Compute kinetic energy transport and error
Ef_flux = compute_kinetic_energy_flux(fname)
Ef_flux_err = compute_Ef_flux_err(fname)

Eth_flux = compute_thermal_energy_flux(fname)
Eth_flux_err = compute_Eth_flux_err(fname)

PW = compute_pressure_work(fname)
PW_err = compute_pressure_work_err(fname)

S = compute_Poynting_flux(fname)
S_err = compute_Poynting_flux_err(fname)

#Calculating norms of fluxes and errors
Ef_flux_norm = norm(Ef_flux)
Ef_flux_err_norm = norm(Ef_flux_err)

Eth_flux_norm = norm(Eth_flux)
Eth_flux_err_norm = norm(Eth_flux_err)

PW_norm = norm(PW)
PW_err_norm = norm(PW_err)

S_norm = norm(S)
S_err_norm = norm(S_err)

fig, axs = plt.subplots(4, 1, sharex=True)

# axs[0].plot(Ef_transport, 'b')
# axs[0].plot(Ef_transport_err, 'gray')

axs[0].plot(Ef_flux_norm * 1e3, 'b', label=r'$|E^f \mathbf{u}|$')
axs[0].fill_between(Ef_flux_norm.index, (Ef_flux_norm - Ef_flux_err_norm)*1e3, (Ef_flux_norm + Ef_flux_err_norm)*1e3, color='gray', alpha=0.5)

axs[1].plot(Eth_flux_norm * 1e3, 'b', label=r'$|E^{th} \mathbf{u}|$')
axs[1].fill_between(Eth_flux_norm.index, (Eth_flux_norm - Eth_flux_err_norm)*1e3, (Eth_flux_norm + Eth_flux_err_norm)*1e3, color='gray', alpha=0.5)

axs[2].plot(PW_norm * 1e3, 'b', label=r'$|\mathbf{P}\cdot\mathbf{u}|$')
axs[2].fill_between(PW_norm.index, (PW_norm - PW_err_norm)*1e3, (PW_norm + PW_err_norm)*1e3, color='gray', alpha=0.5)

axs[3].plot(S_norm * 1e3, 'b', label=r'$|\mathbf{S}|$')
axs[3].fill_between(S_norm.index, (S_norm - S_err_norm)*1e3, (S_norm + S_err_norm)*1e3, color='gray', alpha=0.7)

# axs[0].set_ylabel(, fontsize=15)
for ax in axs:
    ax.legend(fontsize=15)
    ax.set_ylabel(r'[mW/m$^2]$', fontsize=12)

plt.tight_layout()

plt.savefig(r'/home/sroy/Documents/Codes/Plots/Flux_transport_errors/20180415_flux_err.png')
# plt.show()