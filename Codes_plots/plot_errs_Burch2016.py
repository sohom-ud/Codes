import matplotlib.pyplot as plt
from src.compute_routines.compute_transport import *
from src.compute_routines.compute_transport_errors import *
from src.compute_routines.compute_PiD_functions import *
from src.compute_routines.compute_PiD_errors import *

fname = r'/home/sroy/Documents/MMS_events/MP_Reconnection/20180415_043241_20180415_043246.h5'
# fname = r'/home/sroy/Documents/MMS_events/Tail_Reconnection/20170617_202403_20170617_202411.h5'
# fname = r'/home/sroy/Documents/MMS_events/MP_Reconnection/20160214_204150_20160214_204210.h5'

df_dict = hdf_to_df(fname)

# Compute kinetic energy transport and error
Ef_transport = compute_kinetic_energy_transport(fname, species='elc')
Ef_transport_err = compute_Ef_transport_err(fname, species='elc')

Ef_flux = compute_kinetic_energy_flux(fname)
Ef_flux_err = compute_Ef_flux_err(fname)

Ef_flux_norm = np.sqrt(Ef_flux['x']**2 + Ef_flux['y']**2 + Ef_flux['z']**2)

Ef_flux_norm_err = np.sqrt((Ef_flux['x']/Ef_flux_norm)**2*Ef_flux_err['x']**2 + (Ef_flux['y']/Ef_flux_norm)**2*Ef_flux_err['y']**2 + (Ef_flux['z']/Ef_flux_norm)**2*Ef_flux_err['z']**2)

# Eth_flux = compute_thermal_energy_flux(fname)
# Eth_flux_err = compute_Eth_flux_err(fname)

Eth_transport = compute_thermal_energy_transport(fname, species='elc')
Eth_transport_err = compute_Eth_transport_err(fname, species='elc')

PW_transport = compute_pressure_work_transport(fname, species='elc')
PW_transport_err = compute_pressure_work_transport_err(fname, species='elc')

# S_transport = compute_div_Poynting_flux(fname)
# S_transport_err = compute_Poynting_flux_transport_err(fname)

PS = compute_PS(fname, species='elc')
PS_err = compute_PS_err(fname, species='elc')

fig, axs = plt.subplots(4, 1, sharex=True)

# axs[0].plot(Ef_flux_norm, 'b')
# axs[0].fill_between(Ef_flux_norm.index, -Ef_flux_norm_err, Ef_flux_norm_err, color='gray', alpha=0.5)

axs[0].plot(Ef_transport, 'b')
# axs[1].plot(Ef_transport_err, 'gray')
# axs[1].plot(-Ef_transport_err, 'gray')
axs[0].fill_between(Ef_transport.index, -Ef_transport_err, Ef_transport_err, color='gray', alpha=0.5)

axs[1].plot(Eth_transport, 'b')
axs[1].fill_between(Eth_transport.index, -Eth_transport_err, Eth_transport_err, color='gray', alpha=0.5)

axs[2].plot(PW_transport, 'b')
axs[2].fill_between(PW_transport.index, -PW_transport_err, PW_transport_err, color='gray', alpha=0.5)

axs[3].plot(PS, 'b')
axs[3].plot(PS_err, color='gray', alpha=0.5)
# axs[3].plot(S_transport, 'b')
# axs[3].plot(S_transport_err, 'gray')

plt.show()