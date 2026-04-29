import matplotlib.pyplot as plt
from src.compute_routines.compute_transport import *
from src.compute_routines.compute_transport_errors import *
from datetime import datetime

# fname = r'/home/sroy/Documents/MMS_events/MP_Reconnection/20180415_043241_20180415_043246.h5'
# fname = r'/home/sroy/Documents/MMS_events/Tail_Reconnection/20170617_202403_20170617_202411.h5'
fname = r'/home/sohom/Documents/MMS_events/MP_Reconnection/20160214_204150_20160214_204210.h5'

start_time = datetime(2016, 2, 14, 20, 41, 50)
end_time = datetime(2016, 2, 14, 20, 42)

# Compute kinetic energy transport and error
Ef_flux_ion = compute_kinetic_energy_flux(fname, species='ion', probe=2)
Ef_flux_elc = compute_kinetic_energy_flux(fname, species='elc', probe=2)

Ef_flux_ion_err = compute_Ef_flux_err(fname, species='ion', probe=2)
Ef_flux_elc_err = compute_Ef_flux_err(fname, species='elc', probe=2)

Ef_flux_ion_norm = np.sqrt(Ef_flux_ion['x']**2 + Ef_flux_ion['y']**2 + Ef_flux_ion['z']**2)
Ef_flux_ion_err_norm = np.sqrt(Ef_flux_ion_err['x']**2 + Ef_flux_ion_err['y']**2 + Ef_flux_ion_err['z']**2)

Ef_flux_elc_norm = np.sqrt(Ef_flux_elc['x']**2 + Ef_flux_elc['y']**2 + Ef_flux_elc['z']**2)
Ef_flux_elc_err_norm = np.sqrt(Ef_flux_elc_err['x']**2 + Ef_flux_elc_err['y']**2 + Ef_flux_elc_err['z']**2)

#Compute thermal energy flux and error
Eth_flux_ion = compute_thermal_energy_flux(fname, species='ion', probe=2)
Eth_flux_ion_err = compute_Eth_flux_err(fname, species='ion', probe=2)

Eth_flux_ion_norm = np.sqrt(Eth_flux_ion['x']**2 + Eth_flux_ion['y']**2 + Eth_flux_ion['z']**2)
Eth_flux_ion_err_norm = np.sqrt(Eth_flux_ion_err['x']**2 + Eth_flux_ion_err['y']**2 + Eth_flux_ion_err['z']**2)

Eth_flux_elc = compute_thermal_energy_flux(fname, species='elc', probe=2)
Eth_flux_elc_err = compute_Eth_flux_err(fname, species='elc', probe=2)

Eth_flux_elc_norm = np.sqrt(Eth_flux_elc['x']**2 + Eth_flux_elc['y']**2 + Eth_flux_elc['z']**2)
Eth_flux_elc_err_norm = np.sqrt(Eth_flux_elc_err['x']**2 + Eth_flux_elc_err['y']**2 + Eth_flux_elc_err['z']**2)

plt.plot(Ef_flux_ion_norm, 'r')
plt.plot(Ef_flux_ion_err_norm, 'orange')

plt.plot(Ef_flux_elc_norm, 'b')
plt.plot(Ef_flux_elc_err_norm, 'gray')

# plt.plot(Eth_flux_elc_norm, 'b')
# plt.plot(Eth_flux_elc_err_norm, 'gray')

plt.xlim(start_time, end_time)

plt.show()