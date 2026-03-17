import matplotlib.pyplot as plt
from src.compute_routines.compute_transport import *
from src.compute_routines.compute_transport_errors import *

# fname = r'/home/sroy/Documents/MMS_events/MP_Reconnection/20180415_043241_20180415_043246.h5'
# fname = r'/home/sroy/Documents/MMS_events/Tail_Reconnection/20170617_202403_20170617_202411.h5'
fname = r'/home/sroy/Documents/MMS_events/MP_Reconnection/20160214_204150_20160214_204210.h5'

# Compute kinetic energy transport and error
Ef_transport = compute_kinetic_energy_transport(fname)
Ef_transport_err = compute_Ef_transport_err(fname)

Eth_flux = compute_thermal_energy_flux(fname)
Eth_flux_err = compute_Eth_flux_err(fname)

Eth_transport = compute_thermal_energy_transport(fname)
Eth_transport_err = compute_Eth_transport_err(fname)

PW_transport = compute_pressure_work_transport(fname)
PW_transport_err = compute_pressure_work_transport_err(fname)

S_transport = compute_div_Poynting_flux(fname)
S_transport_err = compute_Poynting_flux_transport_err(fname)

fig, axs = plt.subplots(4, 1)

axs[0].plot(Ef_transport, 'b')
axs[0].plot(Ef_transport_err, 'gray')

axs[1].plot(Eth_transport, 'b')
axs[1].plot(Eth_transport_err, 'gray')

axs[2].plot(PW_transport, 'b')
axs[2].plot(PW_transport_err, 'gray')

axs[3].plot(S_transport, 'b')
axs[3].plot(S_transport_err, 'gray')

plt.show()