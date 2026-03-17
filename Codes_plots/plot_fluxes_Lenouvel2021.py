import matplotlib.pyplot as plt
from src.compute_routines.compute_transport import *
from src.compute_routines.compute_transport_errors import *

# fname = r'/home/sroy/Documents/MMS_events/MP_Reconnection/20180415_043241_20180415_043246.h5'
# fname = r'/home/sroy/Documents/MMS_events/Tail_Reconnection/20170617_202403_20170617_202411.h5'
fname = r'/home/sroy/Documents/MMS_events/MP_Reconnection/20160214_204150_20160214_204210.h5'

# Compute kinetic energy transport and error
Ef_flux = compute_kinetic_energy_flux(fname, species='elc', probe=2)
Ef_flux_err = compute_Ef_flux_err(fname, species='elc', probe=2)

plt.plot(Ef_flux, 'b')
plt.plot(Ef_flux_err, 'gray')

plt.show()