import matplotlib.pyplot as plt
from src.compute_routines.compute_transport import *
from src.compute_routines.compute_transport_errors import *

fname = r'/home/sroy/Documents/MMS_events/MP_Reconnection/20180415_043241_20180415_043246.h5'

# Compute kinetic energy transport and error
KE_transport = compute_kinetic_energy_transport(fname)
KE_transport_err = compute_Ef_transport_err(fname)

plt.plot(KE_transport, 'b')
plt.plot(KE_transport_err, 'gray')

plt.show()