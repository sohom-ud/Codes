import matplotlib.pyplot as plt
from src.compute_routines.compute_transport import *
from src.compute_routines.compute_transport_errors import *
from src.compute_routines.compute_PiD_functions import *
from src.compute_routines.compute_PiD_errors import *

# fname = r'/home/sroy/Documents/MMS_events/Tail_Reconnection/20170617_202403_20170617_202411.h5'
fname = r'/home/sroy/Documents/MMS_events/MSH_Reconnection/20151130_002340_20151130_002415.h5'

species='ion'

Ef_transport = compute_kinetic_energy_transport(fname, species=species)*1e9
Ef_transport_err = compute_Ef_transport_err(fname, species=species)*1e9

PS = compute_PS(fname, species=species)*1e9
PS_err = compute_PS_err(fname, species=species)*1e9

PS_NSR = PS_err**2/PS['PS_ion']**2

Eth_transport = compute_thermal_energy_transport(fname, species=species)*1e9
Eth_transport_err = compute_Eth_transport_err(fname, species=species)*1e9

Ef_NSR = Ef_transport_err**2 / Ef_transport**2
Eth_NSR = Eth_transport_err**2 / Eth_transport**2

# Ef_NSR = Ef_NSR[np.logical_and(Ef_NSR>=0, Ef_NSR<=5)]

fig, axs = plt.subplots(2, 1)

axs[0].plot(Ef_transport[Ef_NSR<1], 'b.')
axs[1].plot(Ef_NSR, color='k')

axs[1].axhline(1, color='r')

plt.show()

fig, axs = plt.subplots(2, 1)

axs[0].plot(-PS[PS_NSR<1], 'b.')
axs[1].plot(PS_err, 'k')

plt.show()

# axs[0].hist(np.log10(Ef_NSR), bins=100)
# axs[1].hist(np.log10(Eth_NSR), bins=100)

# plt.show()