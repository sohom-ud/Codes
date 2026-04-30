from src.compute_routines.compute_transport import *
from src.compute_routines.compute_transport_errors import *
from src.compute_routines.compute_PiD_functions import *
from src.compute_routines.compute_PiD_errors import *
from src.utils.hdf_to_df import hdf_to_df
from datetime import datetime
import matplotlib.pyplot as plt

fname = r'/home/sroy/Documents/MMS_events/Tail_Reconnection/20170711_223300_20170711_223500.h5'

df_dict = hdf_to_df(fname)

v1 = df_dict['v_spincorr_ion_1']

EDR_start = datetime(2017, 7, 11, 22, 34)
EDR_end = datetime(2017, 7, 11, 22, 34, 6)

v_xl_X = np.mean(v1['x'][EDR_start:EDR_end])
v_xl_Y = np.mean(v1['y'][EDR_start:EDR_end])
v_xl_Z = np.mean(v1['z'][EDR_start:EDR_end])

v_xl_XYZ = [v_xl_X, v_xl_Y, v_xl_Z]

F_K_ion = compute_kinetic_energy_transport(fname, species='ion', reselectron=True, v_xl=v_xl_XYZ)
F_K_elc = compute_kinetic_energy_transport(fname, species='elc', reselectron=True, v_xl=v_xl_XYZ)

F_T_ion = compute_thermal_energy_transport(fname, species='ion', reselectron=True, v_xl=v_xl_XYZ)
F_T_elc = compute_thermal_energy_transport(fname, species='elc', reselectron=True, v_xl=v_xl_XYZ)

F_P_ion = compute_pressure_work_transport(fname, species='ion', reselectron=True, v_xl=v_xl_XYZ)
F_P_elc = compute_pressure_work_transport(fname, species='elc', reselectron=True, v_xl=v_xl_XYZ)

divS = compute_div_Poynting_flux(fname)

F_K_ion_err = compute_Ef_transport_err(fname, species='ion', reselectron=True, v_xl=v_xl_XYZ)
F_K_elc_err = compute_Ef_transport_err(fname, species='elc', reselectron=True, v_xl=v_xl_XYZ)

F_T_ion_err = compute_Eth_transport_err(fname, species='ion', reselectron=True, v_xl=v_xl_XYZ)
F_T_elc_err = compute_Eth_transport_err(fname, species='elc', reselectron=True, v_xl=v_xl_XYZ)

F_P_ion_err = compute_pressure_work_transport_err(fname, species='ion', reselectron=True, v_xl=v_xl_XYZ)
F_P_elc_err = compute_pressure_work_transport_err(fname, species='elc', reselectron=True, v_xl=v_xl_XYZ)

divS_err = compute_Poynting_flux_transport_err(fname, res='e')

PSi = compute_PS(fname, species='ion')
PSe = compute_PS(fname, species='elc')

PSi_err = compute_PS_err(fname, species='ion')
PSe_err = compute_PS_err(fname, species='elc')

jdotE = compute_jdotE(fname)

jpart = compute_jpart(fname)
jcurl = compute_jcurl(fname)

jpartmag = np.sqrt(jpart['x']**2 + jpart['y']**2 + jpart['z']**2)
jcurlmag = np.sqrt(jcurl['x']**2 + jcurl['y']**2 + jcurl['z']**2)

xlim_start = datetime(2017, 7, 11, 22, 33, 50)
xlim_end = datetime(2017, 7, 11, 22, 34, 10)

fig, axs = plt.subplots(7, 1, sharex=True)

axs[0].plot(F_K_elc, 'k')
axs[0].fill_between(F_K_elc.index, -F_K_elc_err, F_K_elc_err, color='gray')

axs[1].plot(F_T_elc, 'k')
axs[1].fill_between(F_T_elc.index, -F_T_elc_err, F_T_elc_err, color='gray')

axs[2].plot(F_P_elc, 'k')
axs[2].fill_between(F_P_elc.index, -F_P_elc_err, F_P_elc_err, color='gray')

axs[3].plot(divS, 'k')
axs[3].fill_between(divS.index, -divS_err, divS_err, color='gray')

# axs[4].plot(-PSi, 'k')
# axs[4].fill_between(PSi.index, -PSi_err, PSi_err, color='gray')

axs[4].plot(-PSe, 'k')
axs[4].fill_between(PSe.index, -PSe_err, PSe_err, color='gray')

# axs[6].plot(jcurlmag, 'k')
# axs[6].plot(jpartmag, 'b')
# axs[6].fill_between(F_K_ion.index, -F_K_ion_err, F_K_ion_err, color='gray')

plt.xlim(xlim_start, xlim_end)

plt.show()