from src.compute_routines.compute_transport import *
from src.compute_routines.compute_transport_errors import *
from src.compute_routines.compute_PiD_functions import *
from src.compute_routines.compute_PiD_errors import *
from src.utils.hdf_to_df import hdf_to_df
from datetime import datetime
import matplotlib.pyplot as plt

fname = r'/home/sroy/Documents/MMS_events/MP_Reconnection/20160214_204150_20160214_204210.h5'

df_dict = hdf_to_df(fname)

v1 = df_dict['v_spincorr_ion_1']

# v_xl_LMN = [-137, -142, -3] # X-line velocity in km/s (obtained from Eastwood 2020 PRL)

#LMN coordinates obtained from Burch 2016 paper
# L = np.array([0.3665, -0.1201, 0.9226])
# M = np.array([0.5694, -0.7553, -0.3245])
# N = np.array([0.7358, 0.6443, -0.2084])

# L = np.array([0.36407441, 0.071052083, 0.92865571])
# M = np.array([-0.022889708, -0.99610209, 0.085186237])
# N = np.array([0.93108855, -0.052270787, -0.36102892])

# R = np.array([L, M, N]) #Transformation matrix

# v_xl_XYZ = R.T @ v_xl_LMN

EDR_time = datetime(2016, 2, 14, 20, 41, 56)

EDR_start = datetime(2016, 2, 14, 20, 41, 56, 80000)
EDR_end = datetime(2016, 2, 14, 20, 41, 56, 340000)

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

xlim_start = datetime(2016, 2, 14, 20, 41, 54)
xlim_end = datetime(2016, 2, 14, 20, 41, 58)

fig, axs = plt.subplots(7, 1, sharex=True)

axs[0].plot(F_K_ion, 'k')
axs[0].fill_between(F_K_ion.index, -F_K_ion_err, F_K_ion_err, color='gray')

axs[1].plot(F_T_ion, 'k')
axs[1].fill_between(F_T_ion.index, -F_T_ion_err, F_T_ion_err, color='gray')

axs[2].plot(F_P_ion, 'k')
axs[2].fill_between(F_P_ion.index, -F_P_ion_err, F_P_ion_err, color='gray')

axs[3].plot(divS, 'k')
axs[3].fill_between(divS.index, -divS_err, divS_err, color='gray')

axs[4].plot(-PSi, 'k')
axs[4].fill_between(PSi.index, -PSi_err, PSi_err, color='gray')

axs[5].plot(-PSe, 'k')
axs[5].fill_between(PSe.index, -PSe_err, PSe_err, color='gray')

axs[6].plot(jdotE, 'k')
# axs[6].fill_between(F_K_ion.index, -F_K_ion_err, F_K_ion_err, color='gray')

plt.xlim(xlim_start, xlim_end)

plt.show()