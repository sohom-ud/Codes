from src.compute_routines.compute_transport import *
from src.compute_routines.compute_transport_errors import *
from src.compute_routines.compute_PiD_functions import *
from src.compute_routines.compute_PiD_errors import *
from src.utils.hdf_to_df import hdf_to_df
from datetime import datetime
import matplotlib.pyplot as plt
# from src.compute_routines.compute_MVA import MVA

plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15

fname = r'/home/sroy/Documents/MMS_events/MP_Reconnection/20160214_204150_20160214_204210.h5'

df_dict = hdf_to_df(fname)

v1 = df_dict['v_spincorr_ion_1']
B1 = df_dict['b_gse_2']

B1 = B1.drop('mag', axis=1)

# v_xl_LMN = [-137, -142, -3] # X-line velocity in km/s (obtained from Eastwood 2020 PRL)

#LMN coordinates obtained from Fargette 2024 paper

L = np.array([0.09, -0.02, 1])
M = np.array([-0.41, -0.91, 0.02])
N = np.array([0.91, -0.41, -0.09])

# L = np.array([0.36407441, 0.071052083, 0.92865571])
# M = np.array([-0.022889708, -0.99610209, 0.085186237])
# N = np.array([0.93108855, -0.052270787, -0.36102892])

R = np.array([L, M, N]) #Transformation matrix

B_LMN = pd.DataFrame(B1.values @ R, index=B1.index, columns=['L', 'M', 'N'])

# v_xl_XYZ = R.T @ v_xl_LMN

EDR_time = datetime(2016, 2, 14, 20, 41, 56)

EDR_start = datetime(2016, 2, 14, 20, 41, 56, 80000)
EDR_end = datetime(2016, 2, 14, 20, 41, 56, 340000)

v_xl_X = np.mean(v1['x'][EDR_start:EDR_end])
v_xl_Y = np.mean(v1['y'][EDR_start:EDR_end])
v_xl_Z = np.mean(v1['z'][EDR_start:EDR_end])

v_xl_XYZ = [v_xl_X, v_xl_Y, v_xl_Z]

F_K_ion = compute_kinetic_energy_transport(fname, species='ion', reselectron=True, v_xl=v_xl_XYZ)*1e9
F_K_elc = compute_kinetic_energy_transport(fname, species='elc', reselectron=True, v_xl=v_xl_XYZ)*1e9

F_T_ion = compute_thermal_energy_transport(fname, species='ion', reselectron=True, v_xl=v_xl_XYZ)*1e9
F_T_elc = compute_thermal_energy_transport(fname, species='elc', reselectron=True, v_xl=v_xl_XYZ)*1e9

F_P_ion = compute_pressure_work_transport(fname, species='ion', reselectron=True, v_xl=v_xl_XYZ)*1e9
F_P_elc = compute_pressure_work_transport(fname, species='elc', reselectron=True, v_xl=v_xl_XYZ)*1e9

divS = compute_div_Poynting_flux(fname)*1e9

F_K_ion_err = compute_Ef_transport_err(fname, species='ion', reselectron=True, v_xl=v_xl_XYZ)*1e9
F_K_elc_err = compute_Ef_transport_err(fname, species='elc', reselectron=True, v_xl=v_xl_XYZ)*1e9

F_T_ion_err = compute_Eth_transport_err(fname, species='ion', reselectron=True, v_xl=v_xl_XYZ)*1e9
F_T_elc_err = compute_Eth_transport_err(fname, species='elc', reselectron=True, v_xl=v_xl_XYZ)*1e9

F_P_ion_err = compute_pressure_work_transport_err(fname, species='ion', reselectron=True, v_xl=v_xl_XYZ)*1e9
F_P_elc_err = compute_pressure_work_transport_err(fname, species='elc', reselectron=True, v_xl=v_xl_XYZ)*1e9

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

xlim_start = datetime(2016, 2, 14, 20, 41, 55)
xlim_end = datetime(2016, 2, 14, 20, 41, 57, 500000)

fig, axs = plt.subplots(5, 1, figsize=(8,10), sharex=True)
start_time = datetime(2017, 6, 8, 11, 20, 42)
end_time = datetime(2015, 12, 8, 11, 20, 45)

axs[0].plot(B_LMN.index, B_LMN['L'], 'r', label=r'$B_L$')
axs[0].plot(B_LMN.index, B_LMN['M'], 'g', label=r'$B_M$')
axs[0].plot(B_LMN.index, B_LMN['N'], 'b', label=r'$B_N$')

axs[1].plot(F_K_elc, 'k', label=r'$\nabla\cdot\mathbf{F}_{K, e}$')
axs[1].fill_between(F_K_elc.index, -F_K_elc_err, F_K_elc_err, color='gray', alpha=0.5)
axs[1].axhline(0, color='r', ls='--')

axs[2].plot(F_T_elc, 'k', label=r'$\nabla\cdot\mathbf{F}_{T, e}$')
axs[2].fill_between(F_T_elc.index, -F_T_elc_err, F_T_elc_err, color='gray', alpha=0.5)
axs[2].axhline(0, color='r', ls='--')

axs[3].plot(F_P_elc, 'k', label=r'$\nabla\cdot\mathbf{F}_{P, e}$')
axs[3].fill_between(F_P_elc_err.index, -F_P_elc_err, F_P_elc_err, color='gray', alpha=0.5)
axs[3].axhline(0, color='r', ls='--')

# axs[4].plot(h_transport, 'k', label=r'$\nabla\cdot\mathbf{q}_e$')
# axs[4].fill_between(h_transport.index, -h_transport_err, h_transport_err, color='gray', alpha=0.5)
# axs[5].axhline(0, color='r', ls='--')

axs[4].plot(divS, 'k', label=r'$\nabla\cdot\mathbf{S} $')
axs[4].fill_between(divS.index, -divS_err, divS_err, color='gray', alpha=0.5)
axs[4].axhline(0, color='r', ls='--')

# plt.xlim(start_time, end_time)

plt.xlabel(r'Datetime(UTC)', fontsize=15)

axs[0].set_ylabel(r'[nT]', fontsize=15)

for ax in axs:
    ax.legend(fontsize=15, loc='upper left')
    ax.grid(which='both')

for ax in axs[1:]:
    ax.set_ylabel(r'[nW/m$^3$]', fontsize=15)

plt.xlim(xlim_start, xlim_end)

plt.subplots_adjust(hspace=0.05)

plt.savefig(r'transport_errs_Lenouvel_2021_magnetotail_elc.png', dpi=600)