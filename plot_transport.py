from src.compute_routines.compute_transport import *
from src.compute_routines.compute_transport_errors import *
from src.compute_routines.compute_PiD_functions import *
from src.compute_routines.compute_PiD_errors import *
import matplotlib.pyplot as plt
from pyspedas.projects.mms import fpi
from datetime import datetime
from src.utils.hdf_to_df import hdf_to_df
from src.utils.methods import dot

kB = 1.380649e-23
mu_0 = 4*np.pi*1e-7

##----EVENT 1----#

# fname = r'/home/sroy/Documents/MMS_events/MSH_Turbulence/20180423_075010_20180423_075240.h5'

# start_time = datetime(2018, 4, 23, 7, 51, 55)
# end_time = datetime(2018, 4, 23, 7, 52, 5)

# trange = ['2018-04-23/07:50:10', '2018-04-23/07:52:40']

#----EVENT 2----#

# fname = r'/home/sroy/Documents/MMS_events/MSH_Turbulence/20170928_063150_20170928_063330.h5'

# start_time = datetime(2017, 9, 28, 6, 31, 50)
# end_time = datetime(2017, 9, 28, 6, 33, 30)

# trange = ['2017-09-28/06:31:50', '2017-09-28/06:33:30']

##----EVENT 3----#

interval = '20171226_061243_20171226_065223'
fname = rf'/home/sroy/Documents/MMS_events/MSH_Turbulence/{interval}.h5'

# interval = '20151130_002400_20151130_002700'
# fname = rf'/home/sroy/Documents/MMS_events/MSH_Turbulence/{interval}.h5'

start_time = datetime(2017, 12, 26, 6, 12, 43)
end_time = datetime(2017, 12, 26, 6, 52, 23)

# start_time = datetime(2017, 12, 26, 6, 14, 20)
# end_time = datetime(2017, 12, 26, 6, 15, 00)

# start_time = datetime(2015, 11, 30, 0, 25, 30)
# end_time = datetime(2015, 11, 30, 0, 26, 30)

# trange = ['2015-11-30/00:24:00', '2015-11-30/00:27:00']
trange = ['2017-12-26/06:12:43', '2017-12-26/06:52:23']

probe = 1

df_dict = hdf_to_df(fname)

data = fpi(
    trange=trange,
    probe=probe,
    data_rate='brst',
    datatype=f"des-moms",
    center_measurement=True,
    notplot=True
)
B = df_dict['b_gse_1']
Bmag = B['mag']
B = B.drop('mag', axis=1)

ni = df_dict['N_ion_1']
ne = df_dict['N_elc_1']

vi1 = df_dict['v_spincorr_ion_1']
ve1 = df_dict['v_spincorr_elc_1']

vi1_mean = vi1.mean()
ve1_mean = ve1.mean()

Te_para = pd.DataFrame(data=data['mms1_des_temppara_brst']['y'], index=data['mms1_des_temppara_brst']['x'])

Te_perp = pd.DataFrame(data=data['mms1_des_tempperp_brst']['y'], index=data['mms1_des_tempperp_brst']['x'])

Te_perp = Te_perp[start_time:end_time]
Te_para = Te_para[start_time:end_time]

Te_para = resample(Te_para, ne)
Bmag = resample(Bmag, ne)

beta = ne.mul(Te_para[0], axis=0).div(Bmag**2, axis=0) * 2 * mu_0 * 1.6*1e6*1e-19*1e18
# PiDi = compute_PiD(fname, species='ion')
PiDe = compute_PiD(fname, species='elc')

# jdotEprime = compute_jdotEprime(fname)

# F_K_e = compute_kinetic_energy_flux(fname, species='elc', v_xl=vi1_mean.values)

# F_T_e = compute_thermal_energy_flux(fname, species='elc', v_xl=vi1_mean.values)

# F_P_e = compute_pressure_work(fname, species='elc', v_xl=vi1_mean.values)

# S = compute_Poynting_flux(fname, probe=1, res='B')

A = np.abs(Te_perp/Te_para-1)

# SdotB = dot(S, B)

# PiDi_err = compute_PiD_err(fname, species='ion')
# PiDe_err = compute_PiD_err(fname, species='elc')

# J_K_e_err = compute_Ef_transport_err(fname, species='elc', v_xl=vi1_mean.values)
# J_T_e_err = compute_Eth_transport_err(fname, species='elc', v_xl=vi1_mean.values)
# J_P_e_err = compute_pressure_work_transport_err(fname, species='elc', v_xl=vi1_mean.values)

# J_K_e_err = compute_Ef_transport_err(fname, species='elc', v_xl=vi1_mean.values)
# J_T_e_err = compute_Eth_transport_err(fname, species='elc', v_xl=vi1_mean.values)
# J_P_e_err = compute_pressure_work_transport_err(fname, species='elc', v_xl=vi1_mean.values)

# fig, axs = plt.subplots(9, 1, figsize=(10, 18), sharex=True)

# axs[0].plot(B['x'], 'r', label='x')
# axs[0].plot(B['y'], 'g', label='y')
# axs[0].plot(B['z'], 'b', label='z')

# axs[1].plot(ni, 'r', label='ion')
# axs[1].plot(ne, 'b', label='elc')

# axs[2].plot(-PiDi, 'k')
# axs[3].plot(-PiDe, 'r')

# # axs[0].fill_between(PiDi.index, -PiDi_err, PiDi_err, color='gray', alpha=0.5)
# # axs[1].fill_between(PiDe.index, -PiDe_err, PiDe_err, color='gray', alpha=0.5)

# axs[4].plot(A, 'k')
# axs[4].axhline(1, color='r')

# # axs[3].plot(J_K_i, 'k')
# axs[5].plot(F_K_e['x']*1e6, 'r', label='x')
# axs[5].plot(F_K_e['y']*1e6, 'g', label='y')
# axs[5].plot(F_K_e['z']*1e6, 'b', label='z')
# # # axs[3].fill_between(F_K_e.index, -F_K_e_err*1e9, F_K_e_err*1e9, color='gray', alpha=0.5)

# # # axs[4].plot(J_T_i, 'k')
# axs[6].plot(F_T_e['x']*1e6, 'r', label='x')
# axs[6].plot(F_T_e['y']*1e6, 'g', label='y')
# axs[6].plot(F_T_e['z']*1e6, 'b', label='z')
# # # axs[4].fill_between(J_T_e.index, -J_T_e_err*1e9, J_T_e_err*1e9, color='gray', alpha=0.5)

# # # axs[5].plot(J_P_i, 'k')
# axs[7].plot(F_P_e['x']*1e6, 'r', label='x')
# axs[7].plot(F_P_e['y']*1e6, 'g', label='y')
# axs[7].plot(F_P_e['z']*1e6, 'b', label='z')
# # # axs[5].fill_between(J_P_e.index, -J_P_e_err*1e9, J_P_e_err*1e9, color='gray', alpha=0.5)

# axs[8].plot(S['x']*1e6, 'r', label='x')
# axs[8].plot(S['y']*1e6, 'g', label='y')
# axs[8].plot(S['z']*1e6, 'b', label='z')
# # axs[6].plot(SdotB, 'k')

# axs[0].set_ylabel(r'B (nT)')
# axs[1].set_ylabel(r'n (cm$^{-3}$)')
# axs[2].set_ylabel(r'PiDi (nW.m$^{-3}$)')
# axs[3].set_ylabel(r'PiDe (nW.m$^{-3}$)')
# axs[4].set_ylabel(r'$T_{e, \perp}/T_{e, \parallel}$')
# axs[5].set_ylabel(r'KE flux ($\mu$W.m$^{-2}$)')
# axs[6].set_ylabel(r'TE flux ($\mu$W.m$^{-2}$)')
# axs[7].set_ylabel(r'PW ($\mu$W.m$^{-2}$)')
# axs[8].set_ylabel(r'S ($\mu$W. m$^{-2}$)')

# # axs[5].set_ylim(-0.4, 0.4)
# # axs[8].set_ylim(-100, 100)

# for ax in axs:
#     ax.legend(loc='upper right', bbox_to_anchor=(1.1, 1))

# plt.xlim(start_time, end_time)

# plt.savefig(rf'transport_{interval}.png', dpi=300)