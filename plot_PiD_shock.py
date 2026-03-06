import os
import matplotlib.pyplot as plt
from src.compute_routines.compute_PiD_functions import *
from src.utils.hdf_to_df import hdf_to_df
from datetime import datetime, timedelta

base_dir = rf'/home/sroy/Documents/MMS_events/Shock'

interval = '20171102_042623_20171102_042730'

fname = os.path.join(base_dir, f'{interval}.h5')

df_dict = hdf_to_df(fname)

plt.rcParams['text.usetex'] = True

B = df_dict['b_gse_1']
ni = df_dict['N_ion_1']
ne = df_dict['N_elc_1']
E = df_dict['edp_dce_gse_1']

PiDi = compute_PiD(fname, species='ion', reselectron=True)
PiDe = compute_PiD(fname, species='elc')

pthi = compute_ptheta(fname, species='ion', reselectron=True)
pthe = compute_ptheta(fname, species='elc')

PSi = compute_PS(fname, species='ion', reselectron=True)
PSe = compute_PS(fname, species='elc')

PStotal = pd.DataFrame(index=PSi.index)
PStotal['PS'] = PSi['PS_ion'] + PSe['PS_elc'] 

jdotE = compute_jdotE(fname)

start_time = datetime(2017, 11, 2, 4, 26, 35)
end_time = datetime(2017, 11, 2, 4, 26, 55)

plt.rcParams['font.family'] = 'serif'

fig, axs = plt.subplots(5, 1, sharex=True, figsize=(6, 8))

axs[0].plot(PiDe.index, -PiDe.values, color='#1f77b4', lw=1, label='ion')
axs[0].plot(PiDi.index, -PiDi.values, color='#d62728', label='electron')

axs[1].plot(pthe.index, -pthe.values, color='#1f77b4', lw=1, label='ion')
axs[1].plot(pthi.index, -pthi.values, color='#d62728', label='electron')

axs[2].plot(PSe.index, -PSe.values, color='#1f77b4', lw=1, label='ion')
axs[2].plot(PSi.index, -PSi.values, color='#d62728', label='electron')

axs[3].plot(PStotal.index, -PStotal.values, color='k', lw=1)

axs[4].plot(jdotE.index, jdotE.values, color='k', lw=1)

upstream_start = datetime(2017, 11, 2, 4, 26, 35)
upstream_end = datetime(2017, 11, 2, 4, 26, 36)

foot_start = datetime(2017, 11, 2, 4, 26, 42, 500000)
foot_end = datetime(2017, 11, 2, 4, 26, 44, 500000)

ramp_start = datetime(2017, 11, 2, 4, 26, 44, 500000)
ramp_end = datetime(2017, 11, 2, 4, 26, 47, 800000)

downstream_start = datetime(2017, 11, 2, 4, 26, 53)
downstream_end = datetime(2017, 11, 2, 4, 26, 55)

axs[0].set_ylabel(r'$-\Pi_{ij}^{(\alpha)}\mathrm{D}_{ij}^{(\alpha)}$', fontsize=12)
axs[1].set_ylabel(r'$-p^{(\alpha)}\theta^{(\alpha)}$', fontsize=12)
axs[2].set_ylabel(r'$-\left(\mathbf{P}^{(\alpha)}.\nabla\right)\cdot\mathbf{u}^{(\alpha)}$', fontsize=12)
axs[3].set_ylabel(r'$-\displaystyle\sum_\alpha \left(\mathbf{P}^{(\alpha)}.\nabla\right)\cdot\mathbf{u}^{(\alpha)}$', fontsize=12)
axs[4].set_ylabel(r'$\mathbf{j}\cdot\mathbf{E}$', fontsize=12)

for ax in axs[0:3]:
    ax.legend(fontsize=10, bbox_to_anchor=(1, 1), loc='upper left')

for ax in axs:
    ax.axvspan(upstream_start, upstream_end, color='#5b7fa6', alpha=0.3)
    ax.axvspan(foot_start, foot_end, color='#6a5f8f', alpha=0.3)
    ax.axvspan(ramp_start, ramp_end, color='#a5632f', alpha=0.3)
    ax.axvspan(downstream_start, downstream_end, color='#8a7b6d', alpha=0.3)
    ax.axhline(0, color='gray')
plt.xlim(start_time, end_time)

axs[0].text(upstream_start, 1, 'Upstream', fontsize=12, color='#5b7fa6')
axs[0].text(foot_start, 1, 'Foot', fontsize=12, color='#6a5f8f')
axs[0].text(ramp_start + timedelta(seconds=0.5), 1, 'Ramp', fontsize=12, color='#a5632f')
axs[0].text(downstream_end - timedelta(seconds=5), 1, 'Downstream', fontsize=12, color='#8a7b6d')

axs[0].text(-0.2, 1,'[nW/m$^3$]', color='k', fontsize=9, transform=axs[0].transAxes)

plt.tight_layout()

plt.savefig('PiD_shock.png', bbox_inches='tight', dpi=300)

plt.close()