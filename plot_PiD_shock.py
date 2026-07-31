import os
import matplotlib.pyplot as plt
from matplotlib import rc
from src.compute_routines.compute_PiD_functions import *
from src.utils.hdf_to_df import hdf_to_df
from datetime import datetime, timedelta


base_dir = rf'/home/sroy/Documents/MMS_events/Shock'

interval = '20171102_042623_20171102_042730'
# interval = '20151007_113000_20151007_114000'
# interval = '20151021_074000_20151021_075000'
# interval = '20151025_075500_20151025_081000'

fs = 30

fname = os.path.join(base_dir, f'{interval}.h5')

df_dict = hdf_to_df(fname)

plt.rcParams['xtick.labelsize'] = fs
plt.rcParams['ytick.labelsize'] = fs

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "text.latex.preamble": (
        r"\usepackage[T1]{fontenc}"
        r"\usepackage{sansmathfonts}"
        r"\renewcommand{\familydefault}{\sfdefault}"
        r"\usepackage{bm}"
    ),
})

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

jdotEprime = compute_jdotEprime(fname)

smooth_window = '500ms'

PiDi = PiDi.rolling(window=smooth_window).mean()
PiDe = PiDe.rolling(window=smooth_window).mean()

pthi = pthi.rolling(window=smooth_window).mean()
pthe = pthe.rolling(window=smooth_window).mean()

PSi = PSi.rolling(window=smooth_window).mean()
PSe = PSe.rolling(window=smooth_window).mean()

PStotal = PStotal.rolling(window=smooth_window).mean()
jdotEprime = jdotEprime.rolling(window=smooth_window).mean()

start_time = datetime(2017, 11, 2, 4, 26, 35)
end_time = datetime(2017, 11, 2, 4, 26, 52)

# start_time = datetime(2015, 10, 7, 11, 37)
# end_time = datetime(2015, 10, 7, 11, 37, 35)

# start_time = datetime(2015, 10, 21, 7, 45)
# end_time = datetime(2015, 10, 21, 7, 47)

# start_time = datetime(2015, 10, 25, 8)
# end_time = datetime(2015, 10, 25, 8, 1)

# plt.rcParams['font.family'] = 'serif'

fig, axs = plt.subplots(5, 1, sharex=True, figsize=(15, 15))

axs[0].plot(PiDi.index, -PiDi.values, color='#d62728', label='ion', lw=3)
axs[0].plot(PiDe.index, -PiDe.values, color='#1f77b4', label='elc', lw=3)

axs[1].plot(pthi.index, -pthi.values, color='#d62728', label='ion', lw=3)
axs[1].plot(pthe.index, -pthe.values, color='#1f77b4', label='elc', lw=3)

axs[2].plot(PSi.index, -PSi.values, color='#d62728', label='ion', lw=3)
axs[2].plot(PSe.index, -PSe.values, color='#1f77b4', label='elc', lw=3)

axs[3].plot(PStotal.index, -PStotal.values, color='k', lw=3)

axs[4].plot(jdotEprime.index, jdotEprime.values, color='k', lw=3)

# axs[0].plot(PiDe.index, np.cumsum(-PiDe.values), color='#1f77b4', lw=1, label='electron')
# axs[0].plot(PiDi.index, np.cumsum(-PiDi.values), color='#d62728', label='ion')

# axs[1].plot(pthe.index, np.cumsum(-pthe.values), color='#1f77b4', lw=1, label='electron')
# axs[1].plot(pthi.index, np.cumsum(-pthi.values), color='#d62728', label='ion')

# axs[2].plot(PSe.index, np.cumsum(-PSe.values), color='#1f77b4', lw=1, label='electron')
# axs[2].plot(PSi.index, np.cumsum(-PSi.values), color='#d62728', label='ion')

upstream_start = datetime(2017, 11, 2, 4, 26, 35)
upstream_end = datetime(2017, 11, 2, 4, 26, 36)

foot_start = datetime(2017, 11, 2, 4, 26, 42, 500000)
foot_end = datetime(2017, 11, 2, 4, 26, 44, 500000)

ramp_start = datetime(2017, 11, 2, 4, 26, 44, 500000)
ramp_end = datetime(2017, 11, 2, 4, 26, 47, 800000)

downstream_start = datetime(2017, 11, 2, 4, 26, 50)
downstream_end = datetime(2017, 11, 2, 4, 26, 52)

shock_crossing = datetime(2015, 10, 25, 8, 0, 20)

axs[0].set_ylabel(r'$-\Pi_{ij}^{(\alpha)}\mathrm{D}_{ij}^{(\alpha)}$', fontsize=fs)
axs[1].set_ylabel(r'$-p^{(\alpha)}\theta^{(\alpha)}$', fontsize=fs)
axs[2].set_ylabel(r'$PS_\alpha$', fontsize=fs)
axs[3].set_ylabel(r'$\displaystyle\sum_\alpha PS_\alpha$', fontsize=fs)
axs[4].set_ylabel(r'$\bm{j}\cdot\bm{E}^\prime$', fontsize=fs)

for ax in axs[0:3]:
    ax.legend(fontsize=fs, bbox_to_anchor=(1, 1), loc='upper left')

for ax in axs:
    ax.axvspan(upstream_start, upstream_end, color='#5b7fa6', alpha=0.3)
    ax.axvspan(foot_start, foot_end, color='#6a5f8f', alpha=0.3)
    ax.axvspan(ramp_start, ramp_end, color='#a5632f', alpha=0.3)
    ax.axvspan(downstream_start, downstream_end, color='#8a7b6d', alpha=0.3)
    ax.axhline(0, color='gray')
    # ax.axvline(shock_crossing, color='brown', ls='--')

plt.xlim(start_time, end_time)

plt.xlabel('Datetime (UTC)', fontsize=fs)

trans = axs[0].get_xaxis_transform()

axs[0].text(upstream_start, 1.05, r'Upstream', fontsize=fs, color='#5b7fa6', transform=trans)
axs[0].text(foot_start + timedelta(seconds=0.5), 1.05, 'Foot', fontsize=fs, color='#6a5f8f', transform=trans)
axs[0].text(ramp_start + timedelta(seconds=0.8), 1.05, 'Ramp', fontsize=fs, color='#a5632f', transform=trans)
axs[0].text(downstream_end-timedelta(seconds=2), 1.05, 'Downstream', fontsize=fs, color='#8a7b6d', transform=trans)

axs[0].text(-0.2, 1,'[nW/m$^3$]', color='k', fontsize=fs/1.2, transform=axs[0].transAxes)

# axs[4].set_ylim(-5, 5)

plt.suptitle('MMS Observation', fontsize=fs*1.5, y=0.99)

plt.tight_layout()

plt.savefig(f'PiD_shock_20171102_smooth_{smooth_window}.png', bbox_inches='tight', dpi=300)

plt.close()