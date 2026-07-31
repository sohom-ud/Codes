import os
import matplotlib.pyplot as plt
from src.compute_routines.compute_PiD_functions import *
from src.utils.hdf_to_df import hdf_to_df
from datetime import datetime, timedelta

base_dir = rf'/home/sroy/Documents/MMS_events/Shock'

interval = '20171102_042623_20171102_042730'

fname = os.path.join(base_dir, f'{interval}.h5')

df_dict = hdf_to_df(fname)

fs = 15

# plt.rcParams['text.usetex'] = True
# plt.rcParams['font.family'] = 'serif'
plt.rcParams['xtick.labelsize'] = fs
plt.rcParams['ytick.labelsize'] = fs

B = df_dict['b_gse_3']
ni = df_dict['N_ion_3']
ne = df_dict['N_elc_3']
vi = df_dict['v_spincorr_ion_3']
ve = df_dict['v_spincorr_elc_3']
Ti = 1/3. * (df_dict['Temptensor_ion_3']['xx'] + df_dict['Temptensor_ion_3']['yy'] + df_dict['Temptensor_ion_3']['zz'])
Te = 1/3. * (df_dict['Temptensor_elc_3']['xx'] + df_dict['Temptensor_elc_3']['yy'] + df_dict['Temptensor_elc_3']['zz'])

start_time = datetime(2017, 11, 2, 4, 26, 35)
end_time = datetime(2017, 11, 2, 4, 26, 52)

Bmag = np.sqrt(B['x']**2 + B['y']**2 + B['z']**2)

fig, axs = plt.subplots(5, 1, sharex=True, figsize=(8, 8))

axs[0].plot(B.index, B['x'], color='#d62728', lw=2, label='x')
axs[0].plot(B.index, B['y'], color='#2ca02c', lw=2, label='y')
axs[0].plot(B.index, B['z'], color='#1f77b4', lw=2, label='z')
# axs[0].plot(B.index, Bmag, color='k', lw=3, label='mag')

axs[1].plot(ni.index, ni.values, color='#d62728', lw=2, label='ion')
axs[1].plot(ne.index, ne.values, color='#1f77b4', lw=2, label='elc')

axs[2].plot(vi.index, vi['x'], color='#d62728', lw=2, label='x')
axs[2].plot(vi.index, vi['y'], color='#2ca02c', lw=2, label='y')
axs[2].plot(vi.index, vi['z'], color='#1f77b4', lw=2, label='z')

axs[3].plot(ve.index, ve['x'], color='#d62728', lw=2, label='x')
axs[3].plot(ve.index, ve['y'], color='#2ca02c', lw=2, label='y')
axs[3].plot(ve.index, ve['z'], color='#1f77b4', lw=2, label='z')

axs[4].plot(Ti.index, Ti.values, color='#d62728', lw=2, label='ion')
axs[4].plot(Te.index, Te.values, color='#1f77b4', lw=2, label='elc')

upstream_start = datetime(2017, 11, 2, 4, 26, 35)
upstream_end = datetime(2017, 11, 2, 4, 26, 36)

foot_start = datetime(2017, 11, 2, 4, 26, 42, 500000)
foot_end = datetime(2017, 11, 2, 4, 26, 44, 500000)

ramp_start = datetime(2017, 11, 2, 4, 26, 44, 500000)
ramp_end = datetime(2017, 11, 2, 4, 26, 47, 800000)

downstream_start = datetime(2017, 11, 2, 4, 26, 50)
downstream_end = datetime(2017, 11, 2, 4, 26, 52)

axs[0].set_ylabel(r'B[nT]', fontsize=fs)
axs[1].set_ylabel(r'n[cm$^{-3}$]', fontsize=fs)
axs[2].set_ylabel(r'v$_\mathrm{i}$[km/s]', fontsize=fs)
axs[3].set_ylabel(r'v$_\mathrm{e}$[km/s]', fontsize=fs)
axs[4].set_ylabel(r'T[eV]', fontsize=fs)

for ax in axs:
    ax.legend(fontsize=fs, bbox_to_anchor=(1, 1), loc='upper left')

for ax in axs:
    ax.axvspan(upstream_start, upstream_end, color='#5b7fa6', alpha=0.3)
    ax.axvspan(foot_start, foot_end, color='#6a5f8f', alpha=0.3)
    ax.axvspan(ramp_start, ramp_end, color='#a5632f', alpha=0.3)
    ax.axvspan(downstream_start, downstream_end, color='#8a7b6d', alpha=0.3)
    # ax.axhline(0, color='gray')
    # ax.axvline(shock_crossing, color='brown', ls='--')

plt.xlim(start_time, end_time)

x_annotate = -0.12
y_annotate = 1.1

trans = axs[0].get_xaxis_transform()

axs[0].text(upstream_start, 1.05, 'Upstream', fontsize=fs, color='#5b7fa6', transform=trans)
axs[0].text(foot_start + timedelta(seconds=0.5), 1.05, 'Foot', fontsize=fs, color='#6a5f8f', transform=trans)
axs[0].text(ramp_start + timedelta(seconds=0.7), 1.05, 'Ramp', fontsize=fs, color='#a5632f', transform=trans)
axs[0].text(downstream_start - timedelta(seconds=0.5), 1.05, 'Downstream', fontsize=fs, color='#8a7b6d', transform=trans)

axs[0].text(x_annotate, y_annotate, '(a)', fontsize=fs, transform=axs[0].transAxes)
axs[1].text(x_annotate, y_annotate, '(b)', fontsize=fs, transform=axs[1].transAxes)
axs[2].text(x_annotate, y_annotate, '(c)', fontsize=fs, transform=axs[2].transAxes)
axs[3].text(x_annotate, y_annotate, '(d)', fontsize=fs, transform=axs[3].transAxes)
axs[4].text(x_annotate, y_annotate, '(e)', fontsize=fs, transform=axs[4].transAxes)

plt.xlabel('Datetime', fontsize=fs)

plt.suptitle('MMS Observation', fontsize=fs*1.5, y=0.99)

# axs[0].text(-0.2, 1,'[nW/m$^3$]', color='k', fontsize=9, transform=axs[0].transAxes)

# axs[4].set_ylim(-5, 5)

plt.tight_layout()

plt.savefig(f'PiD_shock_20171102_overview.png', bbox_inches='tight', dpi=300)

plt.close()