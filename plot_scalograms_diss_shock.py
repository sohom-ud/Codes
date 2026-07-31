import os
import h5py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.colors as colors
from datetime import datetime
from datetime import timedelta
from src.utils.hdf_to_df import hdf_to_df
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Rectangle

base_dir = r'/home/sroy/Documents/MMS_events/Shock'

interval = r'20171102_042623_20171102_042730'

fname = os.path.join(base_dir, f'{interval}.h5')
filtered_fname = os.path.join(base_dir, 'scale_filtering', f'{interval}_filtered.h5')

filtered_data = h5py.File(filtered_fname, 'r')

df_dict = hdf_to_df(fname)

ni1 = df_dict['N_ion_1_reselectron']
vi1 = df_dict['v_spincorr_ion_1']
B1 = df_dict['b_gse_1']
Ti1 = df_dict['Temptensor_ion_1_reselectron']
Ti = (Ti1['xx'] + Ti1['yy'] + Ti1['zz']) / 3.

epoch = ni1.index
vi = np.sqrt(np.nanmean(vi1['x']**2 + vi1['y']**2 + vi1['z']**2))
B = np.sqrt(np.nanmean(B1['x']**2 + B1['y']**2 + B1['z']**2))

di_loc = 2.28e7/np.sqrt(ni1) * 1e-5 * 2 * np.pi
rhoi_loc = 1.02e2 * np.sqrt(Ti)/B * 2 * np.pi

PSi = -filtered_data['PSi']['Values'][...]
PSe = -filtered_data['PSe']['Values'][...]
PS = PSi + PSe
Lambda_ub_i = -filtered_data['Lambda_ub_i']['Values']['1'][...]
Lambda_ub_e = -filtered_data['Lambda_ub_e']['Values']['1'][...]
Lambda_ub_i_e = Lambda_ub_i + Lambda_ub_e

start_time = datetime(2017, 11, 2, 4, 26, 35)
end_time = datetime(2017, 11, 2, 4, 26, 52)

start_time_idx = np.abs(epoch - start_time).argmin()
end_time_idx = np.abs(epoch - end_time).argmin()

tscales = filtered_data['PSi']['tscales'][...]

fs = 30

plt.rcParams['xtick.labelsize'] = fs
plt.rcParams['ytick.labelsize'] = fs
# plt.rcParams['font.family'] = 'serif'
# plt.rcParams['font.weight'] = 'bold'
# plt.rcParams['text.usetex'] = True

vmax_PSi = np.percentile(PSi, 99.5)
vmax_PSe = np.percentile(PSe, 99)
vmax_PS = np.percentile(PS, 99)
vmax_Lambda_ub = np.percentile(Lambda_ub_i_e, 99)

fig, axs = plt.subplots(4, 1, figsize=(16, 16), sharex=True)

PSi_im = axs[0].imshow(PSi, origin='lower', cmap='RdBu_r', aspect='auto', 
           extent=[epoch[0], epoch[-1], tscales[0]/1000, tscales[-1]/1000],
           vmin=-1, vmax=1,
           interpolation='gaussian')

PSe_im = axs[1].imshow(PSe, origin='lower', cmap='RdBu_r', aspect='auto', 
           extent=[epoch[0], epoch[-1], tscales[0]/1000, tscales[-1]/1000],
           vmin=-1, vmax=1,
           interpolation='gaussian')

PS_im = axs[2].imshow(PS, origin='lower', cmap='RdBu_r', aspect='auto', 
           extent=[epoch[0], epoch[-1], tscales[0]/1000, tscales[-1]/1000],
           vmin=-1, vmax=1,
           interpolation='gaussian')

# Lambda_ub_i_im = axs[3].imshow(Lambda_ub_i, origin='lower', cmap='RdBu_r', aspect='auto', 
#            extent=[epoch[0], epoch[-1], tscales[0]/1000, tscales[-1]/1000],
#            vmin=-2, vmax=2,
#            interpolation='gaussian')

# Lambda_ub_e_im = axs[4].imshow(Lambda_ub_e, origin='lower', cmap='RdBu_r', aspect='auto', 
#            extent=[epoch[0], epoch[-1], tscales[0]/1000, tscales[-1]/1000],
#            vmin=-30, vmax=30,
#            interpolation='gaussian')

Lambda_ub_i_e_im = axs[3].imshow(Lambda_ub_i_e, origin='lower', cmap='RdBu_r', aspect='auto', 
           extent=[epoch[0], epoch[-1], tscales[0]/1000, tscales[-1]/1000],
           vmin=-3, vmax=3,
           interpolation='gaussian')


locator = mdates.SecondLocator(interval=5)
formatter = mdates.ConciseDateFormatter(locator)
plt.gca().xaxis.set_major_locator(locator)
plt.gca().xaxis.set_major_formatter(formatter)

plt.xlim(start_time, end_time)

#------------- Adding colorbars -------------------------------------

cax_size = '2%'
cax_pad = 0.5
cax_x = 5
cax_y = 0.5

divider = make_axes_locatable(axs[0])
cax = divider.append_axes('right', size=cax_size, pad=cax_pad)
fig.colorbar(PSi_im, cax=cax, orientation='vertical')
cax.text(cax_x, cax_y, r'$PS_\mathrm{i}^<$', transform=cax.transAxes, va='center', rotation='horizontal', fontsize=fs)
# cax.set_title('[nW/m$^3]$', fontsize=fs/1.5, y=1.2)
cax.text(2, 1.1, '[nW/m$^3]$', fontsize=fs//1.4)

divider = make_axes_locatable(axs[1])
cax = divider.append_axes('right', size=cax_size, pad=cax_pad)
fig.colorbar(PSe_im, cax=cax, orientation='vertical')
cax.text(cax_x, cax_y, r'$PS_\mathrm{e}^<$', transform=cax.transAxes, va='center', rotation='horizontal', fontsize=fs)

divider = make_axes_locatable(axs[2])
cax = divider.append_axes('right', size=cax_size, pad=cax_pad)
fig.colorbar(PS_im, cax=cax, orientation='vertical')
# cax.text(cax_x, cax_y, r'$\sum\limits_\alpha PS_\alpha^<$', transform=cax.transAxes, va='center', rotation='horizontal', fontsize=fs)
cax.text(cax_x, cax_y, r'$PS_\mathrm{total}^<$', transform=cax.transAxes, va='center', rotation='horizontal', fontsize=fs)

# divider = make_axes_locatable(axs[3])
# cax = divider.append_axes('right', size=cax_size, pad=cax_pad)
# fig.colorbar(Lambda_ub_i_im, cax=cax, orientation='vertical')
# cax.text(cax_x, cax_y, r'$W_\mathrm{i}^<$', transform=cax.transAxes, va='center', rotation='horizontal', fontsize=fs)

# divider = make_axes_locatable(axs[4])
# cax = divider.append_axes('right', size=cax_size, pad=cax_pad)
# fig.colorbar(Lambda_ub_e_im, cax=cax, orientation='vertical')
# cax.text(cax_x, cax_y, r'$W_\mathrm{e}^<$', transform=cax.transAxes, va='center', rotation='horizontal', fontsize=fs)

divider = make_axes_locatable(axs[3])
cax = divider.append_axes('right', size=cax_size, pad=cax_pad)
fig.colorbar(Lambda_ub_i_e_im, cax=cax, orientation='vertical')
cax.text(cax_x, cax_y, r'$W_\mathrm{total}^<$', transform=cax.transAxes, va='center', rotation='horizontal', fontsize=fs)

upstream_start = datetime(2017, 11, 2, 4, 26, 35)
upstream_end = datetime(2017, 11, 2, 4, 26, 36)

foot_start = datetime(2017, 11, 2, 4, 26, 42, 500000)
foot_end = datetime(2017, 11, 2, 4, 26, 44, 500000)

ramp_start = datetime(2017, 11, 2, 4, 26, 44, 500000)
ramp_end = datetime(2017, 11, 2, 4, 26, 47, 800000)

downstream_start = datetime(2017, 11, 2, 4, 26, 50)
downstream_end = datetime(2017, 11, 2, 4, 26, 52)

axs[0].add_patch(Rectangle((foot_start, 25), (foot_end - foot_start), 15, color='#6a5f8f', clip_on=False))
axs[0].add_patch(Rectangle((ramp_start, 25), (ramp_end - ramp_start), 15, color='#a5632f', clip_on=False))
axs[0].add_patch(Rectangle((downstream_start, 25), (downstream_end - downstream_start), 15, color='#8a7b6d', clip_on=False))
axs[0].add_patch(Rectangle((upstream_start, 25), (upstream_end - upstream_start), 15, color='#5b7fa6', clip_on=False))

shock_crossing_time = datetime(2015, 10, 7, 11, 35, 8)

for ax in axs:
    ax.set_ylabel(r'$\tau$(s)', fontsize=fs)
    ax.set_ylim(bottom=0.06, top=20) #0.3s - 60s
    ax.set_yscale('log')
    ax.plot(di_loc.index, di_loc.values/vi, 'black', ls='--', lw=2)
    ax.plot(rhoi_loc.index, rhoi_loc.values/vi, 'black', lw=2)
    ax.text(start_time + timedelta(milliseconds=20), 3, r'$d_\mathrm{i}$', fontsize=fs)
    ax.text(start_time + timedelta(milliseconds=20), 0.15, r'$\rho_\mathrm{i}$', fontsize=fs)

    ax.axvline(foot_end, color='#6a5f8f', ls='--')
    ax.axvline(foot_start, color='#6a5f8f', ls='--')
    ax.axvline(ramp_end, color='#a5632f', ls='--')
    ax.axvline(ramp_start, color='#a5632f', ls='--')
    ax.axvline(downstream_start, color='#8a7b6d', ls='--')
    ax.axvline(upstream_end, color='#5b7fa6', ls='--')

trans = axs[0].get_xaxis_transform()
axs[0].text(upstream_start, 1.2, 'Upstream', fontsize=fs//1.2, color='#5b7fa6', transform=trans)
axs[0].text(foot_start + timedelta(seconds=0.5), 1.2, 'Foot', fontsize=fs//1.2, color='#6a5f8f', transform=trans)
axs[0].text(ramp_start + timedelta(seconds=0.7), 1.2, 'Ramp', fontsize=fs//1.2, color='#a5632f', transform=trans)
axs[0].text(downstream_start, 1.2, 'Downstream', fontsize=fs//1.2, color='#8a7b6d', transform=trans)

axs[-1].set_xlabel('Datetime', fontsize=fs)
axs[-1].xaxis.set_minor_locator(mdates.MicrosecondLocator(1000000))
plt.subplots_adjust(hspace=0.1)

x_annotate = -0.06
y_annotate = 1

axs[0].text(x_annotate, y_annotate, '(a)', fontsize=fs, transform=axs[0].transAxes)
axs[1].text(x_annotate, y_annotate, '(b)', fontsize=fs, transform=axs[1].transAxes)
axs[2].text(x_annotate, y_annotate, '(c)', fontsize=fs, transform=axs[2].transAxes)
axs[3].text(x_annotate, y_annotate, '(d)', fontsize=fs, transform=axs[3].transAxes)
# axs[4].text(x_annotate, y_annotate, '(e)', fontsize=fs, transform=axs[4].transAxes)
# axs[5].text(x_annotate, y_annotate, '(f)', fontsize=fs, transform=axs[5].transAxes)

plt.suptitle('MMS Observation', fontsize=fs*1.5, y=0.99)

plt.tight_layout()

plt.savefig(r'/home/sroy/Documents/Codes/20171102_042623_20171102_042730_diss.png', dpi=300)
