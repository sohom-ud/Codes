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

base_dir = r'/home/sroy/Documents/MMS_events/Shock'

# interval = r'20151007_113000_20151007_114000'
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

PiDi = -filtered_data['PiDi']['Values'][...]
PiDe = -filtered_data['PiDe']['Values'][...]

pthi = -filtered_data['pthi']['Values'][...]
pthe = -filtered_data['pthe']['Values'][...]
# PS = PSi + PSe
# Lambda_ub_i = -filtered_data['Lambda_ub_i']['Values']['1'][...]
# Lambda_ub_e = -filtered_data['Lambda_ub_e']['Values']['1'][...]
# Lambda_ub_i_e = Lambda_ub_i + Lambda_ub_e

# start_time = datetime(2015, 10, 7, 11, 34, 50)
# end_time = datetime(2015, 10, 7, 11, 35, 30)

start_time = datetime(2017, 11, 2, 4, 26, 35)
end_time = datetime(2017, 11, 2, 4, 26, 55)

start_time_idx = np.abs(epoch - start_time).argmin()
end_time_idx = np.abs(epoch - end_time).argmin()

tscales = filtered_data['PSi']['tscales'][...]

fs = 150

plt.rcParams['xtick.labelsize'] = fs
plt.rcParams['ytick.labelsize'] = fs
plt.rcParams['font.family'] = 'serif'
plt.rcParams['text.usetex'] = True

fig, axs = plt.subplots(4, 1, figsize=(60, 56), sharex=True)

PiDi_im = axs[0].imshow(PiDi, origin='lower', cmap='RdBu_r', aspect='auto', 
           extent=[epoch[0], epoch[-1], tscales[0]/1000, tscales[-1]/1000],
           vmin=-0.5, vmax=0.5,
           interpolation='gaussian')

PiDe_im = axs[1].imshow(PiDe, origin='lower', cmap='RdBu_r', aspect='auto', 
           extent=[epoch[0], epoch[-1], tscales[0]/1000, tscales[-1]/1000],
           vmin=-0.1, vmax=0.1,
           interpolation='gaussian')

pthi_im = axs[2].imshow(pthi, origin='lower', cmap='RdBu_r', aspect='auto', 
           extent=[epoch[0], epoch[-1], tscales[0]/1000, tscales[-1]/1000],
           vmin=-1, vmax=1,
           interpolation='gaussian')

pthe_im = axs[3].imshow(pthe, origin='lower', cmap='RdBu_r', aspect='auto', 
           extent=[epoch[0], epoch[-1], tscales[0]/1000, tscales[-1]/1000],
           vmin=-1, vmax=1,
           interpolation='gaussian')

locator = mdates.SecondLocator(interval=5)
formatter = mdates.ConciseDateFormatter(locator)
plt.gca().xaxis.set_major_locator(locator)
plt.gca().xaxis.set_major_formatter(formatter)

plt.xlim(start_time, end_time)

#------------- Adding colorbars -------------------------------------

cax_size = '2%'
cax_pad = 0.07
cax_x = 5
cax_y = 0.5

divider = make_axes_locatable(axs[0])
cax = divider.append_axes('right', size=cax_size, pad=cax_pad)
fig.colorbar(PiDi_im, cax=cax, orientation='vertical')
cax.text(cax_x, cax_y, r'$PiD_\mathrm{i}^<$', transform=cax.transAxes, va='center', rotation='horizontal', fontsize=fs)
cax.set_title('[nW/m$^3]$', fontsize=fs, y=1.1)

divider = make_axes_locatable(axs[1])
cax = divider.append_axes('right', size=cax_size, pad=cax_pad)
fig.colorbar(PiDe_im, cax=cax, orientation='vertical')
cax.text(cax_x, cax_y, r'$PiD_\mathrm{e}^<$', transform=cax.transAxes, va='center', rotation='horizontal', fontsize=fs)

divider = make_axes_locatable(axs[2])
cax = divider.append_axes('right', size=cax_size, pad=cax_pad/2)
fig.colorbar(pthi_im, cax=cax, orientation='vertical')
# cax.text(cax_x, cax_y, r'$\sum\limits_\alpha PS_\alpha^<$', transform=cax.transAxes, va='center', rotation='horizontal', fontsize=fs)
cax.text(cax_x, cax_y, r'$p\theta_\mathrm{i}^<$', transform=cax.transAxes, va='center', rotation='horizontal', fontsize=fs)

divider = make_axes_locatable(axs[3])
cax = divider.append_axes('right', size=cax_size, pad=cax_pad)
fig.colorbar(pthe_im, cax=cax, orientation='vertical')
cax.text(cax_x, cax_y, r'$p\theta_\mathrm{e}^<$', transform=cax.transAxes, va='center', rotation='horizontal', fontsize=fs)

shock_crossing_time = datetime(2015, 10, 7, 11, 35, 8)

for ax in axs:
    ax.set_ylabel(r'$\tau$(s)', fontsize=fs)
    ax.set_ylim(bottom=0.06, top=50) #0.3s - 60s
    ax.set_yscale('log')
    ax.plot(di_loc.index, di_loc.values/vi, 'black', ls='--', lw=12)
    ax.plot(rhoi_loc.index, rhoi_loc.values/vi, 'black', lw=12)
    ax.text(start_time + timedelta(milliseconds=20), 1.5, r'$d_\mathrm{i}$', fontsize=fs)
    ax.text(start_time + timedelta(milliseconds=20), 0.1, r'$\rho_\mathrm{i}$', fontsize=fs)
    ax.axvline(shock_crossing_time, color='green', ls='--', lw=12)

axs[-1].set_xlabel('Datetime', fontsize=fs)
plt.subplots_adjust(hspace=0.1)

x_annotate = -0.06
y_annotate = 1

axs[0].text(x_annotate, y_annotate, '(a)', fontsize=fs, transform=axs[0].transAxes)
axs[1].text(x_annotate, y_annotate, '(b)', fontsize=fs, transform=axs[1].transAxes)
axs[2].text(x_annotate, y_annotate, '(c)', fontsize=fs, transform=axs[2].transAxes)
axs[3].text(x_annotate, y_annotate, '(d)', fontsize=fs, transform=axs[3].transAxes)

plt.tight_layout()

plt.savefig(r'/home/sroy/Documents/Codes/20171102_042623_20171102_042730_PiD_ptheta.png', dpi=300)
