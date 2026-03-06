import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from src.utils.hdf_to_df import hdf_to_df
from src.utils.resample import resample
from src.compute_routines.compute_PiD_functions import *

base_dir = r'/home/sroy/Documents/IWF_research/MMS_events/MSH_Turbulence'
interval = r'20180417_183843_20180417_184543.h5'

fname = os.path.join(base_dir, interval)

df_dict = hdf_to_df(fname)

ui = df_dict['v_spincorr_ion_1']
ue = df_dict['v_spincorr_elc_1']
B = df_dict['b_gse_1']

B = B.drop('mag', axis=1)

Bi = resample(B, ui)
Be = resample(B, ue)

sigmac_i = ui['x'] * Bi['x'] + ui['y'] * Bi['y'] + ui['z'] * Bi['z']
sigmac_e = ue['x'] * Be['x'] + ue['y'] * Be['y'] + ue['z'] * Be['z']

uimag = np.sqrt(ui['x']**2 + ui['y']**2 + ui['z']**2)
uemag = np.sqrt(ue['x']**2 + ue['y']**2 + ue['z']**2)

Bimag = np.sqrt(Bi['x']**2 + Bi['y']**2 + Bi['z']**2)
Bemag = np.sqrt(Be['x']**2 + Be['y']**2 + Be['z']**2)

costheta_i = sigmac_i/(uimag * Bimag)
costheta_e = sigmac_e/(uemag * Bemag)

# PSi = compute_PS(fname, species='ion', reselectron=False)
# PSe = compute_PS(fname, species='elc')

# PiDi = compute_PiD(fname, species='ion', reselectron=False)
# PiDe = compute_PiD(fname, species='elc')

pthi = compute_ptheta(fname, species='ion', reselectron=False)
pthe = compute_ptheta(fname, species='elc')

thresholds = np.linspace(0, 1, 1000)

cond_df = pd.DataFrame(columns=['PSi', 'PSe'], index=thresholds)

for i, angle in enumerate(thresholds):
    cond_df['PSi'].iloc[i] = np.mean(np.abs(pthi[np.abs(costheta_i)>angle]))
    cond_df['PSe'].iloc[i] = np.mean(np.abs(pthe[np.abs(costheta_e)>angle]))

fig, axs = plt.subplots(2, 1, sharex=True, figsize=(6, 4))
axs[0].plot(np.rad2deg(np.acos(cond_df.index)), cond_df['PSi'], 'k')
axs[1].plot(np.rad2deg(np.acos(cond_df.index)), cond_df['PSe'], 'k')

plt.xlabel(r'$\theta$ (in degrees)')
# axs[0].set_ylabel(r'$\langle|PiDi|\, |\, |\cos \theta_i |\rangle$')
# axs[1].set_ylabel(r'$\langle|PiDe|\, |\, |\cos \theta_e |\rangle$')
axs[0].set_ylabel(r'$\langle|pthi|\, |\, |\cos \theta_i |\rangle$')
axs[1].set_ylabel(r'$\langle|pthe|\, |\, |\cos \theta_e |\rangle$')
plt.tight_layout()

plt.show()