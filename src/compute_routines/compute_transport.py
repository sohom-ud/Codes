import pandas as pd
import numpy as np
from src.utils.hdf_to_df import hdf_to_df
from src.utils.resample import resample
from src.utils.methods import dot

m = {'ion': 1.67e-27, 'elc': 9.1e-31}
q = {'ion': 1.6e-19, 'elc': -9.1e-31}

def compute_kinetic_energy(fname, species='ion', probe=1, reselectron=True):

    use_reselectron = (species == 'ion') and reselectron
    suffix = f'_{species}_{probe}' + ('_reselectron' if use_reselectron else '')
    n_var = f'N{suffix}'
    v_var = f'v_spincorr{suffix}'

    df_dict = hdf_to_df(fname, vars=[n_var, v_var])

    n = df_dict[n_var]
    v = df_dict[v_var]

    Ef = 0.5* m[species] * n.values * 1e6 * (v**2).sum(axis=1) * 1e6

    return Ef

def compute_kinetic_energy_flux(fname, species='ion', probe=1, reselectron=True):

    use_reselectron = (species == 'ion') and reselectron
    suffix = f'_{species}_{probe}' + ('_reselectron' if use_reselectron else '')
    v_var = f'v_spincorr{suffix}'

    df_dict = hdf_to_df(fname, vars=[v_var])

    v = df_dict[v_var]

    KE = compute_kinetic_energy(fname, species, probe, reselectron)

    F_K = (v*1e3).mul(KE, axis=0)

    return F_K

def compute_kinetic_energy_transport(fname, species='ion', reselectron=True):

    # Generate variable list
    varlist = []
    for probe in [1, 2, 3, 4]:
        if reselectron:
            varlist.extend([f'k_{probe}_res_ve'])
        else:
            varlist.extend([f'k_{probe}_res_v{species[0]}'])

    df_dict = hdf_to_df(fname, vars=varlist)

    k = dict()
    F_K = dict()

    for probe in [1, 2, 3, 4]:

        if species=='ion':
            if reselectron:
                k[probe] = df_dict[f'k_{probe}_res_ve']
                F_K[probe] = compute_kinetic_energy_flux(fname, species, probe, True)
            else:
                k[probe] = df_dict[f'k_{probe}_res_vi']
                F_K[probe] = compute_kinetic_energy_flux(fname, species, probe, False)
        elif species=='elc':
                k[probe] = df_dict[f'k_{probe}_res_ve']
                F_K[probe] = compute_kinetic_energy_flux(fname, species, probe, reselectron)

    J_K = 0.0

    for probe in [1, 2, 3, 4]:

        J_K += dot(k[probe], F_K[probe]).squeeze()

    return J_K


## THERMAL ENERGY
##-------------------------------------------

def compute_thermal_energy(fname, species='ion', probe=1, reselectron=True):

    # Generate variable list
    if species=='ion':
        if reselectron:
            varlist = [f'Ptensor_{species}_{probe}_reselectron']
        else:
            varlist = [f'Ptensor_{species}_{probe}']
    elif species=='elc':
        varlist = [f'Ptensor_{species}_{probe}']

    df_dict = hdf_to_df(fname, vars=varlist)

    P = df_dict[varlist[0]]

    E_th = 0.5 * (P['xx'] + P['yy'] + P['zz']) * 1e-9 # Pressure tensor is measured in nPa

    return E_th

def compute_thermal_energy_flux(fname, species='ion', probe=1, reselectron=True):

    # Generate variable list
    if species=='ion':
        if reselectron:
            varlist = [f'v_spincorr_{species}_{probe}_reselectron']
        else:
            varlist = [f'v_spincorr_{species}_{probe}']
    elif species=='elc':
        varlist = [f'v_spincorr_{species}_{probe}']

    df_dict = hdf_to_df(fname, vars=varlist)

    v = df_dict[varlist[0]]

    E_th = compute_thermal_energy(fname, species, probe, reselectron)

    F_th = (v*1e3).mul(E_th, axis=0)

    return F_th

def compute_thermal_energy_transport(fname, species='ion', probe=1, reselectron=True):

    # Generate variable list
    varlist = []
    for probe in [1, 2, 3, 4]:
        if reselectron:
            varlist.extend([f'k_{probe}_res_ve'])
        else:
            varlist.extend([f'k_{probe}_res_v{species[0]}'])

    df_dict = hdf_to_df(fname, vars=varlist)

    k = dict()
    F_th = dict()

    for probe in [1, 2, 3, 4]:

        if species=='ion':
            if reselectron:
                k[probe] = df_dict[f'k_{probe}_res_ve']
                F_th[probe] = compute_thermal_energy_flux(fname, species, probe, True)
            else:
                k[probe] = df_dict[f'k_{probe}_res_vi']
                F_th[probe] = compute_thermal_energy_flux(fname, species, probe, False)
        elif species=='elc':
                k[probe] = df_dict[f'k_{probe}_res_ve']
                F_th[probe] = compute_thermal_energy_flux(fname, species, probe, reselectron)

    J_th = 0.0

    for probe in [1, 2, 3, 4]:

        J_th += dot(k[probe], F_th[probe]).squeeze()

    return J_th

## Pressure work and transport
##-----------------------------------------------------

def compute_pressure_work(fname, species='ion', probe=1, reselectron=True):

    use_reselectron = (species == 'ion') and reselectron
    suffix = f'_{species}_{probe}' + ('_reselectron' if use_reselectron else '')
    P_var = f'Ptensor{suffix}'
    v_var = f'v_spincorr{suffix}'

    # Generate variable list
    df_dict = hdf_to_df(fname, vars=[P_var, v_var])

    P = df_dict[P_var]
    v = df_dict[v_var]

    comps = ['x', 'y', 'z']
    F_P = pd.DataFrame(0.0, columns=['x', 'y', 'z'], index=v.index)

    for i in comps:
        for j in comps:
            F_P[i] += P[f'{i}{j}'] * v[j]

    return F_P

def compute_pressure_work_transport(fname, species='ion', probe=1, reselectron=True):

