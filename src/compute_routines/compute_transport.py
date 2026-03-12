import pandas as pd
import numpy as np
from src.utils.hdf_to_df import hdf_to_df
from src.utils.resample import resample
from src.utils.methods import dot, cross
from src.compute_routines.compute_PiD_functions import compute_jcurl

m = {'ion': 1.67e-27, 'elc': 9.1e-31}
q = {'ion': 1.6e-19, 'elc': -9.1e-31}

mu_0 = 4*np.pi*1e-7

PROBES = [1, 2, 3, 4]

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

    F_K_dict = dict()
    k_dict = dict()

    for probe in PROBES:
            
        use_reselectron = (species == 'ion') and reselectron
        k_suffix = f'_{probe}_res_v' + ('e' if use_reselectron else f'{species[0]}')
        k_var = f'k{k_suffix}'

        df_dict= hdf_to_df(fname, vars=[k_var])

        F_K_dict[probe] = compute_kinetic_energy_flux(fname, species, reselectron)
        k_dict[probe] = df_dict[k_var]

    J_K = 0.0

    for probe in [1, 2, 3, 4]:

        J_K += dot(k_dict[probe], F_K_dict[probe]).squeeze()

    return J_K


## THERMAL ENERGY
##-------------------------------------------

def compute_thermal_energy(fname, species='ion', probe=1, reselectron=True):

    use_reselectron = (species == 'ion') and reselectron
    suffix = f'_{species}_{probe}' + ('_reselectron' if use_reselectron else '')

    P_var = 'Ptensor{suffix}'

    df_dict = hdf_to_df(fname, vars=[P_var])

    P = df_dict[P_var]

    E_th = 0.5 * (P['xx'] + P['yy'] + P['zz']) * 1e-9 # Pressure tensor is measured in nPa

    return E_th

def compute_thermal_energy_flux(fname, species='ion', probe=1, reselectron=True):

    use_reselectron = (species == 'ion') and reselectron
    suffix = f'_{species}_{probe}' + ('_reselectron' if use_reselectron else '')
    v_var = f'v_spincorr{suffix}'

    df_dict = hdf_to_df(fname, vars=[v_var])

    v = df_dict[v_var]

    E_th = compute_thermal_energy(fname, species, probe, reselectron)

    F_th = (v*1e3).mul(E_th, axis=0)

    return F_th

def compute_thermal_energy_transport(fname, species='ion', probe=1, reselectron=True):

    F_th_dict = dict()
    k_dict = dict()

    for probe in PROBES:
            
        use_reselectron = (species == 'ion') and reselectron
        k_suffix = f'_{probe}_res_v' + ('e' if use_reselectron else f'{species[0]}')
        k_var = f'k{k_suffix}'

        df_dict= hdf_to_df(fname, vars=[k_var])

        F_th_dict[probe] = compute_thermal_energy_flux(fname, species, reselectron)
        k_dict[probe] = df_dict[k_var]

    J_th = 0.0

    for probe in [1, 2, 3, 4]:

        J_th += dot(k_dict[probe], F_th_dict[probe]).squeeze()

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

    F_P_dict = dict()
    k_dict = dict()

    for probe in PROBES:
            
        use_reselectron = (species == 'ion') and reselectron
        k_suffix = f'_{probe}_res_v' + ('e' if use_reselectron else f'{species[0]}')
        k_var = f'k{k_suffix}'

        df_dict= hdf_to_df(fname, vars=[k_var])

        F_P_dict[probe] = compute_pressure_work(fname, species, reselectron)
        k_dict[probe] = df_dict[k_var]

    J_P = pd.Series(0.0, index=k_dict[1].index)

    for probe in PROBES:

        J_P += dot(k_dict[probe], F_P_dict[probe]).squeeze()

    return J_P

def compute_div_q(fname, species='ion', reselectron=True):

    k_dict = dict()
    heatflux_dict = dict()

    for probe in PROBES:

        use_reselectron = (species == 'ion') and reselectron
        suffix = f'_{species}_{probe}' + ('_reselectron' if use_reselectron else '')
        k_suffix = f'_{probe}_res_v' + ('e' if use_reselectron else f'{species[0]}')
        k_var = f'k{k_suffix}'
        heatflux_var = f'heatflux{suffix}'

        df_dict= hdf_to_df(fname, vars=[k_var, heatflux_var])

        heatflux_dict[probe] = df_dict[heatflux_var]
        k_dict[probe] = df_dict[k_var]

    divq = pd.Series(0.0, index=k_dict[1].index)

    for probe in PROBES:
        divq += dot(k_dict[probe], heatflux_dict[probe])

    return divq

def compute_Poynting_flux(fname, probe=1, reselectron=True):

    E_var = f'edp_dce_gse_{probe}'
    B_var = f'b_gse_{probe}'

    df_dict = hdf_to_df(fname, vars=[E_var, B_var])

    E = df_dict[E_var]
    B = df_dict[B_var]
    B = B.drop('mag', axis=1)

    #Convert E and B to SI

    E = E * 1e-3
    B = B * 1e-9

    E = resample(E, B)

    S = (1/mu_0) * cross(E, B) # W.m^-2

    return S