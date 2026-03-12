from src.utils.hdf_to_df import hdf_to_df
from src.utils.methods import dot
from src.utils.resample import resample
import pandas as pd
import numpy as np

m = {'ion': 1.67e-27, 'elc': 9.1e-31}
q = {'ion': 1.6e-19, 'elc': -9.1e-31}

mu_0 = 4*np.pi*1e-7

comps = ['x', 'y', 'z']
PROBES = [1, 2, 3, 4]

def compute_Ef_flux_err(fname, species='ion', probe=1, reselectron=True):

    use_reselectron = (species == 'ion') and reselectron
    suffix = f'_{species}_{probe}' + ('_reselectron' if use_reselectron else '')
    n_var = f'N{suffix}'
    v_var = f'v_spincorr{suffix}'
    n_err_var = f'N_err{suffix}'
    v_err_var = f'v_err{suffix}'

    df_dict = hdf_to_df(fname, vars=[n_var, v_var, n_err_var, v_err_var])

    n = df_dict[n_var].squeeze() * 1e6
    v = df_dict[v_var] * 1e3
    n_err = df_dict[n_err_var].squeeze() * 1e6
    v_err = df_dict[v_err_var] * 1e3

    v2 = (v**2).sum(axis=1)

    dFK_drho = pd.DataFrame(0.0, columns=comps, index=v.index)

    for i in ['x', 'y', 'z']:
        dFK_drho[i] = 0.5 * v2 * v[i]

    dFK_du = pd.DataFrame(0.0, columns=['xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy', 'zz'], index=v.index)

    for i in comps:
        for j in comps:
            dFK_du[f'{i}{j}'] = m[species] * n * v[i] * v[j]
            dFK_du[f'{i}{i}'] += 0.5 * m[species] * n * v2

    err = pd.DataFrame(0.0, columns=comps, index=v.index)

    for i in comps:
        err[i] = np.sqrt(dFK_drho[i]**2 * m[species]**2 * n_err**2 + dFK_du[f'{i}x']**2 * v_err['x']**2 + dFK_du[f'{i}y']**2 * v_err['y']**2 + dFK_du[f'{i}z']**2 * v_err['z']**2)

    # Error units: kg.s^-3 or W.m^-2 

    return err

def compute_Ef_transport_err(fname, species='ion', reselectron=True):

    k_dict = dict()
    F_K_err_dict = dict()

    for probe in PROBES:
            
        use_reselectron = (species == 'ion') and reselectron
        k_suffix = f'_{probe}_res_v' + ('e' if use_reselectron else f'{species[0]}')
        k_var = f'k{k_suffix}'

        df_dict= hdf_to_df(fname, vars=[k_var])

        k_dict[probe] = df_dict[k_var] * 1e-3 #Reciprocal vectors in m^-1

        F_K_err_dict[probe] = compute_Ef_flux_err(fname, species, probe, reselectron) # Flux errors in W.m^-2

    err = pd.Series(0.0, index=k_dict[1].index)

    for i in comps:
        for j in comps:
            for probe in PROBES:
                err += dot(k_dict[probe]**2 , F_K_err_dict[probe]**2).squeeze()

    err = np.sqrt(err)

    return err

def compute_Eth_flux_err(fname, species='ion', probe=1, reselectron=True):

    use_reselectron = (species == 'ion') and reselectron
    suffix = f'_{species}_{probe}' + ('_reselectron' if use_reselectron else '')
    P_var = f'Ptensor{suffix}'
    v_var = f'v_spincorr{suffix}'
    v_err_var = f'v_err{suffix}'
    P_err_var = f'Ptensor_err{suffix}'

    df_dict = hdf_to_df(fname, vars=[v_var, P_var, v_err_var, P_err_var])

    v = df_dict[v_var] * 1e3
    v_err = df_dict[v_err_var] * 1e3

    P = df_dict[P_var] * 1e-9 #Pressure tensor in Pa
    P_err = df_dict[P_err_var] * 1e-9 #Pressure tensor error in Pa

    dFT_dP = (1/4.) * v[i]**2
    dFT_du = (1/4.) * (P['xx'] + P['yy'] + P['zz'])

    err = pd.DataFrame(0.0, columns=comps, index=v.index)

    for i in comps:
        err[i] = np.sqrt(dFT_dP * (P_err['xx']**2 + P_err['yy']**2 + P_err['zz']**2) + dFT_du * v_err[i]**2)

    return err

def compute_Eth_transport_err(fname, species='ion', reselectron=True):

    k_dict = dict()
    F_T_err_dict = dict()

    for probe in PROBES:
            
        use_reselectron = (species == 'ion') and reselectron
        k_suffix = f'_{probe}_res_v' + ('e' if use_reselectron else f'{species[0]}')
        k_var = f'k{k_suffix}'

        df_dict= hdf_to_df(fname, vars=[k_var])

        k_dict[probe] = df_dict[k_var] * 1e-3 #Reciprocal vectors in m^-1

        F_T_err_dict[probe] = compute_Eth_flux_err(fname, species, probe, reselectron) # Flux errors in W.m^-2

    err = pd.Series(0.0, index=k_dict[1].index)

    for i in comps:
        for j in comps:
            for probe in PROBES:
                err += dot(k_dict[probe]**2 , F_T_err_dict[probe]**2).squeeze()

    err = np.sqrt(err)

    return err

def compute_pressure_work_err(fname, species='ion', probe=1, reselectron=True):

    use_reselectron = (species == 'ion') and reselectron
    suffix = f'_{species}_{probe}' + ('_reselectron' if use_reselectron else '')
    P_var = f'Ptensor{suffix}'
    v_var = f'v_spincorr{suffix}'
    v_err_var = f'v_err{suffix}'
    P_err_var = f'Ptensor_err{suffix}'

    df_dict = hdf_to_df(fname, vars=[v_var, P_var, v_err_var, P_err_var])

    v = df_dict[v_var] * 1e3
    v_err = df_dict[v_err_var] * 1e3

    P = df_dict[P_var] * 1e-9 #Pressure tensor in Pa
    P_err = df_dict[P_err_var] * 1e-9 #Pressure tensor error in Pa
    
    err = pd.DataFrame(0.0, columns=comps, index=v.index)

    for i in comps:

        err[i] = v['x']**2 * P_err[f'{i}x']**2 + v['y']**2 * P_err[f'{i}y']**2 + v['z']**2 * P_err[f'{i}z']**2 + \
                 P[f'{i}x']**2 * v_err['x']**2 + P[f'{i}y']**2 * v_err['y']**2 + P[f'{i}z']**2 * v_err['z']**2
        
        err[i] = np.sqrt(err[i])

    return err

def compute_pressure_work_transport_err(fname, species='ion', reselectron=True):

    k_dict = dict()
    P_W_err_dict = dict()

    for probe in PROBES:
            
        use_reselectron = (species == 'ion') and reselectron
        k_suffix = f'_{probe}_res_v' + ('e' if use_reselectron else f'{species[0]}')
        k_var = f'k{k_suffix}'

        df_dict= hdf_to_df(fname, vars=[k_var])

        k_dict[probe] = df_dict[k_var] * 1e-3 #Reciprocal vectors in m^-1

        P_W_err_dict[probe] = compute_pressure_work_err(fname, species, probe, reselectron) # Flux errors in W.m^-2

    err = pd.Series(0.0, index=k_dict[1].index)

    for i in comps:
        for j in comps:
            for probe in PROBES:
                err += dot(k_dict[probe]**2 , P_W_err_dict[probe]**2).squeeze()

    err = np.sqrt(err)

    return err

def compute_heatflux_transport_err(fname, species='ion', reselectron=True):

    k_dict = dict()
    heatflux_err_dict = dict()

    for probe in PROBES:
            
        use_reselectron = (species == 'ion') and reselectron
        k_suffix = f'_{probe}_res_v' + ('e' if use_reselectron else f'{species[0]}')
        suffix = f'_{species}_{probe}' + ('_reselectron' if use_reselectron else '')
        k_var = f'k{k_suffix}'
        heatflux_err_var = f'heatflux_err{suffix}'

        df_dict= hdf_to_df(fname, vars=[k_var, heatflux_err_var])

        k_dict[probe] = df_dict[k_var] * 1e-3 #Reciprocal vectors in m^-1

        heatflux_err_dict[probe] = df_dict[heatflux_err_var]

    err = pd.Series(0.0, index=k_dict[1].index)

    for i in comps:
        for j in comps:
            for probe in PROBES:
                err += dot(k_dict[probe]**2 , heatflux_err_dict[probe]**2).squeeze()

    err = np.sqrt(err)

    return err

def compute_Poynting_flux_err(fname, probe=1, reselectron=True):

    E_var = f'edp_dce_gse_{probe}'
    B_var = f'b_gse_{probe}'

    df_dict = hdf_to_df(fname, vars=[E_var, B_var])

    E = df_dict[E_var] * 1e-3 # Electric field in V/m
    B = df_dict[B_var] * 1e-9 # Magnetic field in T

    B = B.drop('mag', axis=1)

    E = resample(E, B)

    B2 = (B**2).sum(axis=1)
    E2 = (E**2).sum(axis=1)

    sigma_B = 0.1 * 1e-9 # Error in magnetic field components in T
    sigma_E = 0.5 * 1e-3 # Error in electric field components in V/m

    S_err = pd.DataFrame(0.0, columns=comps, index=B.index)

    err1 = E2 * sigma_B**2 + B2 * sigma_E**2

    for i in comps:
        S_err[i] = (1/mu_0) * np.sqrt(err1 - E[i]**2 * sigma_B**2 - B[i]**2 * sigma_E**2)

    return S_err