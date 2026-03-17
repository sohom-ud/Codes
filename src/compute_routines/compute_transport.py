import pandas as pd
import numpy as np
from src.utils.hdf_to_df import hdf_to_df
from src.utils.resample import resample
from src.utils.methods import dot, cross
from src.compute_routines.compute_PiD_functions import compute_jcurl

m = {'ion': 1.67e-27, 'elc': 9.1e-31}
q = {'ion': 1.6e-19, 'elc': -1.6e-19}

mu_0 = 4*np.pi*1e-7

PROBES = [1, 2, 3, 4]

#========================================
#             KINETIC ENERGY             
#========================================

def compute_kinetic_energy(fname, species='ion', probe=1, reselectron=True):

    use_reselectron = (species == 'ion') and reselectron
    suffix = f'_{species}_{probe}' + ('_reselectron' if use_reselectron else '')
    n_var = f'N{suffix}'
    v_var = f'v_spincorr{suffix}'

    df_dict = hdf_to_df(fname, vars=[n_var, v_var])

    n = df_dict[n_var] * 1e6 # Density in m^-3
    v = df_dict[v_var] * 1e3 # Velocity in m.s^-1

    n = resample(n, v)

    Ef = 0.5* m[species] * n.squeeze() * (v**2).sum(axis=1) # Units: kg.m^-1.s^-2

    return Ef

def compute_kinetic_energy_flux(fname, species='ion', probe=1, reselectron=True):

    use_reselectron = (species == 'ion') and reselectron
    suffix = f'_{species}_{probe}' + ('_reselectron' if use_reselectron else '')
    v_var = f'v_spincorr{suffix}'

    df_dict = hdf_to_df(fname, vars=[v_var])

    v = df_dict[v_var] * 1e3 # Velocity in m.s^-1

    KE = compute_kinetic_energy(fname, species, probe, reselectron) # Kinetic energy in kg.m^-1.s^-2

    F_K = v.mul(KE, axis=0) # Units: kg.s^-3 or W.m^-2

    return F_K

def compute_kinetic_energy_transport(fname, species='ion', reselectron=True):

    F_K_dict = dict()
    k_dict = dict()

    for probe in PROBES:
            
        use_reselectron = (species == 'ion') and reselectron
        k_suffix = f'_{probe}_res_v' + ('e' if use_reselectron else f'{species[0]}')
        k_var = f'k{k_suffix}'

        df_dict= hdf_to_df(fname, vars=[k_var])

        F_K_dict[probe] = compute_kinetic_energy_flux(fname, species, probe, reselectron) # Units: W.m^-2
        k_dict[probe] = df_dict[k_var] * 1e-3 # Reciprocal vectors in m^-1

    J_K = 0.0

    for probe in [1, 2, 3, 4]:

        J_K += dot(k_dict[probe], F_K_dict[probe]).squeeze()

    # Units: W.m^-3

    return J_K

#========================================
#             THERMAL ENERGY             
#========================================

def compute_thermal_energy(fname, species='ion', probe=1, reselectron=True):

    use_reselectron = (species == 'ion') and reselectron
    suffix = f'_{species}_{probe}' + ('_reselectron' if use_reselectron else '')

    P_var = f'Ptensor{suffix}'

    df_dict = hdf_to_df(fname, vars=[P_var])

    P = df_dict[P_var] * 1e-9 # Pressure tensor in Pa

    E_th = 0.5 * (P['xx'] + P['yy'] + P['zz']) # Units: Pa

    return E_th

def compute_thermal_energy_flux(fname, species='ion', probe=1, reselectron=True):

    use_reselectron = (species == 'ion') and reselectron
    suffix = f'_{species}_{probe}' + ('_reselectron' if use_reselectron else '')
    v_var = f'v_spincorr{suffix}'

    df_dict = hdf_to_df(fname, vars=[v_var])

    v = df_dict[v_var] * 1e3 # Velocity in m.s^-1

    E_th = compute_thermal_energy(fname, species, probe, reselectron) # Thermal energy in Pa

    F_th = v.mul(E_th, axis=0) # Units: Pa.m.s^-1 = kg.s^-3 or W.m^-2

    return F_th

def compute_thermal_energy_transport(fname, species='ion', probe=1, reselectron=True):

    F_th_dict = dict()
    k_dict = dict()

    for probe in PROBES:
            
        use_reselectron = (species == 'ion') and reselectron
        k_suffix = f'_{probe}_res_v' + ('e' if use_reselectron else f'{species[0]}')
        k_var = f'k{k_suffix}'

        df_dict= hdf_to_df(fname, vars=[k_var])

        F_th_dict[probe] = compute_thermal_energy_flux(fname, species, probe, reselectron)
        k_dict[probe] = df_dict[k_var] * 1e-3 # Reciprocal vectors in m^-1

    J_th = 0.0

    for probe in [1, 2, 3, 4]:

        J_th += dot(k_dict[probe], F_th_dict[probe]).squeeze()

        # Units: W.m^-3

    return J_th
#=========================================
#              PRESSURE WORK              
#=========================================
def compute_pressure_work(fname, species='ion', probe=1, reselectron=True):

    use_reselectron = (species == 'ion') and reselectron
    suffix = f'_{species}_{probe}' + ('_reselectron' if use_reselectron else '')
    P_var = f'Ptensor{suffix}'
    v_var = f'v_spincorr{suffix}'

    # Generate variable list
    df_dict = hdf_to_df(fname, vars=[P_var, v_var])

    P = df_dict[P_var] * 1e-9 # Pressure tensor in Pa
    v = df_dict[v_var] * 1e3 # Velocity in m.s^-1

    comps = ['x', 'y', 'z']
    F_P = pd.DataFrame(0.0, columns=['x', 'y', 'z'], index=v.index)

    for i in comps:
        for j in comps:
            F_P[i] += P[f'{i}{j}'] * v[j] # Units: kg.s^-3 or W.m^-2    

    return F_P

def compute_pressure_work_transport(fname, species='ion', probe=1, reselectron=True):

    F_P_dict = dict()
    k_dict = dict()

    for probe in PROBES:
            
        use_reselectron = (species == 'ion') and reselectron
        k_suffix = f'_{probe}_res_v' + ('e' if use_reselectron else f'{species[0]}')
        k_var = f'k{k_suffix}'

        df_dict= hdf_to_df(fname, vars=[k_var])

        F_P_dict[probe] = compute_pressure_work(fname, species, probe, reselectron)
        k_dict[probe] = df_dict[k_var] * 1e-3 # Reciprocal vectors in m^-1

    J_P = pd.Series(0.0, index=k_dict[1].index)

    for probe in PROBES:

        J_P += dot(k_dict[probe], F_P_dict[probe]).squeeze() # Units: W.m^-3

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

        heatflux_dict[probe] = df_dict[heatflux_var] * 1e-3 # Units: W.m^-2
        k_dict[probe] = df_dict[k_var] * 1e-3 # Reciprocal vectors in m^-1

    divq = pd.Series(0.0, index=k_dict[1].index)

    for probe in PROBES:
        divq += dot(k_dict[probe], heatflux_dict[probe]) # Units: W.m^-3

    return divq

def compute_Poynting_flux(fname, probe=1, res='e'): # Resolution can be 'i', 'e' or 'B'

    species = ['ion', 'elc']

    if res in ['i', 'e']:
        for i in range(len(species)):
            if res == species[i][0]:
                s = species[i]     
    else:
        s = 'B'

    E_var = f'edp_dce_gse_{probe}'
    B_var = f'b_gse_{probe}'
    v_var = f'v_spincorr_{s}_{probe}'

    if s!='B':
        df_dict = hdf_to_df(fname, vars=[E_var, B_var, v_var])
    else:
        df_dict = hdf_to_df(fname, vars=[E_var, B_var])

    E = df_dict[E_var] * 1e-3 # Electric field in V.m^-1
    B = df_dict[B_var] * 1e-9 # Magnetic field in T
    B = B.drop('mag', axis=1)

    if s!='B':
        v = df_dict[v_var]
        E = resample(E, v)
        B = resample(B, v)
    else:
        E = resample(E, B)

    S = (1/mu_0) * cross(E, B) # Units: W.m^-2

    return S

def compute_div_Poynting_flux(fname, probe=1, res='e'): # Resolution can be 'i', 'e' or 'B'

    k_dict = dict()
    S_dict = dict()

    for probe in PROBES:

        k_suffix = (f'_{probe}_res_v{res}' if res in ['i', 'e'] else f'_{probe}_res_{res}')
        k_var = f'k{k_suffix}'

        S_dict[probe] = compute_Poynting_flux(fname, probe, res) # Units: W.m^-2

        df_dict= hdf_to_df(fname, vars=[k_var])

        k_dict[probe] = df_dict[k_var] * 1e-3 # Reciprocal vectors in m^-1

    divS = pd.Series(0.0, index=k_dict[1].index)

    for probe in PROBES:
        divS += dot(k_dict[probe], S_dict[probe]).squeeze() # Units: W.m^-3

    return divS