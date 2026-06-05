import pyspedas
import pandas as pd
from src.utils.hdf_to_df import hdf_to_df

PROBES = [1, 2, 3, 4]
COMPS = ['x', 'y', 'z']

m_p = 1.67e-27 # Mass of proton in kg

def correct_moments_bg(fname):

    n_bg_corr = dict()
    v_bg_corr = dict()
    P_bg_corr = dict()
    T_bg_corr = dict()

    # Read file

    df_dict = hdf_to_df(fname)

    for probe in PROBES:

        n = df_dict[f'N_ion_{probe}']
        v = df_dict[f'v_spincorr_ion_{probe}']
        P = df_dict[f'Ptensor_ion_{probe}']

        n_bg = df_dict[f'N_bg_ion_{probe}']
        P_bg = df_dict[f'Ptensor_bg_ion_{probe}']

        n_bg_corr[probe] = df_dict[f'N_ion_{probe}'] - df_dict[f'N_bg_ion_{probe}']
        
        v_bg_corr[probe] = pd.DataFrame(columns=['x', 'y', 'z'])

        for comp in COMPS:
            v_bg_corr[probe][comp] = n.squeeze() * v[comp] / n_bg_corr[probe].squeeze()

        for comp1 in COMPS:
            for comp2 in COMPS:
                P_bg_corr[probe][f'{comp1}{comp2}'] = (P[f'{comp1}{comp2}'] + m_p * n.squeeze() * v[comp1] * v[comp2]) - \
                                                       P_bg[f'{comp1}{comp2}'] - m_p * n_bg_corr.squeeze() * v_bg_corr[comp1] * v_bg_corr[comp2] 
                
                T_bg_corr[probe][f'{comp1}{comp2}'] = P_bg_corr[f'{comp1}{comp2}'] / n.squeeze()