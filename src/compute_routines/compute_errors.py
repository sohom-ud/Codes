from src.compute_routines.compute_PiD_functions import *
from src.utils.methods import dot

PROBES = [1, 2, 3, 4]

def compute_theta_err(fname, species='ion', probe=1, reselectron=True):

    v_err_dict = dict()
    k_dict = dict()

    for probe in PROBES:
            
        use_reselectron = (species == 'ion') and reselectron
        suffix = f'_{species}_{probe}' + ('_reselectron' if use_reselectron else '')
        err_suffix = f'_{species}_{probe}'
        k_suffix = f'_{probe}_res_v' + ('e' if use_reselectron else f'{species[0]}')
        P_var = f'Ptensor{suffix}'
        P_err_var = f'Ptensor_err{err_suffix}'
        v_var = f'v_spincorr{suffix}'
        v_err_var = f'v_err{err_suffix}'
        k_var = f'k{k_suffix}'

        df_dict= hdf_to_df(fname, vars=[P_var, v_var, k_var, v_err_var, P_err_var])

        v_err_dict[probe] = df_dict[v_err_var]
        k_dict[probe] = df_dict[k_var]

    theta_err = pd.Series(0, index=k_dict[1].index)

    for probe in PROBES:
        theta_err += dot(k_dict[probe]**2, v_err_dict[probe]**2)
    
    return np.sqrt(theta_err)
