from src.compute_routines.compute_PiD_functions import *
from src.utils.methods import dot

PROBES = [1, 2, 3, 4]

def compute_theta_err(fname, species='ion', probe=1, reselectron=True):

    v_err_dict = dict()
    k_dict = dict()

    for probe in PROBES:
            
        use_reselectron = (species == 'ion') and reselectron
        suffix = f'_{species}_{probe}' + ('_reselectron' if use_reselectron else '')
        k_suffix = f'_{probe}_res_v' + ('e' if use_reselectron else f'{species[0]}')
        v_err_var = f'v_err{suffix}'
        k_var = f'k{k_suffix}'

        df_dict= hdf_to_df(fname, vars=[v_err_var, k_var])

        v_err_dict[probe] = df_dict[v_err_var]
        k_dict[probe] = df_dict[k_var]

    theta_err = pd.Series(0, index=k_dict[1].index)

    for probe in PROBES:
        theta_err += dot(k_dict[probe]**2, v_err_dict[probe]**2)
    
    theta_err = np.sqrt(theta_err)

    return np.sqrt(theta_err)

def compute_P_av_err(fname, species='ion', reselectron=True):

    P_err_dict = dict()

    for probe in PROBES:

        use_reselectron = (species == 'ion') and reselectron
        suffix = f'_{species}_{probe}' + ('_reselectron' if use_reselectron else '')
        P_err_var = f'Ptensor_err{suffix}'

        df_dict = hdf_to_df(fname, vars=[P_err_var])

        P_err_dict[probe] = df_dict[P_err_var]

    P_av_err = pd.DataFrame(0.0, columns=['xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy', 'zz'], index=P_err_dict[1].index)

    for probe in PROBES:
        P_av_err += P_err_dict[probe]**2

    P_av_err = np.sqrt(P_av_err)

    return P_av_err

def compute_p_err(fname, species='ion', reselectron=True):

    P_av_err = compute_P_av_err(fname, species, reselectron)

    p_err = 1/3. * np.sqrt(P_av_err['xx']**2 + P_av_err['yy']**2 + P_av_err['zz']**2)

    return p_err

def compute_ptheta_err(fname, species='ion', reselectron=True):

    p = compute_p(fname, species, reselectron)
    theta = compute_theta(fname, species, reselectron)
    ptheta = compute_ptheta(fname, species, reselectron)

    p_err = compute_p_err(fname, species, reselectron)
    theta_err = compute_theta_err(fname, species, reselectron)

    ptheta_err = ptheta * np.sqrt((theta_err/theta)**2 + (p_err/p)**2)

    return ptheta_err

def compute_gradv_err(fname, species='ion', reselectron=True):

    v_err_dict = dict()
    k_dict = dict()

    for probe in PROBES:
            
        use_reselectron = (species == 'ion') and reselectron
        suffix = f'_{species}_{probe}' + ('_reselectron' if use_reselectron else '')
        k_suffix = f'_{probe}_res_v' + ('e' if use_reselectron else f'{species[0]}')
        v_err_var = f'v_err{suffix}'
        k_var = f'k{k_suffix}'

        df_dict= hdf_to_df(fname, vars=[v_err_var, k_var])

        v_err_dict[probe] = df_dict[v_err_var]
        k_dict[probe] = df_dict[k_var]

    gradv_err = pd.DataFrame(0, columns=['xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy', 'zz'], index=k_dict[1].index)

    for i in ['x', 'y', 'z']:
        for j in ['x', 'y', 'z']:
            for probe in PROBES:
                gradv_err[f'{i}{j}'] += (k_dict[probe][i]**2 * v_err_dict[probe][j]**2)
    
    gradv_err = np.sqrt(gradv_err)

    return gradv_err