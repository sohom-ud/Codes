'''
Compute the reciprocal vectors of MMS tetrahedron (required for estimation of gradients)
'''
import numpy as np
import pandas as pd
from .interpolate_r import interpolate_r

__all__ = ['compute_k']

def compute_k(data, res='B'):

    interpolate_r(data, res)

    R1 = pd.DataFrame(data[f'r_gse_1_res_{res}']['Values'][:, :3], index=data[f'r_gse_1_res_{res}']['Epoch'], columns=['x', 'y', 'z'])
    R2 = pd.DataFrame(data[f'r_gse_2_res_{res}']['Values'][:, :3], index=data[f'r_gse_2_res_{res}']['Epoch'], columns=['x', 'y', 'z'])
    R3 = pd.DataFrame(data[f'r_gse_3_res_{res}']['Values'][:, :3], index=data[f'r_gse_3_res_{res}']['Epoch'], columns=['x', 'y', 'z'])
    R4 = pd.DataFrame(data[f'r_gse_4_res_{res}']['Values'][:, :3], index=data[f'r_gse_4_res_{res}']['Epoch'], columns=['x', 'y', 'z'])

    r12 = R2 - R1
    r13 = R3 - R1
    r14 = R4 - R1

    r23 = R3 - R2
    r24 = R4 - R2

    r34 = R4 - R3

    n1 = np.cross(r23.values, r24.values)
    d1 = (-r12.values * n1).sum(axis=1)
    k1 = np.divide(n1.T, d1).T

    n2 = np.cross(r34.values, -r13.values)
    d2 = (-r23.values * n2).sum(axis=1)
    k2 = np.divide(n2.T, d2).T

    n3 = np.cross(-r14.values, -r24.values)
    d3 = (-r34.values * n3).sum(axis=1)
    k3 = np.divide(n3.T, d3).T

    n4 = np.cross(r12.values, r13.values)
    d4 = (r14.values * n4).sum(axis=1)
    k4 = np.divide(n4.T, d4).T

    k1 = pd.DataFrame(k1, columns=['x', 'y', 'z'])
    k2 = pd.DataFrame(k2, columns=['x', 'y', 'z'])
    k3 = pd.DataFrame(k3, columns=['x', 'y', 'z'])
    k4 = pd.DataFrame(k4, columns=['x', 'y', 'z'])

    k1.index = r12.index
    k2.index = r12.index
    k3.index = r12.index
    k4.index = r12.index

    k = dict()

    k[f'k_1_res_{res}'] = dict()
    k[f'k_2_res_{res}'] = dict()
    k[f'k_3_res_{res}'] = dict()
    k[f'k_4_res_{res}'] = dict()

    k[f'k_1_res_{res}']['Epoch'] = k1.index
    k[f'k_1_res_{res}']['Values'] = k1.values
    k[f'k_2_res_{res}']['Epoch'] = k2.index
    k[f'k_2_res_{res}']['Values'] = k2.values
    k[f'k_3_res_{res}']['Epoch'] = k3.index
    k[f'k_3_res_{res}']['Values'] = k3.values
    k[f'k_4_res_{res}']['Epoch'] = k4.index
    k[f'k_4_res_{res}']['Values'] = k4.values

    return k