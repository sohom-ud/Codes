import numpy as np
from pyspedas.projects.mms import fpi
from scipy.interpolate import RegularGridInterpolator
import plotly.graph_objects as go

e = 1.6e-19
mp = 1.67e-27

times = ['2016-12-09/09:03:54.2', '2016-12-09/09:03:54.4']

dist = fpi(trange=times, probe=1, data_rate='brst', level='l2', datatype=['des-dist', 'des-moms'],
           time_clip=True, notplot=True)

t = 20

E_centers = dist['mms1_des_energy_brst']['y'][t]
theta_centers = np.deg2rad(dist['mms1_des_theta_brst']['y'])
phi_centers = np.deg2rad(dist['mms1_des_phi_brst']['y'][t])

F_grid = dist['mms1_des_dist_brst']['y'][0]

v_centers = np.sqrt(2 * E_centers * e / mp)

dtheta = np.gradient(theta_centers)
dphi = np.gradient(phi_centers)

v3, theta3, phi3 = np.meshgrid(v_centers, theta_centers, phi_centers, indexing='ij')

F_pad = np.pad(F_grid, ((1, 1), (0, 0), (0, 0)), mode='edge')
F_pad = np.pad(F_pad, ((0, 0), (1, 1), (0, 0)), mode='edge')
F_pad = np.pad(F_pad, ((0, 0), (0, 0), (1, 1)), mode='edge')
F_pad[-1] = 0.0

v_pad = np.concatenate([[v_centers[0]*1e-3], v_centers, [v_centers[-1] * 10.0]])
theta_pad = np.concatenate([[0.0], theta_centers, [np.pi]])
phi_pad = np.concatenate([[0.0], phi_centers, [2*np.pi]])

floor = 1e-40
LOGF_pad = np.log10(np.maximum(F_pad, floor))

interp = RegularGridInterpolator(
    (np.log(v_pad), theta_pad, phi_pad), LOGF_pad,
    method="linear", bounds_error=False, fill_value=np.log10(floor),
)

N_cart = 64
lim = 1000
g = np.linspace(-lim, lim, N_cart)

VX, VY, VZ = np.meshgrid(g, g, g, indexing='ij')

V = np.sqrt(VX**2 + VY**2 + VZ**2)
V = np.maximum(V, v_pad[0])
THETA = np.arccos(np.clip(VZ/V, -1.0, 1.0))
PHI = np.arctan2(VY, VX) % (2.0 * np.pi)

pts = np.stack([np.log(V).ravel(), THETA.ravel(), PHI.ravel()], axis=-1)
LOGF = interp(pts).reshape(V.shape)

vmax = LOGF.max()
DECADES = 4.0

fig = go.Figure(go.Volume(
    x=(VX/1e3).ravel(), y=(VY.ravel()/1e3).ravel(), z=(VZ/1e3).ravel(),
    value=LOGF.ravel(),
    isomin=vmax - DECADES,
    isomax=vmax,
    opacity=0.08,
    surface_count=21,
    colorscale="Inferno",    
))

fig.show()