"""
A Maxwellian defined on an MMS/FPI-like velocity-space grid, volume-rendered
with plotly.

Instrument grid:  32 log-spaced energy bins x 16 polar bins x 32 azimuth bins
                  (= 16384 bins, 512 look directions of 11.25 deg each)

plotly's go.Volume needs a *regular Cartesian* grid, so the workflow is:

    1. build the spherical (E, theta, phi) grid and evaluate f on it
    2. pad it so it covers the whole Cartesian cube
    3. for every Cartesian node, invert to (v, theta, phi) and interpolate
    4. render log10(f), because f spans ~20 orders of magnitude

Requires: numpy, scipy, plotly
"""

import numpy as np
import plotly.graph_objects as go
from scipy.interpolate import RegularGridInterpolator

# --------------------------------------------------------------------------
# constants (SI)
# --------------------------------------------------------------------------
QE = 1.602176634e-19        # C, also J per eV
MP = 1.67262192369e-27      # kg
ME = 9.1093837015e-31       # kg

# --------------------------------------------------------------------------
# plasma parameters  -- edit these
# --------------------------------------------------------------------------
MASS   = MP                 # ME for the electron instrument (DES)
N_DENS = 5.0e6              # m^-3   (5 cm^-3)
T_EV   = 100.0              # eV

# --------------------------------------------------------------------------
# FPI-like grid
# --------------------------------------------------------------------------
N_E, N_TH, N_PH = 32, 16, 32
E_MIN, E_MAX = 10.0, 3.0e4  # eV, approximate FPI range

# bin centres
e_ctr  = np.logspace(np.log10(E_MIN), np.log10(E_MAX), N_E)
v_ctr  = np.sqrt(2.0 * e_ctr * QE / MASS)               # m/s
d_th   = np.pi / N_TH
d_ph   = 2.0 * np.pi / N_PH
th_ctr = (np.arange(N_TH) + 0.5) * d_th                 # 5.625 .. 174.375 deg
ph_ctr = (np.arange(N_PH) + 0.5) * d_ph                 # 5.625 .. 354.375 deg

# --------------------------------------------------------------------------
# the distribution, as a phase-space density  [s^3 m^-6]
#   int f d3v = n
# --------------------------------------------------------------------------
VTH = np.sqrt(2.0 * T_EV * QE / MASS)                   # sqrt(2kT/m)


def f_maxwell(v, theta, phi):
    """Isotropic Maxwellian. theta/phi unused, but kept in the signature so
    anisotropic or drifting forms drop straight in."""
    return N_DENS * (MASS / (2.0 * np.pi * T_EV * QE)) ** 1.5 * np.exp(-(v / VTH) ** 2)


V3, TH3, PH3 = np.meshgrid(v_ctr, th_ctr, ph_ctr, indexing="ij")
F_grid = f_maxwell(V3, TH3, PH3)                        # shape (32, 16, 32)

# --------------------------------------------------------------------------
# pad the spherical grid so the Cartesian cube is fully covered
#   radial  : flat fill below E_MIN (see note), zero above E_MAX
#   polar   : nearest-value extension to theta = 0 and pi
#   azimuth : wrap, so the phi = 0 seam interpolates continuously
# --------------------------------------------------------------------------
F_pad = np.pad(F_grid, ((1, 1), (0, 0), (0, 0)), mode="edge")
F_pad = np.pad(F_pad,  ((0, 0), (1, 1), (0, 0)), mode="edge")
F_pad = np.pad(F_pad,  ((0, 0), (0, 0), (1, 1)), mode="wrap")
F_pad[-1] = 0.0                                          # beyond the top energy

v_pad  = np.concatenate([[v_ctr[0] * 1e-3], v_ctr, [v_ctr[-1] * 10.0]])
th_pad = np.concatenate([[0.0], th_ctr, [np.pi]])
ph_pad = np.concatenate([[ph_ctr[0] - d_ph], ph_ctr, [ph_ctr[-1] + d_ph]])

# interpolate log10(f) against log(v): both the grid and the distribution are
# closer to linear in that space, so this is far more faithful than linear-f
FLOOR = 1e-40
LOGF_pad = np.log10(np.maximum(F_pad, FLOOR))

interp = RegularGridInterpolator(
    (np.log(v_pad), th_pad, ph_pad), LOGF_pad,
    method="linear", bounds_error=False, fill_value=np.log10(FLOOR),
)

# --------------------------------------------------------------------------
# resample onto a regular Cartesian grid
# --------------------------------------------------------------------------
N_CART = 64                                              # 48 if the browser lags
LIM = 3.5 * VTH
g = np.linspace(-LIM, LIM, N_CART)
VX, VY, VZ = np.meshgrid(g, g, g, indexing="ij")

V = np.sqrt(VX**2 + VY**2 + VZ**2)
V = np.maximum(V, v_pad[0])                              # keep log(V) finite
TH = np.arccos(np.clip(VZ / V, -1.0, 1.0))               # clip: guards arccos
PH = np.arctan2(VY, VX) % (2.0 * np.pi)                  # arctan2, not arctan

pts = np.stack([np.log(V).ravel(), TH.ravel(), PH.ravel()], axis=-1)
LOGF = interp(pts).reshape(V.shape)

# sanity check: the volume integral should return the density
dv = g[1] - g[0]
print(f"density recovered / n = {(10.0**LOGF).sum() * dv**3 / N_DENS:.4f}")

# --------------------------------------------------------------------------
# volume render
# --------------------------------------------------------------------------
vmax = LOGF.max()
DECADES = 4.0                                            # dynamic range shown

fig = go.Figure(go.Volume(
    x=(VX / 1e3).ravel(), y=(VY / 1e3).ravel(), z=(VZ / 1e3).ravel(),
    value=LOGF.ravel(),
    isomin=vmax - DECADES,
    isomax=vmax,
    opacity=0.08,
    surface_count=21,
    colorscale="Inferno",
    colorbar=dict(title="log<sub>10</sub> f<br>[s<sup>3</sup> m<sup>-6</sup>]"),
    caps=dict(x_show=False, y_show=False, z_show=False),
))

fig.update_layout(
    title=(f"Maxwellian on a {N_E}x{N_TH}x{N_PH} FPI-like grid  |  "
           f"n = {N_DENS/1e6:.0f} cm<sup>-3</sup>, T = {T_EV:.0f} eV"),
    scene=dict(
        xaxis_title="v<sub>x</sub> [km/s]",
        yaxis_title="v<sub>y</sub> [km/s]",
        zaxis_title="v<sub>z</sub> [km/s]",
        aspectmode="data",                               # never 'cube'
    ),
    width=850, height=750,
)

fig.show()
# fig.write_html("maxwellian_fpi.html")


# --------------------------------------------------------------------------
# optional: the raw instrument sampling, one marker per bin centre.
# Worth a look -- the volume render above is interpolated onto 64^3 nodes and
# so looks much smoother than the 16384 bins it was built from.
# --------------------------------------------------------------------------
def plot_bin_centres(decades=6.0):
    x = (V3 * np.sin(TH3) * np.cos(PH3) / 1e3).ravel()
    y = (V3 * np.sin(TH3) * np.sin(PH3) / 1e3).ravel()
    z = (V3 * np.cos(TH3) / 1e3).ravel()
    c = np.log10(np.maximum(F_grid, FLOOR)).ravel()
    keep = c > c.max() - decades

    go.Figure(go.Scatter3d(
        x=x[keep], y=y[keep], z=z[keep], mode="markers",
        marker=dict(size=2, color=c[keep], colorscale="Inferno", opacity=0.5,
                    colorbar=dict(title="log<sub>10</sub> f")),
    )).update_layout(
        title="FPI bin centres (nested spherical shells)",
        scene=dict(xaxis_title="v<sub>x</sub> [km/s]",
                   yaxis_title="v<sub>y</sub> [km/s]",
                   zaxis_title="v<sub>z</sub> [km/s]", aspectmode="data"),
    ).show()
