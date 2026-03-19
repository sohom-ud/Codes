import numpy as np
from pyspedas.projects.mms import fpi

# FPI uses a lower bound E0 to rescale the energy bins
E0 = 100 # Units: eV

times = ['2015-10-16/13:05:00', '2015-10-16/13:07:00']

dist = fpi(trange=times, probe=1, data_rate='brst', level='l2', datatype='dis-dist', time_clip=True, center_measurement=True, notplot=True)

# Loading the first record for energy bins, theta bins, phi bins, and distribution

E_bins = dist['mms1_dis_energy_brst']['y'][0] # Energy bin centers
theta_bins = dist['mms1_dis_theta_brst']['y'] # Theta bin centers (elevation)
phi_bins = dist['mms1_dis_phi_brst']['y'][0] # Phi bin centers (azimuth)theta_bins = dist['mms1_dis_theta_brst']['y'] # Theta bin centers (elevation)

f = dist['mms1_dis_dist_brst']['y'][0]

# Converting the energy bins to dimensionless variable U

U_bins = E_bins/(E_bins+100)

# For ease of trapezoidal integration, f(phi=0) is copied to f(phi=360)

phi_bins = np.append(phi_bins, phi_bins[0] + 360.0)
f = np.concatenate([f, f[0:1, :, :]], axis=0)

# Adding f(theta=0)=0 and f(theta=180)=0 data points

theta_bins = np.concatenate([[0.0], theta_bins, [180.0]])

f0_theta = np.zeros((33, 1, 32))

f = np.concatenate([f0_theta, f, f0_theta], axis=1)

# Adding f(U=0)=0 and f(U=1)=0 data points

U_bins = np.concatenate([[0.0], U_bins, [1.0]])

f0_U = np.zeros((33, 18, 1))

f = np.concatenate([f0_U, f, f0_U], axis=2)

# Converting angles from degrees to radians

phi_bins = np.deg2rad(phi_bins)
theta_bins = np.deg2rad(theta_bins)

# Integrating along azimuth (phi)

I_phi = np.trapezoid(f, x=phi_bins, axis=0)

# Integrating along elevation (theta)

sin_theta = np.sin(theta_bins)

integrand_theta = I_phi * sin_theta[:, np.newaxis]

I_theta = np.trapezoid(integrand_theta, x=theta_bins, axis=0)

# Integrating along energy (U)
energy_weight = np.zeros_like(U_bins)
mask = (U_bins>0) & (U_bins<1)
energy_weight[mask] = np.sqrt(U_bins[mask])/(1.0 - U_bins[mask])**2.5

integrand_U = I_theta * energy_weight
I_total = np.trapezoid(integrand_U, x=U_bins)

E0_erg = 100.0 * 1.602e-12
m_p = 1.6726e-24 # mass of proton in grams

prefactor = np.sqrt(2.0) * E0_erg**1.5/ m_p**1.5
n = prefactor * I_total