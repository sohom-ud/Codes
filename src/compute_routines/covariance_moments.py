import numpy as np
from pyspedas.projects.mms import fpi

E0 = 100

times = ['2015-10-16/13:05:00', '2015-10-16/13:07:00']

dist = fpi(trange=times, probe=1, data_rate='brst', level='l2', datatype='dis-dist', time_clip=True, center_measurement=True, notplot=True)

E_bins = dist['mms1_dis_energy_brst']['y'][0]
theta_bins = dist['mms1_dis_theta_brst']['y']
phi_bins = dist['mms1_dis_phi_brst']['y'][0]

f = dist['mms1_dis_dist_brst']['y'][0]
f_err = dist['mms1_dis_disterr_brst']['y'][0]

U_bins = E_bins / (E_bins + 100)

