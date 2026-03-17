from pyspedas.projects.mms import mms_part_getspec
from pytplot import get_data, tplot_names

times = ['2015-10-16/13:05:00', '2015-10-16/13:07:00']

mms_part_getspec(instrument='fpi', probe='1', data_rate='brst', species='i',
                 trange=times, output=['moments'], units='eflux')

tplot_names()  # verify they exist

density  = get_data('mms1_dis_dist_brst_density')
velocity = get_data('mms1_dis_dist_brst_velocity')