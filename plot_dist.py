from pyspedas.mms import fpi

trange = ['2017-06-24/23:58:22', '2017-06-24/23:58:37']

data = fpi(trange, notplot=True, data_rate='brst')