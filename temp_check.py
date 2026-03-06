from pyspedas.mms import fpi
import matplotlib.pyplot as plt
from datetime import datetime

trange = ['2017-07-11/22:33:30', '2017-07-11/22:34:30']

probe = 3

data = fpi(
    trange=trange,
    probe=probe,
    data_rate='brst',
    datatype=f"des-moms",
    center_measurement=True,
    notplot=True
)

Te_perp = data['mms3_des_temppara_brst']

Te_para = data['mms3_des_tempperp_brst']

fig, axs = plt.subplots(1, 1)

plt.plot(Te_perp['x'], Te_perp['y'])
plt.plot(Te_para['x'], Te_para['y'])

start_time = datetime(2017)

plt.show()
