import os
import pandas as pd
import numpy as np

data_dir = r'/home/sroy/Documents/MMS_events/Tail_Reconnection'

parsed_list = os.path.join(data_dir, 'Hubbert2022_traditional_recx_interval_list.csv')
start_times = []
end_times = []
# Parse the Hubbert 2022 interval list
interval_list = pd.read_csv(os.path.join(data_dir, 'Hubbert2022_traditional_recx.csv'), header=None)

for interval in interval_list[1]:

    date, time = interval.split('/')
    month, day, year = date.split('-')
    hh, mm_start = time[:5].split(':')
    mm_end = time[6:]
    ss = '00'

    if len(month) == 1:
        month = f'0{month}'

    start_time = rf'20{year}-{month}-{day}/{hh}:{mm_start}:{ss}'
    end_time = rf'20{year}-{month}-{day}/{hh}:{mm_end}:{ss}'

    start_times.append(start_time)
    end_times.append(end_time)

df = pd.DataFrame(data=np.column_stack([start_times, end_times]))

df.to_csv(parsed_list, index=None, header=None)
