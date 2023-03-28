import xarray as xr 
import numpy as np 
import pandas as pd 
import glob, os
import numpy as np
from datetime import datetime, timedelta
from tqdm import tqdm
from pathlib import Path



def transfert_opar_database(filepath):
    data = xr.open_dataset(filepath)
    time_arr = data['time'].values
    range_arr = data['range'].values
    latitude_sol = data['signal'].attrs['latitude']
    longitude_sol = data['signal'].attrs['longitude']
    alt_station = data['signal'].attrs['alt_station(km)']
    start_time = data['signal'].attrs['start_time']
    end_time = data['signal'].attrs['end_time']
#     attrs_dict = data.to_dict(data=False)['data_vars']['signal']['attrs']    
    data_dict = {
        "coords": {
            "time": {"dims": "time", "data": time_arr},
            "range": {"dims": "range", "data": range_arr*1e3, "attrs": {"units": "m"}},
            },
        "attrs": {
            "lat": float(latitude_sol),
            "lon": float(longitude_sol),
            "start_time": start_time,
            "end_time": end_time,
             },
        "dims": {"time": len(time_arr), "range": len(range_arr)},
        "data_vars": {
            "altitude": {"dims": "range", "data": (range_arr+alt_station)*1e3, "attrs": {"units": "m"}},
            },
    }
    return data_dict 

def transfert_ipral_database(filepath):
    data = xr.open_dataset(filepath)
    if (pd.to_datetime(data['time'].values[0]).day == pd.to_datetime(Path(filepath).stem.split('_')[4]).day):
        start_time = data['time'].values[0].astype('str')
    else:
        start_time = pd.to_datetime(Path(filepath).stem.split('_')[4]).strftime("%m/%d/%Y, %H:%M:%S")
    print(f'type of Start time: {type(start_time)}')
    end_time = data['time'].values[-1].astype('str')
    print(f'type pf End time: {type(end_time)}')
    time_arr = data['time'].values
    range_arr = data['range'].values  
    altitude = data['altitude'].values 
    latitude_sol = data.attrs['geospatial_lat_min']
    longitude_sol = data.attrs['geospatial_lon_min']
    data_dict = {
        "coords": {
            "time": {"dims": "time", "data": time_arr},
            "range": {"dims": "range", "data": range_arr},
            },
        "attrs": {
            "lat": float(latitude_sol),
            "lon": float(longitude_sol),
            "start_time": start_time,
            "end_time": end_time,
             },
        "dims": {"time": len(time_arr), "range": len(range_arr)},
        "data_vars": {
            "altitude": {"dims": "range", "data": altitude+range_arr},
            "rcs355" : {},
            "rcs532" : {}
            },
    }
    return data_dict

# def transfert_er2_database(filepath):
#     data = xr.open_dataset(filepath, group='')
#     start_time = 

def option_database(type_database, filepath):
    if type_database == 'opar':
        output_dataset_dict = transfert_opar_database(filepath)
        return output_dataset_dict
    elif type_database == 'ipral':
        output_dataset_dict = transfert_ipral_database(filepath)
        return output_dataset_dict
#     elif type_database == 'er2-hsrl2':
#     elif type_database == 'lng_hsrl':
    else:
        print('Please entry type_database')
    


from argparse import Namespace, ArgumentParser
parser = ArgumentParser()
parser.add_argument("--path_of_raw_file", "-file", type=str, help="path file", required=True)
parser.add_argument("--instrument", "-inst", type=str, help="intrusment name", required=True)
opts = parser.parse_args()
print(opts)

if __name__ == "__main__":
    dict_data = option_database(opts.instrument, opts.path_of_raw_file)
    tmp_dataset = xr.Dataset.from_dict(dict_data)
    if os.path.exists('/homedata/nmpnguyen/database_lidars/tmp_file.nc'):
        os.remove('/homedata/nmpnguyen/database_lidars/tmp_file.nc')
    tmp_dataset.to_netcdf('/homedata/nmpnguyen/database_lidars/tmp_file.nc', 'w')

