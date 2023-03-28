import numpy as np
import xarray as xr
import pandas as pd
from pathlib import Path
from datetime import date

import sys
sys.path.append('/homedata/nmpnguyen/IPRAL/NETCDF/')
from flag_functions import filter_profile_file, validated_profile, ipral_remove_cloud_profiles

'''
This function is used to create the simple version of netcdf when adding flags into calibrated data 

Input: 
- function flag_functions 
- calibrated data in /homedata/nmpnguyen/RF/Calibrated/zone-3000-4000/*.nc

Output:
- calibrated flagged data in /homedata/nmpnguyen/IPRAL/NETCDF/v-simple/*.nc

Structures of NetCDF:
variables: 
    'calibrated': calibrated data, total attenuated backscatter, 1/m.sr
    'simulated': simulated data, total attenuated molecular backscatter, 1/m.sr
    'flags': flags of calibrated data by profiles: 0=no validated, 1=validated 
dims:
    'time'
    'range'
    'wavelength'
attritbuts: 
'''

def add_attrs_variables(dataset_input, 
                        variable_name, 
                        long_name, 
                        unit, 
                        resolution, 
                        height_of_calibration):
    dict_attrs = {
        "long_name":long_name,
        "unit":unit,
        "resolution":resolution,
        "height_of_calibration_m":height_of_calibration
    }

    for key, value in dict_attrs.items():
        if value is not None:
            dataset_input[variable_name].attrs[key]= value 

    return dataset_input

def adding_flags(calibfile, params_flags):
    # open calibrated file
    data = xr.open_dataset(calibfile)

    # coords
    time = data['time'].values
    range = data['range'].values
    wavelength = data['wavelength'].values

    # Flags
    print('add flags')
    #--------------parameters of flags
    range_limite_top = params_flags['range_limite_top'] #[26000,28000]
    range_limite_bottom = params_flags['range_limite_bottom'] #[2000,3000]
    alt_max = data.attrs['calibration_height'][1]
    limitez = params_flags['limitez']
    rawpath = params_flags['rawpath']
    #-----------------------------------------------------------------------------------
    print(rawpath)
    dataraw = xr.open_dataset(rawpath)
    Range_BckgrdCorr_Signal_355 = (dataraw['rcs_12']/np.square(dataraw['range']) - dataraw['bckgrd_rcs_12'])*np.square(dataraw['range'])
    Range_BckgrdCorr_Signal_532 = (dataraw['rcs_16']/np.square(dataraw['range']) - dataraw['bckgrd_rcs_16'])*np.square(dataraw['range'])

    mask_crit1 = filter_profile_file(Range_BckgrdCorr_Signal_355, 'rcs_12', range_limite_top, range_limite_bottom)
    
    limitez = (range < limitez) #(dataraw['range']<20000)
    # attrs_calibration_heigth = xr.open_dataset(calibpath).attrs['calibration height']
    mask_crit2 = validated_profile(Range_BckgrdCorr_Signal_355.isel(range=limitez))

    try:
        mask_crit3 = ipral_remove_cloud_profiles(alt_max, rawpath)
    except:
        mask_crit3 = np.zeros((mask_crit1.shape), dtype='bool')
        pass
    
    mask_crit355 = mask_crit1.astype('int')*(2**2) + mask_crit2.astype('int')*(2**1) + mask_crit3.astype('int')*(2**0)
    mask_crit1 = filter_profile_file(Range_BckgrdCorr_Signal_532, 'rcs_16', range_limite_top, range_limite_bottom)
    mask_crit2 = validated_profile(Range_BckgrdCorr_Signal_532.isel(range=limitez))
    mask_crit532 = mask_crit1.astype('int')*(2**2) + mask_crit2.astype('int')*(2**1) + mask_crit3.astype('int')*(2**0)

    print('create dataset')
    flags_dataarray = xr.DataArray(
        data = np.vstack((mask_crit355, mask_crit532)).T,
        dims = ['time', 'wavelength'],
        coords = {
            'time': time, 
            'wavelength': wavelength
        },
    )
    data['flags'] = flags_dataarray
    return data

ouputdir = Path('/homedata/nmpnguyen/IPRAL/NETCDF/v_simple/')
rawdir = Path('/bdd/SIRTA/pub/basesirta/1a/ipral/')
inputdir = Path('/homedata/nmpnguyen/IPRAL/RF/Calibrated/zone-3000-4000/')
year = '2020'

for file in sorted(inputdir.glob(f'ipral_1a_Lz1R15mF30sPbck_v01_{year}*.nc')):
    print(file)
    params_flags = {
        'range_limite_top' : [26000,28000],
        'range_limite_bottom' : [2000,3000],
        'limitez' : 20000,
        'rawpath' : sorted(rawdir.glob(f'{year}/**/**/{file.name}'))[0]
    }
    dataout = adding_flags(file, params_flags)
    dataout.to_netcdf(Path(ouputdir, file.name))


