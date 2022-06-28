import numpy as np
import xarray as xr
import pandas as pd
from pathlib import Path
from datetime import date

import sys
sys.path.append('/homedata/nmpnguyen/IPRAL/NETCDF/')
from flag_functions import filter_profile_file, validated_profile, ipral_remove_cloud_profiles

def create_dataset(rawpath, modelpath, calibpath):
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
    print('add coords')
    # coords
    time = xr.open_dataset(rawpath)['time'].values
    range = xr.open_dataset(rawpath)['range'].values
    altitude = xr.open_dataset(rawpath)['altitude'].values + range
    print('add raw signal')
    # Raw
    dataraw = xr.open_dataset(rawpath)
    Raw_AnalogNR_None_355 = xr.open_dataset(rawpath)['rcs_12'].values
    Raw_AnalogNR_None_532 = xr.open_dataset(rawpath)['rcs_16'].values
    Range_BckgrdCorr_Signal_355 = (dataraw['rcs_12']/np.square(dataraw['range']) - dataraw['bckgrd_rcs_12'])*np.square(dataraw['range'])
    Range_BckgrdCorr_Signal_532 = (dataraw['rcs_16']/np.square(dataraw['range']) - dataraw['bckgrd_rcs_16'])*np.square(dataraw['range'])
    latitude = xr.open_dataset(rawpath)['lat'].values
    latitude = np.round(latitude, 3)
    longitude = xr.open_dataset(rawpath)['lon'].values
    longitude = np.round(longitude, 3)
    print('add thermo model')
    # Thermo
    # Pressure = pd.read_pickle(modelpath)['pression'].unstack(level=1).values
    # Temperature = pd.read_pickle(modelpath)['ta'].unstack(level=1).values
    Pressure = xr.open_dataset(modelpath)['pression'].values
    Temperature = xr.open_dataset(modelpath)['ta'].values
    # Model 
    Model_Molecular_Backscatter_532 = pd.read_pickle(modelpath)['beta532mol'].unstack(level=1).values
    Model_Molecular_Backscatter_355 = pd.read_pickle(modelpath)['beta355mol'].unstack(level=1).values
    print('add calibrated signal')
    # Calibration 
    Total_Attenuated_Molecular_Backscatter_532 = xr.open_dataset(calibpath)['simulated'].sel(wavelength = 532).values
    Total_Attenuated_Molecular_Backscatter_355 = xr.open_dataset(calibpath)['simulated'].sel(wavelength = 355).values
    Calibrated_Attenuated_Backscatter_532 = xr.open_dataset(calibpath)['calibrated'].sel(wavelength = 532).values
    Calibrated_Attenuated_Backscatter_355 = xr.open_dataset(calibpath)['calibrated'].sel(wavelength = 355).values
    Total_ScaterringRatio_532 = Calibrated_Attenuated_Backscatter_532 / Total_Attenuated_Molecular_Backscatter_532
    Total_ScaterringRatio_355 = Calibrated_Attenuated_Backscatter_355 / Total_Attenuated_Molecular_Backscatter_355    
    # Flags
    print('add flags')
    range_limite_top = [26000,28000]
    range_limite_bottom = [2000,3000]
    mask_crit1 = filter_profile_file(Range_BckgrdCorr_Signal_355, 'rcs_12', range_limite_top, range_limite_bottom)
    
    limitez = (dataraw['range']<20000)
    attrs_calibration_heigth = xr.open_dataset(calibpath).attrs['calibration height']
    mask_crit2 = validated_profile(Range_BckgrdCorr_Signal_355.isel(range=limitez))

    try:
        mask_crit3 = ipral_remove_cloud_profiles(attrs_calibration_heigth[1], rawpath)
    except:
        mask_crit3 = np.zeros((mask_crit1.shape), dtype='bool')
        pass
    
    mask_crit355 = mask_crit1.astype('int')*(2**2) + mask_crit2.astype('int')*(2**1) + mask_crit3.astype('int')*(2**0)
    mask_crit1 = filter_profile_file(Range_BckgrdCorr_Signal_532, 'rcs_16', range_limite_top, range_limite_bottom)
    mask_crit2 = validated_profile(Range_BckgrdCorr_Signal_532.isel(range=limitez))
    mask_crit532 = mask_crit1.astype('int')*(2**2) + mask_crit2.astype('int')*(2**1) + mask_crit3.astype('int')*(2**0)
    print('create dataset')
    # dataset
    version = "03"
    dataset_to_netcdf = xr.Dataset({
        'Altitude' : (('range'), altitude),
        'RCS_355' : (('time', 'range'), Raw_AnalogNR_None_355),
        'RCS_532' : (('time', 'range'), Raw_AnalogNR_None_532), 
        'RCS_noBckgrd_355' : (('time', 'range'), Range_BckgrdCorr_Signal_355),
        'RCS_noBckgrd_532' : (('time', 'range'), Range_BckgrdCorr_Signal_532),
        'Pressure' : (('time', 'range'), Pressure),
        'Temperature' : (('time', 'range'), Temperature),
        'Model_Molecular_Backscatter_355' : (('time', 'range'),Model_Molecular_Backscatter_355),
        'Model_Molecular_Backscatter_532' : (('time', 'range'),Model_Molecular_Backscatter_532),
        'Attn_Molecular_Backscatter_355' : (('time', 'range'),Total_Attenuated_Molecular_Backscatter_355),
        'Attn_Molecular_Backscatter_532' : (('time', 'range'),Total_Attenuated_Molecular_Backscatter_532),
        'Total_Calib_Attn_Backscatter_355' : (('time', 'range'),Calibrated_Attenuated_Backscatter_355),
        'Total_Calib_Attn_Backscatter_532' : (('time', 'range'),Calibrated_Attenuated_Backscatter_532), 
        'Total_ScattRatio_355' : (('time', 'range'),Total_ScaterringRatio_355),
        'Total_ScattRatio_532' : (('time', 'range'),Total_ScaterringRatio_532),
        # 'C_355' : (('time'), [None]*len(time)),
        # 'C_532' : (('time'), [None]*len(time)),
        'flags_355' : (('time'), mask_crit355),
        'flags_532' : (('time'), mask_crit532),
        },
        coords = {
            'time' : time,
            'range' : range,
        },
        attrs = {
            'Title' : 'SIRTA IPRAL multiwavelength LIDAR L1A data. Range corrected signal with calibration products',
            'Instrument_name' : 'IPSL HiPerformance multi-wavelength Raman Lidar',
            'Station_name' : 'SIRTA',
            'Start_Datetime' : f'{time[0]}',
            'End_Datetime' : f'{time[-1]}',
            'Latitude' : f'{np.round(latitude,3)}',
            'Latitude_unit' : 'degrees_north',
            'Longitude' : f'{np.round(longitude,3)}',
            'Longitude_unit' : 'degrees_east',
            'Responsable_data' : 'N.M.Phuong Nguyen (phuong.nguyen@latmos.ipsl.fr)',
            'Responsable_instrument' : 'Christophe Pietras (christophe.pietras@lmd.ipsl.fr)',
            'Date_data_created' : date.today().strftime('%d/%m/%Y'),
            'Version': version,
        },
    )
    print('add attrs to variables')

    # add attrs
    
    dataset_to_netcdf = add_attrs_variables(dataset_to_netcdf, 'Altitude', 
                                           'Altitude of measurements', 'm', altitude[1]-altitude[0], f'{attrs_calibration_heigth}')
    dataset_to_netcdf =  add_attrs_variables(dataset_to_netcdf, 'RCS_355',
                                            'Range corrected analog signal at 355nm', 'mV.m^2', None, f'{attrs_calibration_heigth}')
    dataset_to_netcdf =  add_attrs_variables(dataset_to_netcdf, 'RCS_532',
                                            'Range corrected analog signal at 532nm', 'mV.m^2', None, f'{attrs_calibration_heigth}')
    dataset_to_netcdf =  add_attrs_variables(dataset_to_netcdf, 'RCS_noBckgrd_355',
                                             'Range corrected signal with background correction at 355nm', 'mV.m^2', None, f'{attrs_calibration_heigth}')
    dataset_to_netcdf =  add_attrs_variables(dataset_to_netcdf, 'RCS_noBckgrd_532',
                                             'Range corrected signal with background correction at 532nm', 'mV.m^2', None, f'{attrs_calibration_heigth}')
    dataset_to_netcdf =  add_attrs_variables(dataset_to_netcdf, 'Pressure',
                                            'Pressure above sea level from ERA5 hourly', 'Pa', None, None)
    dataset_to_netcdf =  add_attrs_variables(dataset_to_netcdf, 'Temperature',
                                            'Temperature from ERA5 hourly', 'K', None, None )
    dataset_to_netcdf =  add_attrs_variables(dataset_to_netcdf, 'Model_Molecular_Backscatter_355',
                                            'Molecular backscatter model at 355nm calculated with Temperature and Pressure fields', 'm^-1.sr^-1', None, None)
    dataset_to_netcdf =  add_attrs_variables(dataset_to_netcdf, 'Model_Molecular_Backscatter_532',
                                            'Molecular backscatter model at 532nm calculated with Temperature and Pressure fields', 'm^-1.sr^-1', None, None)
    dataset_to_netcdf =  add_attrs_variables(dataset_to_netcdf, 'Attn_Molecular_Backscatter_532',
                                            'Attenuated molecular backscatter at 532nm from molecular backscatter model ', 'm^-1.sr^-1', None, None)
    dataset_to_netcdf =  add_attrs_variables(dataset_to_netcdf, 'Attn_Molecular_Backscatter_355',
                                            'Attenuated molecular backscatter at 355nm from molecular backscatter model ', 'm^-1.sr^-1', None, None)
    dataset_to_netcdf =  add_attrs_variables(dataset_to_netcdf, 'Total_Calib_Attn_Backscatter_355',
                                            'Calibrated signal of Attenuated Backscatter at 355nm', 'm^-1.sr^-1', None, f'{attrs_calibration_heigth}')
    dataset_to_netcdf =  add_attrs_variables(dataset_to_netcdf, 'Total_Calib_Attn_Backscatter_532',
                                            'Calibrated signal of Attenuated Backscatter at 532nm', 'm^-1.sr^-1', None, f'{attrs_calibration_heigth}')
    dataset_to_netcdf =  add_attrs_variables(dataset_to_netcdf, 'Total_ScattRatio_355',
                                            'Total scattering ratio of signal at 355nm = Calibrated signal of Attenuated Backscatter / Attenuated molecular backscatter ', None, None, f'{attrs_calibration_heigth}')
    dataset_to_netcdf =  add_attrs_variables(dataset_to_netcdf, 'Total_ScattRatio_532',
                                            'Total scattering ratio of signal at 532nm = Calibrated signal of Attenuated Backscatter / Attenuated molecular backscatter ', None, None, f'{attrs_calibration_heigth}')
    # dataset_to_netcdf =  add_attrs_variables(dataset_to_netcdf, 'C_532', 
    #                                          'Calibration constant for the 532nm signal', None, None, '3000 - 4000')
    # dataset_to_netcdf =  add_attrs_variables(dataset_to_netcdf, 'C_355', 
    #                                          'Calibration constant for the 355nm signal', None, None, '3000 - 4000') 
    dataset_to_netcdf = add_attrs_variables(dataset_to_netcdf, 'flags_355', 'flags: 0=No flags, 1=profile which has a ratio between the invalidated points and the total points greater than the threshold, \n2=profile which detect many points in height altitude than low altitude., \n4=profile which detects the cloud in the calibration area.', 
                                       'AU', None, None)
    dataset_to_netcdf = add_attrs_variables(dataset_to_netcdf, 'flags_532', 'flags: 0=No flags, 1=profile which has a ratio between the invalidated points and the total points greater than the threshold, \n2=profile which detect many points in height altitude than low altitude., \n4=profile which detects the cloud in the calibration area.', 
                                       'AU', None, None)
    return dataset_to_netcdf, version


year = 2020
RAW_DIR = Path('/bdd/SIRTA/pub/basesirta/1a/ipral/', str(year))
MODEL_DIR = Path('/homedata/nmpnguyen/IPRAL/RF/Simul/')
CALIB_DIR = Path('/homedata/nmpnguyen/IPRAL/RF/Calibrated/zone-3000-4000/')

for rawpath in sorted(RAW_DIR.glob('**/**/ipral_1a_Lz1R15mF30sPbck_v01_202009*_000000_1440.nc')): 
    modelpath = sorted(MODEL_DIR.glob(f'{rawpath.stem}*'))[0]
    modelpath = Path('/home/nmpnguyen/tmp_simul.nc')
    calibpath = sorted(CALIB_DIR.glob(f'{rawpath.name}'))[0]
    print(modelpath, calibpath)
    if (not rawpath) | (not modelpath) | (not calibpath):
        print('Not Found Files')
        pass
    dataset_to_write, version = create_dataset(rawpath, modelpath, calibpath)
    '''
    Before saving new file, checking if the old file are still present, remove it.
    '''
    outputpath = Path('/homedata/nmpnguyen/IPRAL/NETCDF/v2/', f'ipral_calib_{version}_{rawpath.name.split("v01_")[1]}')
    import os 
    if os.path.exists(f'/homedata/nmpnguyen/IPRAL/NETCDF/v2/ipral_calib_{version}_{rawpath.name.split("v01_")[1]}'): #outputpath.isfile()
        print('The file exists --> should be removed')
        os.remove(f'/homedata/nmpnguyen/IPRAL/NETCDF/v2/ipral_calib_{version}_{rawpath.name.split("v01_")[1]}')
        print('Save new file')
    else:
        print('The file does not exist --> can write new file')
    
    print(f'Output file: /homedata/nmpnguyen/IPRAL/NETCDF/v2/ipral_calib_{version}_{rawpath.name.split("v01_")[1]}')
    dataset_to_write.to_netcdf(f'/homedata/nmpnguyen/IPRAL/NETCDF/v2/ipral_calib_{version}_{rawpath.name.split("v01_")[1]}', 'w')


