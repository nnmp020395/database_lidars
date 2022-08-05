import xarray as xr 
import numpy as np 
import pandas as pd 
import glob, os
import numpy as np
from datetime import datetime, timedelta
from tqdm import tqdm
from pathlib import Path
import scipy.interpolate as spi
import math

"""
  update 2022-03-28
  -------
  Ojectif du script: 
    - obternir les fichiers ERA5
    - convertir time, altitude
    - obtenir pression, temperature -> calculer backscatter coefficient moléculaire
    - convertir au altitude des observations.
    - enregistrer fichier NetCDF temporaire tmp.nc 
  -------

  Input: fichier Opar, Ipral
  Output: fichier temporaire tmp.nc
   
"""

def variables_from_era(file):
    """
    Le script permet de lire input ERA5 et outpt Pression/Temperature
    input = raw file path
    
    """
    print('-----GET TMP FILE-----')
    d = xr.open_dataset(file)
    time = d.time.values
    YEAR = pd.to_datetime(d.attrs['start_time']).strftime('%Y')
    YEARMONTH = pd.to_datetime(d.attrs['start_time']).strftime('%Y%m')
    lon_site = round(4*d.attrs['lon'])/4#round(float(d.geospatial_lon_min),2)
    lat_site = round(4*d.attrs['lat'])/4
    print(f'longitude near era dataset: {lon_site}')
    print(f'latitude near era dataset: {lat_site}')
    #----
    print('-----GET ERA5 FILE-----')
    ERA_FOLDER = Path("/bdd/ERA5/NETCDF/GLOBAL_025/hourly/AN_PL")
    ERA_FILENAME = YEARMONTH+".ap1e5.GLOBAL_025.nc"
    GEOPT_PATH = ERA_FOLDER / YEAR / Path("geopt."+ERA_FILENAME)
    TA_PATH = ERA_FOLDER / YEAR / Path("ta."+ERA_FILENAME)
    print(f'path of temperature {TA_PATH}')
    print(f'path of geopotential {GEOPT_PATH}')
    geopt = xr.open_dataset(GEOPT_PATH)
    ta = xr.open_dataset(TA_PATH)
    #----
    print('-----CONVERT TIME AND LOCALISATION-----')
    # date_start = pd.to_datetime(time[0])
    # date_end = pd.to_datetime(time[-1])
    time = pd.to_datetime(time).strftime('%Y-%m-%dT%H:00:00.000000000').astype('datetime64[ns]')
    time_unique = np.unique(time)
    LAT = geopt.latitude[np.where(np.abs(geopt.latitude.values - lat_site) <=0.25)[0][1]].values
    LON = geopt.longitude[np.where(np.abs(geopt.longitude.values - lon_site) <=0.25)[0][1]].values
    #----
    from timeit import default_timer as timer
    TIME = timer()
    geopt_for_ipral = geopt.sel(time=time_unique, latitude=LAT, longitude=LON).to_dataframe()#['geopt']
    ta_for_ipral = ta.sel(time=time_unique, latitude=LAT, longitude=LON).to_dataframe()#['ta']
    print(f'Time loading {timer()-TIME}')
    #----
    print('-----GETTING PRESSURE AND TEMPERATURE-----')
    lat_site = np.deg2rad(lat_site)
    acc_gravity = 9.78032*(1+5.2885e-3*(np.sin(lat_site))**2 - 5.9e-6*(np.sin(2*lat_site))**2)
    r0 = 2*acc_gravity/(3.085462e-6 + 2.27e-9*np.cos(2*lat_site) - 2e-12*np.cos(4*lat_site))
    g0 = 9.80665
    geopt_for_ipral['geopt_height'] = geopt_for_ipral["geopt"]/g0
    geopt_for_ipral['altitude'] = (geopt_for_ipral['geopt_height']*r0)/(acc_gravity*r0/g0 - geopt_for_ipral['geopt_height'])
    M = 28.966E-3 
    R = 8.314510
    T = (15 + 273.15)
    const = -(M*g0)/(R*T)
    p0 = 101325
    geopt_for_ipral['pression'] = p0*np.exp(const*geopt_for_ipral['altitude'])
    output_era = pd.merge(geopt_for_ipral, ta_for_ipral['ta'], left_index=True, right_index=True) 
    print('variables_from_era --> end')
    return output_era

def simulate_Rayleight_coef(output_era):
    """
    Input is the output dataframe of variables_from_era function --> don't need anymore
    """
    
    print('-----SIMULATE EXTINCTION & BACKSCATTER COEFF. FROM ERA5-----')
    # compute molecular backscatter 
    #----------------------------------------------------------------------------------
    k = 1.38e-23
    const355 = (5.45e-32/1.38e-23)*((355e-3/0.55)**(-4.09))
    const532 = (5.45e-32/1.38e-23)*((532e-3/0.55)**(-4.09))
    output_era['beta355'] = const355*output_era['pression'].div(output_era['ta'])
    output_era['beta532'] = const532*output_era['pression'].div(output_era['ta'])
    # compute extinction coef
    #----------------------------------------------------------------------------------
    ratio = 8*math.pi/3 # 0.119
    output_era['alpha355'] = output_era['beta355']*ratio
    output_era['alpha532'] = output_era['beta532']*ratio
    output_era = output_era.sort_index()
    level = np.unique(output_era.index.get_level_values(0))
    time = np.unique(output_era.index.get_level_values(1)) 
    # compute transmiitance two-way 
    #----------------------------------------------------------------------------------
    # output_era['tau355'] = 0 
    # output_era['tau532'] = 0
    # A = pd.DataFrame()
    # for t in time:
    #     a = output_era.loc[pd.IndexSlice[:,t],:].sort_index(ascending = False)
    #     for i in range(1, a.shape[0]):
    #         a['tau355'].iloc[i] = a['tau355'].iloc[i-1] + a['alpha355'].iloc[i]*(a['altitude'].iloc[i]-a['altitude'].iloc[i-1])
    #         a['tau532'].iloc[i] = a['tau532'].iloc[i-1] + a['alpha532'].iloc[i]*(a['altitude'].iloc[i]-a['altitude'].iloc[i-1])
    #     A = pd.concat((A, a), axis=0)
    # compute attenuated molecular backscatter 
    #----------------------------------------------------------------------------------
    # A['beta355mol'] = A['beta355']*np.exp(-2*A['tau355'])
    # A['beta532mol'] = A['beta532']*np.exp(-2*A['tau532'])
    print('simulate_Rayleight_coef --> end')
    return output_era


def interpolate_atb_mol(file, era):     
    """
    the Input is the output dataframe of simulate_Rayleight_coef function
    """
    print('-----BEFORE INTERPOLATE-----')
    d = xr.open_dataset(file)
    r = d.range.values 
    alt = d.altitude.values
    timeData = d.time.values
    timeEra = np.unique(era.index.get_level_values(1)) 
    time_tmp = np.array(pd.to_datetime(timeData).strftime('%Y-%m-%dT%H:00:00')).astype('datetime64[ns]')
    if len(time_tmp) != len(timeData):
        print("Time Error")
        sys.exit(1)
    #------
    columns_names = ['altitude', 'pression', 'ta', 'beta355', 'beta532', 'alpha355', 'alpha532']#, 'beta355mol', 'beta532mol', 'tau355', 'tau532', 'beta355mol', 'beta532mol'
    pression_interp, ta_interp, beta355_interp ,beta532_interp ,alpha355_interp ,alpha532_interp = [[] for _ in range(len(columns_names)-1)] #, beta355mol_interp ,beta532mol_interp ,tau355_interp ,tau532_interp , beta355mol_interp ,beta532mol_interp, 
    new_index = pd.MultiIndex.from_product([timeData, r], names = ['time', 'range'])
    # df_new = pd.DataFrame(index = new_index, columns = era.columns)
    print('-----INTERPOLATE ATTENUATED BACKSCATTERING FROM ERA5-----')
    for t1 in tqdm(time_tmp):
        a = era.loc[pd.IndexSlice[:, t1], columns_names]
        # f1 = spi.interp1d(a['altitude'], a['beta355mol'], kind = 'linear', bounds_error=False, fill_value="extrapolate")
        # f2 = spi.interp1d(a['altitude'], a['beta532mol'], kind = 'linear', bounds_error=False, fill_value="extrapolate")
        f5 = spi.interp1d(a['altitude'], a['alpha355'], kind = 'linear', bounds_error=False, fill_value="extrapolate")
        f6 = spi.interp1d(a['altitude'], a['alpha532'], kind = 'linear', bounds_error=False, fill_value="extrapolate")
        f7 = spi.interp1d(a['altitude'], a['beta355'], kind = 'linear', bounds_error=False, fill_value="extrapolate")
        f8 = spi.interp1d(a['altitude'], a['beta532'], kind = 'linear', bounds_error=False, fill_value="extrapolate")
        f9 = spi.interp1d(a['altitude'], a['pression'], kind = 'linear', bounds_error=False, fill_value="extrapolate")
        f10 = spi.interp1d(a['altitude'], a['ta'], kind = 'linear', bounds_error=False, fill_value="extrapolate")
        # beta355mol_interp, beta532mol_interp = np.append(beta355mol_interp, np.array(f1(alt))), np.append(beta532mol_interp, np.array(f2(alt)))
        # tau355_interp, tau532_interp = np.append(tau355_interp, np.array(f3(r))), np.append(tau532_interp, np.array(f4(r)))
        alpha355_interp, alpha532_interp = np.append(alpha355_interp, np.array(f5(r))), np.append(alpha532_interp, np.array(f6(r)))
        beta355_interp, beta532_interp = np.append(beta355_interp, np.array(f7(r))), np.append(beta532_interp, np.array(f8(r)))
        pression_interp, ta_interp = np.append(pression_interp, np.array(f9(alt))), np.append(ta_interp, np.array(f10(alt)))
        #print(str(t1))
    #------
    new_df = pd.DataFrame(index = new_index, 
        data = np.array([pression_interp, ta_interp, alpha355_interp, alpha532_interp, beta355_interp, beta532_interp]).T, 
        columns = columns_names[1:])
    # ,tau355_interp ,tau532_interp, beta355mol_interp ,beta532mol_interp
    # lidar = file.parts[4].split('.')[0]
    # new_df.to_pickle(Path("/homedata/nmpnguyen/OPAR/Processed",lidar.upper(),file.name.split(".")[0]+f'_{ratio}'+"_simul.pkl"))
    # print('interpolate_atb_mol --> end')

    #-----SAVE NETCDF-----
    new_df_dataset = xr.Dataset.from_dataframe(new_df)
    return new_df_dataset


def main_simulate(file):
    variables_from_era_data = variables_from_era(file)
    variables_Rayleight_computed = simulate_Rayleight_coef(variables_from_era_data)
    new_netcdf = interpolate_atb_mol(file, variables_Rayleight_computed)
    if os.path.exists('/homedata/nmpnguyen/database_lidars/tmp_simul.nc'):
        os.remove('/home/nmpnguyen/database_lidars/tmp_simul.nc')
    new_netcdf.to_netcdf('/home/nmpnguyen/database_lidars/tmp_simul.nc', 'w')
    return 0

from argparse import Namespace, ArgumentParser
# def main():
parser = ArgumentParser()
parser.add_argument("--path_of_raw_tmpfile", "-file", type=str, help=" ", required=True)
opts = parser.parse_args()
print(opts)

if __name__ == "__main__":
    main_simulate(Path(opts.path_of_raw_tmpfile))