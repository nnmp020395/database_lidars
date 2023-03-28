import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import sys
sys.path.append('/homedata/nmpnguyen/database_lidars/Codes/')
from outils import get_step

class ReshapedMatrice:
    def __init__(self, data, step, axis=0):
        self.axis = axis
        if (axis == 0):
            if len(data.shape) == 1:
                sub_data = [data[n:n+step] for n in range(0, data.shape[axis], step)]
            else:
                sub_data = [data[n:n+step, :] for n in range(0, data.shape[axis], step)]
        else:
            sub_data = [data[:, n:n+step] for n in range(0, data.shape[axis], step)]
        self.data = sub_data
        self.step = step
    
    def count(self):
        counted_data = [((self.data[i] > 0) | np.isnan(self.data[i])).sum(axis=self.axis) for i in range(len(self.data))] 
        final_data = np.vstack(counted_data).T
        return final_data
    
    def density(self):
        density_data = [((self.data[i] > 0) | np.isnan(self.data[i])).sum(axis=self.axis)/self.step for i in range(len(self.data))] 
        final_data = np.vstack(density_data).T
#         final_data = np.where(density_data < 1, np.nan, 1)
        return final_data

    def mean(self):
        mean_data = [np.nanmean(self.data[i], axis=self.axis) for i in range(len(self.data))] 
        final_data = np.vstack(mean_data).T
        return final_data
    
    def close_top(self):
        close_data = [np.nanmax(self.data[i], axis=self.axis) for i in range(len(self.data))] 
        final_data = np.vstack(close_data).T
        return final_data    
    
    def close_bottom(self):
        close_data = [np.nanmin(self.data[i], axis=self.axis) for i in range(len(self.data))] 
        final_data = np.vstack(close_data).T
        return final_data           

def filter_profile_file(data, altitude, limites):
    '''
    Critere 1: flagger si le signal en haut plus qu'en bas
    Input: 
        data: range & background corrected signal
        altitude: canal utilisé
        limiteTop: la zone à haute altitude pour le filtre
        limiteBottom: la zone à basse altitude pour le filtre
    Output:
        index_mask: boolean xarray data, index indiquant des profils validés/invalidés 
    '''
    # 1. MEAN TOP AND BOTTOM SIGNAL
    limite = ((altitude > limites['top'][0]) & (altitude < limites['top'][1]))
    meanTop = np.nanmean(data[:, limite], axis=1)
#     limite = (filecorrected['range']>limiteTop[0]) & (filecorrected['range']<limiteTop[1])
#     meanTop = filecorrected.isel(range=limite).mean(dim='range')
    limite = (altitude > limites['bottom'][0]) & (altitude < limites['bottom'][1])
    meanBottom = np.nanmean(data[:, limite], axis=1)
    # 2. GETTING GOOD PROFILE #selectionner le profil incorrect
    index_mask = (meanTop-meanBottom) > 0 # attention si meantop-meanBottom vient du raw (sélectionner channel) ou filecorrected (pas selectionner channel) 
    return index_mask

class clouds_height:
    def __init__(self, ipral_file):
        self.ipral_file = ipral_file
#         self.clouds_height_array = clouds_height.get_clouds_height(self.ipral_file)
        print(self.ipral_file)
    
    def get_clouds_height(self):
        date = self.ipral_file.name.split('_')[4] #if maindir = /bdd/SIRTA... else self.ipral_file.name.split('_')[3]
        print(date)
        # read CHM15k file
        # ----------------
        CHM15K_PATH = Path("/bdd/SIRTA/pub/basesirta/1a/chm15k/")
        chm15k_file = sorted(CHM15K_PATH.glob(f'**/**/**/chm15k_1a_z1Ppr2R15mF15s_v01_{date}_000000_1440.nc'))
        if not chm15k_file:
            print("No CHM15k file found.")
            print("Quitting.")
            sys.exit(1)

        chm15k_file = chm15k_file[0]
        print(f"CHM15k file found: {chm15k_file}")
        df_cbh = xr.open_dataset(chm15k_file)["cloud_base_height"][:, 0].to_dataframe()#["cloud_base_height"]
        df_cbh = df_cbh.groupby(df_cbh.index.round(freq='15S'), dropna=True).max()['cloud_base_height']
        print(df_cbh.shape)
        # read IPRAL data
        # ----------------
        ipral_data = xr.open_dataset(self.ipral_file)
        ipral_time_array = pd.to_datetime(ipral_data['time'].values).round(freq='15S', nonexistent='NaT')
        
        # get cloud height
        # ---------------    
        cloud_height_array = np.zeros_like(ipral_time_array, dtype='float')    
        for j in range(len(ipral_time_array)):            
            cloud_height_array[j] = df_cbh.iloc[df_cbh.index.get_loc(ipral_time_array[j], method='nearest')]
        
        return cloud_height_array

        
def remove_cloud_profiles(data, cloud_height_array, alt_max, mask = 'top', return_index = False):
    """
    Remove IPRAL profiles containing cloud below a defined altitude.

    Parameters
    ----------
    Input:
        alt_max : float, the maximum altitude of clouds in meters.
        data : matrice input, data need to remove profiles
    Output : boolean xarray data, profiles index that any clouds under the chosen altitude
    """    
    if (mask == 'top'):
        id_profiles_mask = np.where((cloud_height_array < alt_max)|np.isnan(cloud_height_array))[0]
    else:
        id_profiles_mask = np.where((cloud_height_array > alt_max)|np.isnan(cloud_height_array))[0]
    
    masked_data = data[id_profiles_mask, :]
    if (return_index):
        return masked_data, id_profiles_mask
    else:
        return masked_data


def get_all_flags(pathfile, data, altitude, limites):
    '''
    Main function is used to apply all flag ways
    '''
    #-----------------
    mask_crit1 = filter_profile_file(data, altitude, limites)
    mask_crit1 = np.where(~mask_crit1)[0]
#     mask_crit2 = validated_profile(Range_BckgrdCorr_Signal_355.isel(range=limitez))
    try:
        clouds_height_array = clouds_height(pathfile).get_clouds_height()
        mask_crit3_data, mask_crit3 = remove_cloud_profiles(data, clouds_height_array, limites['cloud_top_max'], mask='top', return_index=True)
    except:
        mask_crit3 = np.zeros((data.shape[0]), dtype='bool')
        pass
#     mask_crit355 = mask_crit1.astype('int')*(2**1) + mask_crit2.astype('int')*(2**2) + mask_crit3.astype('int')*(2**3)

    return mask_crit1, mask_crit3


limites = {
    'top': [26000,28000],
    'bottom' : [2000,3000],
    'cloud_top_max' : 3000,
    'limite_z' : None,
    'remove_background' : False,
    'resolution' : 500.0
}
year = 2018
maindir = f'/bdd/SIRTA/pub/basesirta/1a/ipral/{year}'#'/homedata/nmpnguyen/IPRAL/NETCDF/v2/'
pattern = f'**/**/ipral_1a_Lz1R15mF30sPbck_v01_{year}*_000000_1440.nc'#'ipral_calib_03_2018*_000000_1440.nc'


all_counted_points = pd.DataFrame([])
for file in sorted(Path(maindir).glob(pattern)):
    print(file)
    dt = xr.open_dataset(file)

    if (limites['remove_background']):
        data_no_background = dt['RCS_noBckgrd_355'].values
    else:
        data_no_background = (dt['rcs_12']/np.square(dt['range']) - dt['bckgrd_rcs_12'])*np.square(dt['range']).values

    # new range after reshape
    n_step = get_step(dt['range'].values, limites['resolution'])
    new_range = ReshapedMatrice(dt['range'].values, n_step, axis=0).close_top().ravel()

    # get all flags index 
    a0, a = get_all_flags(file, data_no_background, dt['range'].values, limites)

    # density after subsetting by flags
    # density_array = ReshapedMatrice(data_no_background, n_step, axis=1).density() # no setting flags 
    density_array = ReshapedMatrice(data_no_background[np.intersect1d(a0, a), :], n_step, axis=1).density()
    
    # count validated points
    count_points = pd.DataFrame((~np.isnan(np.where(density_array < 0.5, np.nan, density_array))).sum(axis=0), 
             index = new_range, columns=[dt['time'][-1].values])
    all_counted_points = pd.concat([all_counted_points, count_points], axis='columns')

# monthly convert
all_counted_points = all_counted_points.T
results = all_counted_points.groupby(all_counted_points.index.month).sum()

results.to_pickle(f'/homedata/nmpnguyen/database_lidars/ipral_{year}_count_points_monthly.pkl')