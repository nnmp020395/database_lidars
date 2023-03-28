## Ce script contient spécifiquement des fonctions servant à la détection des flags dans le database Ipral

# Date: _2022.04.12_
# Author: _Phuong Nguyen_ 

import xarray as xr
import numpy as np
from pathlib import Path
import pandas as pd
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import tqdm
import warnings
warnings.filterwarnings("ignore")

'''
tester flag code et verifier avec les quicklooks
enregistrer les quicklooks sur /scratchx/nmpnguyen/IPRAL/raw/
'''

def filter_profile_file(filecorrected, channel, limiteTop, limiteBottom):
    '''
    Critere 1: flagger si le signal en haut plus qu'en bas
    Input: 
        raw: range & background corrected signal
        channel: canal utilisé
        limiteTop: la zone à haute altitude pour le filtre
        limiteBottom: la zone à basse altitude pour le filtre
    Output:
        index_mask: boolean xarray data, index indiquant des profils validés/invalidés 
    '''
    # 1. MEAN TOP AND BOTTOM SIGNAL
    limite = (filecorrected['range']>limiteTop[0]) & (filecorrected['range']<limiteTop[1])
    meanTop = filecorrected.isel(range=limite).mean(dim='range')
    limite = (filecorrected['range']>limiteBottom[0]) & (filecorrected['range']<limiteBottom[1])
    meanBottom = filecorrected.isel(range=limite).mean(dim='range')
    # 2. GETTING GOOD PROFILE #selectionner le profil incorrect
    index_mask = (meanTop-meanBottom) > 0 # attention si meantop-meanBottom vient du raw (sélectionner channel) ou filecorrected (pas selectionner channel) 
    return index_mask


def invalidated_profile(data):
    nb_nonzero = np.count_nonzero(data>0, axis=1)
    nb_points_by_profile = data.shape[1]
    fraction_nonzero = nb_nonzero/nb_points_by_profile
    index_mask = fraction_nonzero > 0.5  #selectionner le profil incorrect    
    return index_mask

def validated_profile(dataipral): 
    '''
    Find the profiles that proportion of validated points compared to total points is higher the given threshold

    Parameters
    ----------
    Input:
        dataipral: range and background corrected dataset
    Output:
        fraction_nonzero: proportion of validated values on total values 
        index_mask: boolean xarray data, mask of validated values following the given threshold
    '''
    nb_nonzero = ((dataipral>0) & ~np.isnan(dataipral)).sum(axis=1)
    nb_points_by_profile = dataipral.shape[1]
    fraction_nonzero = nb_nonzero/nb_points_by_profile
    seuil = np.mean(fraction_nonzero.values)
    index_mask = fraction_nonzero < seuil  #selectionner le profil incorrect    
    return index_mask

import sys
import datetime
def ipral_remove_cloud_profiles(alt_max, ipral_file):
    """
    Remove IPRAL profiles containing cloud below a defined altitude.

    Parameters
    ----------
    Input:
        alt_max : float, the maximum altitude of clouds in meters.
        ipral_file : str or pathlib.Path, the path of the IPRAL file to process.
    Output : boolean xarray data, profiles index that any clouds under the chosen altitude
    """    
    date = ipral_file.name.split('_')[4]
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
    # read IPRAL data
    # ----------------
    ipral_data = xr.open_dataset(ipral_file)
    ipral_time_array = ipral_data['time'].values
    # get cloud height
    # ---------------    
    cloud_height_array = np.zeros_like(ipral_time_array, dtype='float')    
    for j in range(len(ipral_time_array)):
        cloud_height_array[j] = df_cbh.iloc[df_cbh.index.get_loc(ipral_time_array[j], method='nearest')]['cloud_base_height']
    
    cloud_over_altmax = (cloud_height_array > alt_max)|np.isnan(cloud_height_array)
    cloud_under_altmax = (cloud_height_array < alt_max)
    cloud_mask_xarray = xr.DataArray(data=cloud_under_altmax, 
                                     dims=['time'], 
                                     coords=dict(time=ipral_time_array))
    print(f"{len(ipral_time_array)} in IPRAL data")
    profs_to_keep = cloud_over_altmax.astype("i2").sum()
    print(f"{len(ipral_time_array) - profs_to_keep} profiles will be remove")
    return cloud_mask_xarray


def flag_clouds(data, id_Z, rawpath, wave):
    def dispersion_standard_deviation(signal, id_left, seuil, influence):
        '''
        signal : est un vecteur numpy
        id_left : le nombre de valeurs avant la valeur actuelle que nous voulons utiliser pour calculer la ligne de base mobile.
        seuil : le nombre d'écarts types par rapport à la ligne de base mobile qu'un pic doit dépasser pour être compté.
        influence : la quantité d'influence qu'un pic a sur la ligne de base mobile. Celui-ci doit être compris entre 0 et 1.
        '''
        peakIndex = []
        processedSignal = signal[0:id_left]
        for ids in range(id_left, len(signal)):
            y = signal[ids]
            avg = np.nanmean(processedSignal[id_left:ids])
            sd = np.nanstd(processedSignal[id_left:ids])
            if ((y-avg) > (sd*seuil)):
                peakIndex.append(ids)
#                 print(ids, len(processedSignal))
    #             ajustedValued = (influence*y)+((1-influence)*processedSignal[ids-1])
            else:
                processedSignal = np.append(processedSignal, y)
        return peakIndex

    def verification_by_SR(rawpath, wave, id_Z, mask_profiles):
        '''
        rawpath : le chemin du fichier Ipral --> le nom du fichier 
        wave : 532nm et 355nm 
        id_Z : indices des altitudes à étudier les nuages 
        mask_profiles : indices bool des nuages détectés à vérifier
        '''
        # Retrouver le fichier calibré correspondant 
        IPRAL_RF_PATH = Path('/homedata/nmpnguyen/IPRAL/RF/Calibrated')
        IPRAL_RF_FILE = Path(IPRAL_RF_PATH, rawpath.name.split('.')[0]+'.nc')
        print(IPRAL_RF_FILE)
        # Caculer SR, range correspondant aux profiles à vérifier
        datacalib = xr.open_dataset(IPRAL_RF_FILE)
        # datacalib = datacalib.resample(time='15min').mean()
        SR2Darray = (datacalib['calibrated']/datacalib['simulated']).isel(range=id_Z).sel(wavelength = wave).values
        Zlimite2Darray = np.array([datacalib['range'][id_Z].values] * len(datacalib['time']))
        # Retourner des indices indiqués les nuages 
        zcalib_top = 5000 #datacalib.attrs['calibration height'][1]
        selected_indices_profiles = np.where((np.ma.masked_array(SR2Darray, mask=mask_profiles)>1.7) & (np.ma.masked_array(Zlimite2Darray, mask=mask_profiles)<zcalib_top))
        final_indices_profiles = np.unique(selected_indices_profiles[0])
        return final_indices_profiles, zcalib_top
        
    indices_profiles = np.zeros_like(data.isel(range=id_Z).values, dtype=bool)
    for t in range(len(data['time'])):
        indices = dispersion_standard_deviation(data.isel(time=t, range=id_Z).values, 
                                                id_left = 5,
                                                seuil = 4, 
                                                influence = 0.1)
        indices_profiles[t, indices] = True 

    indices_clouds_profiles, zcalib_top = verification_by_SR(rawpath, wave, id_Z, indices_profiles)
    indices_clouds_profiles = np.in1d(np.arange(len(data['time'])), indices_clouds_profiles)
    return indices_clouds_profiles


# main function
def get_all_flags(rawfilepath, limitez, max_calibration_height):
    '''
    Main function is used to apply all flag ways
    '''
    #-----------------
    range_limite_top = [26000,28000]
    range_limite_bottom = [2000,3000]
    dataraw = xr.open_dataset(rawfilepath)
    limitez = (dataraw['range']<limitez)

    # 355nm
    #-----------------
    channel_355 = 'rcs_12' #Analog
    Range_BckgrdCorr_Signal_355 = (dataraw[channel_355]/np.square(dataraw['range']) - dataraw['bckgrd_'+channel_355])*np.square(dataraw['range'])
    mask_crit1 = filter_profile_file(Range_BckgrdCorr_Signal_355, channel_355, range_limite_top, range_limite_bottom)
    mask_crit2 = validated_profile(Range_BckgrdCorr_Signal_355.isel(range=limitez))
    try:
        mask_crit3 = ipral_remove_cloud_profiles(max_calibration_height, rawpath)
    except:
        mask_crit3 = np.zeros((mask_crit1.shape), dtype='bool')
        pass
    mask_crit355 = mask_crit1.astype('int')*(2**1) + mask_crit2.astype('int')*(2**2) + mask_crit3.astype('int')*(2**3)

    # 532nm
    #-----------------
    channel_532 = 'rcs_16' #Analog
    Range_BckgrdCorr_Signal_532 = (dataraw[channel_532]/np.square(dataraw['range']) - dataraw['bckgrd_'+channel_532])*np.square(dataraw['range'])
    mask_crit1 = filter_profile_file(Range_BckgrdCorr_Signal_532, channel_532, range_limite_top, range_limite_bottom)
    mask_crit2 = validated_profile(Range_BckgrdCorr_Signal_532.isel(range=limitez))
    mask_crit532 = mask_crit1.astype('int')*(2**1) + mask_crit2.astype('int')*(2**2) + mask_crit3.astype('int')*(2**3)

    # intersection
    #-----------------
    mask_crits = np.intersect1d(mask_crit355.where(mask_crit355==0)['time'], mask_crit532.where(mask_crit532==0)['time'])
    return mask_crit355, mask_crit532, mask_crits