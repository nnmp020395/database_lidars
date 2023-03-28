'''
Ce script contient des fonctions d'outils qui servent au traitement mineur des données 
'''

import numpy as np 
import pandas as pd 
import math, os, sys
from pathlib import Path 



def convert_gpstime(gpstime, date, convert=False):
    '''
    Cette fonction convertit les jours du format décimal en format date-heure
    '''
    def frmt(decimal_time): # You can rewrite it with 'while' if you wish
        hours = int(decimal_time)#
        minutes = np.round(decimal_time - hours, 4)*60
        seconds = np.round(minutes - int(minutes), 4)*60
        HMS_time = f'{hours}:{int(minutes)}:{int(seconds)}'#"%s:%s:%f"%(hours, int(minutes), int(seconds))
        return HMS_time

    if convert==True:
        list_HMS_time = list(map(frmt, gpstime))
        list_YMD_HMS = list(map(lambda orig_string: date+' '+orig_string, list_HMS_time))
        pd_YMD_HMS = pd.to_datetime(list_YMD_HMS).strftime('%Y-%m-%d %H:%M:%S')
    else:
        list_gpstime_str = list(map(lambda n: '%.3f'%n, gpstime))
        list_YMD_HMS = list(map(lambda orig_string: date+' '+orig_string, list_gpstime_str))
        pd_YMD_HMS = list_YMD_HMS
    return pd.to_datetime(pd_YMD_HMS)


def find_nearest(array, value, time_option=False):
    '''
    Cette fonction permet de retrouver l'index et la valeur de array le plus proche à la valeur référence. 

    Paramètres:
    --------------
    Input: 
        array: numpy array 1D
        value: valeur référence
        time_option: si True, convertir les arrays en format datetime

    Output:
        data
        index
    '''
    if time_option:
        array = pd.to_datetime(array)
        value = pd.to_datetime(value)

    index = np.abs(array - valye).argmin()
    result = array[index]
    return index, result


def get_step(array, new_resolution):
    '''
    Cette fonction permet de retrouver le pas des éléments selon la nouvelle résolution souhaitée.

    Paramètres:
    --------------
    '''
    resolution_init = np.abs(array[3] - array[2])
    step = np.int(new_resolution / resolution_init)
    return step

    