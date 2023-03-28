import xarray as xr 
import numpy as np 
import pandas as pd 
import glob, os
import numpy as np
from datetime import datetime, timedelta
from tqdm import tqdm
from pathlib import Path

''' 
Le script sert à maintenir le processus de calibration principal pour normaliser les signaux mesurés aux signaux théoriques
Utiliser les fichiers de données temporaires homogènes
'''

# read config file 
import sys
sys.path.append('/homedata/nmpnguyen/database_lidars/Codes/')
import config_calib as cfg

# according to the instrument
if cfg.instrument == 'ipral':
    import Ipral_main_process.ensemble as ensemble
    ensemble(cfg)
elif cfg.instrument == 'er2-hsrl2':
    import ER2_main_process.ensemble as ensemble
    ensemble(cfg)

# according to the observation mode
if cfg.mode == 'nadir':
    import nadir_main_process.ensemble as ensemble
    ensemble(cfg)
elif cfg.mode == 'zenith':
    import zenith_main_process.ensemble as ensemble
    ensemble(cfg)
