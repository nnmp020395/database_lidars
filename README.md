# database_lidars
build new database for lidar on Climserv

# Outils pour Calibrer les données Lidar et créer la base de données

## Installation

Ces scripts sont faits pour être utilisés sur le mesocentre IPSL Climserv. Avant de les utiliser, il faut réaliser ces étapes. Les chemins sont spécifiques à Climserv pour une utilisation sur ciclad ou camelot, il faut réverser le dossier `homedata`.

### Téléchargement des scripts

```bash
git clone https://gitlab.in2p3.fr/ipsl/sirta/ipral/ipral-tools.git
```

### Installation de l'environnement python

- Charger une version récente de python

    ```bash
    module load python/3.8-anaconda
    ```

- Créer un dossier pour stocker les environnements

    ```bash
    mkdir /homedata/your_user/python_envs
    ```

- Créer l'environnement

    ```bash
    cd ipral-tools
    conda create -p /homedata/your_user/python_envs/env_tools --file environment.yml python=3.8 xarray=
    ```

- Activer l'environnment python

```bash
module load python/3.8-anaconda
source activate /homedata/your_user/python_envs/env_tools
```

### Installation des packages 

```bash
...
source activate /homedata/your_user/python_envs/env_tools
conda istall numpy pandas matplotlib
conda install xarray dask netCDF4 bottleneck
```

## Scripts

### convert_raw_file.py

Ce script permet de convertir les données des lidars différents en format homogène, qui va stocker les signaux bruts, les infos de canaux, altitude et temps. 

### get_simulate_data.py

Ce script utilise les coordonnées des données pour créer un base des Température et des Pressions simulées à partir des données réanalyses. 

### flag_functions.py

Ce script applique spécifiquement à Ipral pour marquer les flags, en fonction de profils de données. Ces infos vont être enregistrées dans les fichiers de sortie à la fin du processus


