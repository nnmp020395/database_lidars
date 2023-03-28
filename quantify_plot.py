'''
Ce script permet de tracer le scatter dont l'axe y représente l'altitude, l'axe x représente les mois, 
le colorbar représente soit la fraction des profils validés sur le nombre des profils total, 
soit le nombre total des profils validés, soit le nombre total des profils du dataset
'''
import numpy as np 
import pandas as pd 
from pathlib import Path 
import matplotlib.pyplot as plt 

def quantify_plot(x, y, data, colorbar_bounds, lims, labels, title, suptitle, outputfile, figsize):
	'''
	Input 
	___________

	x : array1D pour l'axe x (= mois)
	y : array1D pour l'axe y (= altitude)
	data : array2D pour l'affichage 
	colorbar_bounds : l'intervalle des valeurs pour le colorbar
	lims : dictionnaire de limites 'xlim' et 'ylim'
	labels : dictionnaire des titres des axes 'x', 'y'
	title : titre de la figure
	suptitle : grand titre de la figure
	outputfile : chemin d'enregistrement pour la figure
	figsize : taille souhaité de la figure 

	Output
	___________

	plot 
	
	'''
    cmap = mpl.cm.turbo
#     bounds = np.arange(colorbar_bounds)#[-1, 2, 5, 7, 12, 15]
    norm = mpl.colors.BoundaryNorm(colorbar_bounds, cmap.N, extend='max')

    fig, ax = plt.subplots(figsize=figsize)
    p=ax.pcolormesh(x, y, data.T, norm=norm, cmap=cmap, shading='auto')
    cbar = plt.colorbar(p, ax=ax, orientation='horizontal', label= labels['colorbar'])
    cbar.ax.set_xticklabels(np.round(cbar.get_ticks(), decimals=2), rotation=45)
    
    plt.ylim(lims['ylim'])
    plt.xlim(lims['xlim'])
    plt.xlabel(labels['x'])
    plt.ylabel(labels['y']) 
#     ax.set_xticklabels(np.arange(1,13,1))
    plt.title(title, loc='right', fontsize=11)
    plt.suptitle(suptitle)
#     ax.grid(which='major', color='k', linestyle='-')
#     ax.grid(which='minor', color='k', linestyle='-')
#     ax.grid(True, color='k')
    plt.savefig(outputfile)