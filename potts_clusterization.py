# basic libraries
from re import T
import pandas as pd
import numpy as np

# import funcs
from libs import get_dbscan
from libs import find_eps, plot_eps
from libs import plot_clusters
from libs import get_pca, get_pca_results

# removing warnings
import warnings
warnings.filterwarnings("ignore")

# definitions

potts_samples80 = pd.read_parquet("E:\TCC\Potts/results_ok/q3-df80-final.gzip")
potts_samples802 = pd.read_parquet("E:\TCC\Potts/results_ok/q4-df80-final.gzip")
potts_samples803 = pd.read_parquet("E:\TCC\Potts/results_ok/q7-df80-final.gzip")

N = 10 
Nmin = 25

#---------------------#
#DEFS
#---------------------#

evals80, evecs80, pca80 = get_pca(potts_samples80, N)
pca_r80 = get_pca_results(potts_samples80, pca80, N)
evals80, evecs80, pca80 = get_pca(potts_samples802, N)
pca_r802 = get_pca_results(potts_samples802, pca80, N)
evals80, evecs80, pca80 = get_pca(potts_samples803, N)
pca_r803 = get_pca_results(potts_samples803, pca80, N)

# Defs for Clusterization
C0 = pca_r80['C0']
C1 = pca_r80['C1']
TEMP = pca_r80['T']

C02 = pca_r802['C0']
C12 = pca_r802['C1']

C03 = pca_r803['C0']
C13 = pca_r803['C1']

#---------------------#
#eps
#---------------------#
distances = find_eps(C0, C1, Nmin)
distances2 = find_eps(C02, C12, Nmin)
distances3 = find_eps(C03, C13, Nmin)

from kneed import KneeLocator

x_ax = np.arange(1, len(distances)+1, 1, dtype = 'int')
x_ax2 = np.arange(1, len(distances2)+1, 1, dtype = 'int')
x_ax3 = np.arange(1, len(distances3)+1, 1, dtype = 'int')

kneedle = KneeLocator(x_ax, distances, S=1.0, curve="convex", direction="increasing")#, interp_method='polynomial')
print(round(kneedle.knee, 3), distances[round(kneedle.knee, 3)])

kneedle2 = KneeLocator(x_ax2, distances2, S=1.0, curve="convex", direction="increasing")#, interp_method='polynomial')
print(round(kneedle2.knee, 3), distances2[round(kneedle2.knee, 3)])

kneedle3 = KneeLocator(x_ax3, distances3, S=1.0, curve="convex", direction="increasing")#, interp_method='polynomial')
print(round(kneedle3.knee, 3), distances3[round(kneedle3.knee, 3)])

plot_eps(distances, round(kneedle.knee, 3),
         distances2, round(kneedle2.knee, 3),
         distances3, round(kneedle3.knee, 3),
         'E:\TCC\Potts\plots\eps_plot.png')

#---------------------#
#DBSCAN
#---------------------#

# Passing DBSCAN
eps = distances[round(kneedle.knee, 3)]
dbc_lb, dbc_nc, dbc_nn = get_dbscan(C0, C1, eps, Nmin)
X = pd.DataFrame({'C0': C0, 'C1':C1, 'LABEL': dbc_lb})

#number of clusters
print("Number of clusters:", dbc_nc)

#number of noise
print("Number of noise points:", dbc_nn)

# ploting clusters and noise
plot_clusters(X, 'E:\TCC\Potts\plots\ising_clusters_db_L80-q3.png')


# Passing DBSCAN
eps = distances2[round(kneedle2.knee, 3)]
dbc_lb, dbc_nc, dbc_nn = get_dbscan(C02, C12, eps, Nmin)
X = pd.DataFrame({'C0': C02, 'C1':C12, 'LABEL': dbc_lb})

#number of clusters
print("Number of clusters:", dbc_nc)

#number of noise
print("Number of noise points:", dbc_nn)

# ploting clusters and noise
plot_clusters(X, 'E:\TCC\Potts\plots\ising_clusters_db_L80-q4.png')


# Passing DBSCAN
eps = distances3[round(kneedle3.knee, 3)]
dbc_lb, dbc_nc, dbc_nn = get_dbscan(C03, C13, eps, Nmin)
X = pd.DataFrame({'C0': C03, 'C1':C13, 'LABEL': dbc_lb})

#number of clusters
print("Number of clusters:", dbc_nc)

#number of noise
print("Number of noise points:", dbc_nn)

# ploting clusters and noise
plot_clusters(X, 'E:\TCC\Potts\plots\ising_clusters_db_L80-q7.png')