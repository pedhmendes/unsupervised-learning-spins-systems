#basic imports
from ctypes.wintypes import HANDLE
from re import I
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cm
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from matplotlib.lines import Line2D

#machine learning libraries
from sklearn.decomposition import PCA
from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors

# save figs
save = 1

#show figs
show = 1

###### FUNCTIONS ######

###### FUNCTIONS ######

#With a matrix X return the pca model, eigenvalues and eigenvectors
#https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html
def get_pca(samples, N):
  pca = PCA(n_components=N)

  saux = samples.drop('T', axis = 1)
  saux = saux.drop('S', axis = 1)

  pca.fit(saux)
  evals = pca.explained_variance_ratio_
  evecs = pca.components_

  return evals, evecs, pca

#Returns the projection of matrix X in the eigenvectors of the pca model
def get_pca_results(samples, pca, N):
  saux = samples.drop('T', axis = 1)
  saux = saux.drop('S', axis = 1)
  
  pca_result = pca.transform(saux)
  pca_result = pd.DataFrame(pca_result, columns = ['C'+str(i) for i in range(0,N)])
  pca_result['T'] = samples['T']

  return pca_result

#Getting cluster with DBSCAN
def get_dbscan(A, B, EPS, MINSAMPLES):
  X_df = pd.DataFrame({'C0': A, 'C1':B})
  X = X_df.to_numpy()
  db = DBSCAN(eps=EPS, min_samples=MINSAMPLES)
  db.fit(X)
  clusters_labels = db.labels_

  n_clusters_ = len(set(clusters_labels)) - (1 if -1 in clusters_labels else 0)
  n_noise_ = list(clusters_labels).count(-1)

  return clusters_labels, n_clusters_, n_noise_

def find_eps(C0, C1, n_neigh):
    X_df = pd.DataFrame({'C0': C0, 'C1':C1})
    X = X_df.to_numpy()
    neigh = NearestNeighbors(n_neighbors=n_neigh)
    nbrs = neigh.fit(X)
    distances, indices = nbrs.kneighbors(X)
    distances = np.sort(distances, axis=0)
    distances = distances[:,1]

    return distances

def make_counts(array, bin):
    counts, bins_edge = np.histogram(array, bins=bin)
    bins = []

    for i in range(len(counts)):
        a = bins_edge[i]+bins_edge[i+1]
        a = a/2
        bins.append(int(a))
    
    df = pd.DataFrame({'counts': counts, 'bins':bins})
    df = df[df['counts'] > 0.0001]
    return df



#---------------------------------#
# pca analogue of order parameter #
#---------------------------------#
def pca_order_parameter(s1, r1, s2, r2, s3, r3):
    resultados1 = {}
    resultados2 = {}
    resultados3 = {}
    L1 = 20
    L2 = 40
    L3 = 80
    m = 50

    for temp in s1['T'].unique():

        df_ = r1.loc[r1['T'] == temp].copy()

        df_['gamminha'] = np.sqrt(np.fabs(df_['C0'])*np.fabs(df_['C0']) +
                                  np.fabs(df_['C1'])*np.fabs(df_['C1']))

        # m = len(df_) # quantidade de pontos nessa temp
        gamma = np.sum(df_['gamminha'])/(m*L1)

        resultados1[temp] = gamma

    for temp in s2['T'].unique():

        df_ = r2.loc[r2['T'] == temp].copy()

        df_['gamminha'] = np.sqrt(np.fabs(df_['C0'])*np.fabs(df_['C0']) +
                                  np.fabs(df_['C1'])*np.fabs(df_['C1']))

        # m = len(df_) # quantidade de pontos nessa temp
        gamma = np.sum(df_['gamminha'])/(m*L2)

        resultados2[temp] = gamma

    for temp in s3['T'].unique():

        df_ = r3.loc[r3['T'] == temp].copy()

        df_['gamminha'] = np.sqrt(np.fabs(df_['C0'])*np.fabs(df_['C0']) +
                                  np.fabs(df_['C1'])*np.fabs(df_['C1']))

        # m = len(df_) # quantidade de pontos nessa temp
        gamma = np.sum(df_['gamminha'])/(m*L3)

        resultados3[temp] = gamma

    return resultados1, resultados2, resultados3  


#---------------------------------#
# binder comulant                 #
#---------------------------------#
def binder_comulant(s1, r1, s2, r2, s3, r3):
    u12 = {}
    u14 = {}
    u18 = {}

    L1 = 20
    L2 = 40
    L3 = 80
    
    m = 50
  
    for temp in r1['T'].unique():
        df_ = r1.loc[r1['T'] == temp].copy()

        df_['gamminha'] = np.sqrt(np.fabs(df_['C0'])*np.fabs(df_['C0']) +
                                  np.fabs(df_['C1'])*np.fabs(df_['C1']))
      
        a2 = np.power( np.sum(np.power(df_['gamminha'],2)),2 )/(m*L1)
        a4 = np.sum( np.power(df_['gamminha'],4) )/(m*L1)

        u1 = 1 - a4/(3*a2)
        
        u12[temp] = u1
      

    for temp in r2['T'].unique():
        df_ = r2.loc[r2['T'] == temp].copy()

        df_['gamminha'] = np.sqrt(np.fabs(df_['C0'])*np.fabs(df_['C0']) +
                                  np.fabs(df_['C1'])*np.fabs(df_['C1']))
      
        a2 = np.power(np.sum(np.power(df_['gamminha'],2)),2)/(m*L2)
        a4 = np.sum(np.power(df_['gamminha'],4))/(m*L2)

        u1 = 1 - a4/(3*a2)

        u14[temp] = u1
      
    for temp in r3['T'].unique():
        df_ = r3.loc[r3['T'] == temp].copy()

        df_['gamminha'] = np.sqrt(np.fabs(df_['C0'])*np.fabs(df_['C0']) +
                                  np.fabs(df_['C1'])*np.fabs(df_['C1']))
      
        a2 = np.power(np.sum(np.power(df_['gamminha'],2)),2)/(m*L3)
        a4 = np.sum(np.power(df_['gamminha'],4))/(m*L3)

        u1 = 1 - a4/(3*a2)     

        u18[temp] = u1
        
    return u12, u14, u18




############################
###### PLOT FUNCTIONS ######
############################

def plot_eps(distances1, knee1, distances2, knee2, distances3, knee3, filename):
    #https://stackoverflow.com/questions/13583153/how-to-zoomed-a-portion-of-image-and-insert-in-the-same-plot-in-matplotlib

    fig, ax = plt.subplots(1,3, figsize=(12, 8))


    ax[0].plot(distances1)
    ax[0].plot(knee1,distances1[knee1],'ro')
    ax[0].set_xlim(4500,5100)   
    ax[0].set_xlabel('Pontos', fontsize=18)
    ax[0].set_ylabel('Distância', fontsize=18)
    axins1 = zoomed_inset_axes(ax[0], 15, loc=10) 
    axins1.plot(distances1)
    axins1.plot(knee1,distances1[knee1],'ro')
    x1, x2, y1, y2 = 4990, 5000, 4.1, 4.3
    axins1.set_xlim(x1, x2)
    axins1.set_ylim(y1, y2)
    mark_inset(ax[0], axins1, loc1=1, loc2=3, fc="none", ec="0.5")


    ax[1].plot(distances2)
    ax[1].plot(knee2,distances2[knee2],'ro')
    ax[1].set_xlim(4500,5100)   
    ax[1].set_xlabel('Pontos', fontsize=18)
    ax[1].set_ylabel('')
    axins2 = zoomed_inset_axes(ax[1], 15, loc=10) 
    axins2.plot(distances2)
    axins2.plot(knee2,distances2[knee2],'ro')
    axins2.set_xlim(x1, x2)
    axins2.set_ylim(y1, y2)
    mark_inset(ax[1], axins2, loc1=1, loc2=3, fc="none", ec="0.5")


    ax[2].plot(distances3)
    ax[2].plot(knee3,distances3[knee3],'ro')
    ax[2].set_xlim(4500,5100)   
    ax[2].set_xlabel('Pontos', fontsize=18)
    ax[2].set_ylabel('')
    axins3 = zoomed_inset_axes(ax[2], 15, loc=10) 
    axins3.plot(distances3)
    axins3.plot(knee3,distances3[knee3],'ro')
    axins3.set_xlim(x1, x2)
    axins3.set_ylim(y1, y2)
    mark_inset(ax[2], axins3, loc1=1, loc2=3, fc="none", ec="0.5")


    ax[1].set_yticks([])
    ax[1].set_yticks([])
    ax[2].set_yticklabels([])
    ax[2].set_yticklabels([])

    
    xticks = ax[0].xaxis.get_major_ticks()
    xticks[0].set_visible(False)
    xticks[-1].set_visible(False)

    xticks = ax[1].xaxis.get_major_ticks()
    xticks[0].set_visible(False)
    xticks[-1].set_visible(False)

    xticks = ax[2].xaxis.get_major_ticks()
    xticks[0].set_visible(False)
    xticks[-1].set_visible(False)
    
    plt.subplots_adjust(wspace=0)

    plt.xticks(visible=False)
    if (save==1):
      plt.savefig(filename, bbox_inches='tight')
  
    if (show==1):
      plt.show()

#Plot the eigenvalues
def plot_evals_v2(e1, e2, e3, e4, e5, e6, e7, e8, e9, filename):
  x_ax = np.arange(1, len(e1)+1, 1, dtype = 'int')
  
  fig, ax = plt.subplots(figsize=(12, 8))
  
  plt.plot(x_ax, e1, '-o', color = 'C0',markersize=9)
  plt.plot(x_ax, e2, '-*', color = 'C0',markersize=9)
  plt.plot(x_ax, e3, '-s', color = 'C0',markersize=9)

  plt.plot(x_ax, e4, '-o', color = 'C1',markersize=9)
  plt.plot(x_ax, e5, '-*', color = 'C1',markersize=9)
  plt.plot(x_ax, e6, '-s', color = 'C1',markersize=9)

  plt.plot(x_ax, e7, '-o', color = 'C2',markersize=9)
  plt.plot(x_ax, e8, '-*', color = 'C2',markersize=9)
  plt.plot(x_ax, e9, '-s', color = 'C2',markersize=9)

  plt.xlabel('$i$', fontsize=14)
  plt.xticks(np.arange(1, 11, step=1))
  plt.ylabel(r'$\tilde{\lambda}_i$', fontsize=14)
  plt.yscale("log")

  legend_elements = [Line2D([0], [0], color='C0', label='$q = 3$'),
                     Line2D([0], [0], color='C1', label='$q = 4$'),
                     Line2D([0], [0], color='C2', label='$q = 7$')
                     ]

  plt.legend(handles=legend_elements, loc='upper right', fontsize=18)


  # plt.legend(loc='best', fontsize=18)
  #plt.title("Componentes Principais")
  
  if (save==1):
    plt.savefig(filename, bbox_inches='tight')
  
  if (show==1):
    plt.show()

#Plot pca projections. 
#for more cmaps https://matplotlib.org/stable/tutorials/colors/colormaps.html
#Plot pca projections. 
#for more cmaps https://matplotlib.org/stable/tutorials/colors/colormaps.html
def plot_pca_components(A1, B1, T1, A2, B2, T2, A3, B3, T3, filename):
  #color func
  my_cmap = plt.get_cmap("rainbow")
  rescale = lambda y: (y - np.min(y)) / (np.max(y) - np.min(y))

  fig, axs = plt.subplots(1,3, figsize=(12,8), subplot_kw=dict(box_aspect=1))

  axs[0].scatter(A1, B1, color=my_cmap(rescale(T1)), edgecolor='black')
  axs[1].scatter(A2, B2, color=my_cmap(rescale(T2)), edgecolor='black')
  axs[2].scatter(A3, B3, color=my_cmap(rescale(T3)), edgecolor='black')

  axs[0].annotate('(a)', xy=(-13, 15), fontsize=18)
  axs[1].annotate('(b)', xy=(-27, 30), fontsize=18)
  axs[2].annotate('(c)', xy=(-57, 60), fontsize=18)

  for ax in axs.flat:
    ax.set_xlabel('$y_1$', fontsize=14)
    
  axs[0].set_ylabel('$y_2$', fontsize=14)

  plt.subplots_adjust(hspace=0.1)

  norm = Normalize(vmin=np.min(T1), vmax=np.max(T1)) 
  scalarMap = cm.ScalarMappable(norm=norm, cmap=my_cmap)
  fig.colorbar(scalarMap, location='top', ax=fig.get_axes())
  
  #fig.suptitle("Projeções dos Componentes Principais", x=0.45, y=.93)
  
  if (save==1):
    plt.savefig(filename, bbox_inches='tight')
  
  if (show==1):
    plt.show()

#Plot clusters
def plot_clusters(X_df, filename):
  # the desired legends
  fig, ax = plt.subplots(figsize=(12,8), subplot_kw=dict(box_aspect=1))  

  legends = ['Ruído', 'Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4','Cluster 5',
            'Cluster 5', 'Cluster 6', 'Cluster 7', 'Cluster 8']


  my_cmap = plt.get_cmap("rainbow")
  rescale = lambda y: (y - np.min(X_df['LABEL'])) / (np.max(X_df['LABEL']) - np.min(X_df['LABEL']))

  
  A = X_df[X_df['LABEL'] == -1]
  B = X_df[X_df['LABEL'] ==  0]
  C = X_df[X_df['LABEL'] ==  1]
  D = X_df[X_df['LABEL'] ==  2]
  E = X_df[X_df['LABEL'] ==  3]
  F = X_df[X_df['LABEL'] ==  4]
  G = X_df[X_df['LABEL'] ==  5]
  H = X_df[X_df['LABEL'] ==  6]
  I = X_df[X_df['LABEL'] ==  7]
  J = X_df[X_df['LABEL'] ==  8]

  plt.scatter(A['C0'], A['C1'], color=my_cmap(rescale(A['LABEL'])),
              edgecolor='black', label=legends[0])
  plt.scatter(B['C0'], B['C1'], color=my_cmap(rescale(B['LABEL'])),
              edgecolor='black', label=legends[1])
  plt.scatter(C['C0'], C['C1'], color=my_cmap(rescale(C['LABEL'])),
              edgecolor='black', label=legends[2])
  plt.scatter(D['C0'], D['C1'], color=my_cmap(rescale(D['LABEL'])),
              edgecolor='black', label=legends[3])
  plt.scatter(E['C0'], E['C1'], color=my_cmap(rescale(E['LABEL'])),
              edgecolor='black', label=legends[4])
  plt.scatter(F['C0'], F['C1'], color=my_cmap(rescale(F['LABEL'])),
              edgecolor='black', label=legends[5])
  plt.scatter(G['C0'], G['C1'], color=my_cmap(rescale(G['LABEL'])),
              edgecolor='black', label=legends[6])
  plt.scatter(H['C0'], H['C1'], color=my_cmap(rescale(H['LABEL'])),
              edgecolor='black', label=legends[7])
  plt.scatter(I['C0'], I['C1'], color=my_cmap(rescale(I['LABEL'])),
              edgecolor='black', label=legends[8])
  plt.scatter(J['C0'], J['C1'], color=my_cmap(rescale(J['LABEL'])),
              edgecolor='black', label=legends[9])
  
  ax.annotate('(c)', xy=(-57, 60), fontsize=18)

  
  plt.ylabel(r'$y_2$', fontsize=14)
  plt.xlabel(r'$y_1$', fontsize=14)
    
  #plt.title("Clusters das Projeções dos Componentes")
  # plt.legend(loc='best', fontsize = 18)
 
  if (save==1):
      plt.savefig(filename, bbox_inches='tight')

  if (show==1):
      plt.show()

def multi_plot_snapshot(datafile, N, filename):


  t = []
  t.append(round(datafile.iloc[N[0]]['T'],3))
  t.append(round(datafile.iloc[N[1]]['T'],3))
  t.append(round(datafile.iloc[N[2]]['T'],3))


  datafile = datafile.drop('T', axis = 1)  
  datafile = datafile.drop('S', axis = 1)
  datafile = datafile.drop('Q', axis = 1)

  #plot
  fig, axs = plt.subplots(1,3, figsize=(15,15))
  
  aa = datafile.iloc[N[0]]
  bb = datafile.iloc[N[1]]
  cc = datafile.iloc[N[2]]

  aa = np.array(aa)
  bb = np.array(bb)
  cc = np.array(cc)

  aa = aa.reshape((80,80))
  bb = bb.reshape((80,80))
  cc = cc.reshape((80,80))

  axs[0].set_title('T = '+ str(t[0]), fontsize=18)
  axs[0].imshow(aa, cmap='gist_rainbow')
  axs[0].tick_params(left = False, right = False,
                     labelleft=False, labelbottom = False,
                     bottom = False)

  axs[1].set_title('T = '+ str(t[1]), fontsize=18)
  axs[1].imshow(bb, cmap='gist_rainbow')
  axs[1].tick_params(left = False, right = False,
                     labelleft=False, labelbottom = False,
                     bottom = False)

  axs[2].set_title('T = '+ str(t[2]), fontsize=18)
  axs[2].imshow(cc, cmap='gist_rainbow')
  axs[2].tick_params(left = False, right = False,
                     labelleft=False, labelbottom = False,
                     bottom = False)


  plt.subplots_adjust(hspace=0.1, wspace=0.1)

  if (save==1):
    plt.savefig(filename, bbox_inches='tight')

  if (show==1):
    plt.show()

  return

def series_temp(A, B, ylabel, title, filename, cor):
    fig, axs = plt.subplots(figsize=(12,8))
    
    ymin,ymax = 450000,550000
    plt.xlim(ymin, ymax)

    plt.scatter(A, B, label=title, color=cor, alpha=0.7, s=10)
    
    plt.xlabel('MCS', fontsize = 32)
    plt.ylabel(ylabel, fontsize = 32)
    axs.tick_params(axis='both', which='major', labelsize=18)
    axs.tick_params(axis='both', which='minor', labelsize=18)
    # plt.legend(loc='upper right', fontsize=16)
    # plt.title("Serie Temporal da " + title)
    # axs.xaxis.set_major_formatter(ticker.FuncFormatter(myticks))

    if (save==1):
        plt.savefig(filename, bbox_inches='tight')

    if (show==1):
        plt.show()

    return

def plot_order_parameter(e11, e12, e13, e21, e22, e23, e31, e32, e33, filename):
  fig, ax = plt.subplots(3, 1, figsize=(12, 8))
  
  ax[0].plot(e11.keys(), e11.values(), '-p', label='$L = 20$')
  ax[0].plot(e12.keys(), e12.values(), '-p', label='$L = 40$')
  ax[0].plot(e13.keys(), e13.values(), '-p', label='$L = 80$')

  ax[1].plot(e21.keys(), e21.values(), '-p')
  ax[1].plot(e22.keys(), e22.values(), '-p')
  ax[1].plot(e23.keys(), e23.values(), '-p')

  ax[2].plot(e31.keys(), e31.values(), '-p')
  ax[2].plot(e32.keys(), e32.values(), '-p')
  ax[2].plot(e33.keys(), e33.values(), '-p')

  ax[0].set_xlim(0.895,1.095)
  ax[1].set_xlim(0.81,1.01)
  ax[2].set_xlim(0.673,0.873)

  ax[0].annotate('(a)', xy=(0.95, 0.6), fontsize=18)
  ax[1].annotate('(b)', xy=(0.86, 0.6), fontsize=18)
  ax[2].annotate('(c)', xy=(0.725, 0.6), fontsize=18)

  for a in ax.flat:
    a.set_ylabel(r'$\Gamma(T)$', fontsize=14)
    
  ax[2].set_xlabel('$T$', fontsize=14)

  legend_elements = [Line2D([0], [0], color='C0', label='$L = 20$'),
                     Line2D([0], [0], color='C1', label='$L = 40$'),
                     Line2D([0], [0], color='C2', label='$L = 80$')]

  ax[0].legend(handles=legend_elements, loc='lower left', fontsize=18)
  # ax[0].legend(loc='lower left', fontsize=18)
  # ax[1].legend(loc='lower left', fontsize=18)
  # ax[2].legend(loc='lower left', fontsize=18)

  if (save==1):
    plt.savefig(filename, bbox_inches='tight')
  
  if (show==1):
    plt.show()

def plot_binder(e11, e12, e13, lim1, e21, e22, e23, lim2, e31, e32, e33, lim3, filename):
  fig, ax = plt.subplots(1,3, figsize=(12, 8))
  
  ax[0].plot(e11.keys(), e11.values(), '-^', label='$L = 20$')
  ax[0].plot(e12.keys(), e12.values(), '-^', label='$L = 40$')
  ax[0].plot(e13.keys(), e13.values(), '-^', label='$L = 80$')

  ax[1].plot(e21.keys(), e21.values(), '-^')
  ax[1].plot(e22.keys(), e22.values(), '-^')
  ax[1].plot(e23.keys(), e23.values(), '-^')

  ax[2].plot(e31.keys(), e31.values(), '-^')
  ax[2].plot(e32.keys(), e32.values(), '-^')
  ax[2].plot(e33.keys(), e33.values(), '-^')

  ax[0].annotate('(a)', xy=(0.92, 0.99375), fontsize=18)
  ax[1].annotate('(b)', xy=(0.84, 0.99375), fontsize=18)
  ax[2].annotate('(c)', xy=(0.70, 0.99375), fontsize=18)

  for a in ax.flat:
    a.set_xlabel('$T$', fontsize=14)
   

  xticks = ax[0].xaxis.get_major_ticks()
  xticks[0].set_visible(False)
  xticks[-1].set_visible(False)

  xticks = ax[1].xaxis.get_major_ticks()
  xticks[0].set_visible(False)
  xticks[-1].set_visible(False)

  xticks = ax[2].xaxis.get_major_ticks()
  xticks[0].set_visible(False)
  xticks[-1].set_visible(False)
  
  plt.subplots_adjust(wspace=0)

  ax[0].set_ylabel('$B(T,L)$', fontsize=14)

  ax[1].set_yticks([])
  ax[1].set_yticks([])
  ax[2].set_yticklabels([])
  ax[2].set_yticklabels([])

  legend_elements = [Line2D([0], [0], color='C0', label='$L = 20$'),
                     Line2D([0], [0], color='C1', label='$L = 40$'),
                     Line2D([0], [0], color='C2', label='$L = 80$')]

  ax[0].legend(handles=legend_elements, loc='lower left', fontsize=18)

  # ax[0].legend(loc='best', fontsize=18)
  # ax[1].legend(loc='best', fontsize=18)
  # ax[2].legend(loc='best', fontsize=18)

  ax[0].set_xlim(lim1[0],lim1[1])
  ax[0].set_ylim(lim1[2],lim1[3])
  
  ax[1].set_xlim(lim2[0],lim2[1])
  ax[1].set_ylim(lim2[2],lim2[3])
  
  ax[2].set_xlim(lim3[0],lim3[1])
  ax[2].set_ylim(lim3[2],lim3[3])
  


  if (save==1):
    plt.savefig(filename, bbox_inches='tight')
  
  if (show==1):
    plt.show()

