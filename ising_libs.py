#basic imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cm
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

#machine learning libraries
from sklearn.decomposition import PCA
from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors

# save figs
save = 1  

#show figs
show = 1

###### FUNCTIONS ######

#With a matrix X return the pca model, eigenvalues and eigenvectors
#https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html
def get_pca(samples, N):
  pca = PCA(n_components=N)

  saux = samples.drop('T', axis = 1)
  saux = saux.drop('M', axis = 1)
  saux = saux.drop('E', axis = 1)
  saux = saux.drop('S', axis = 1)

  pca.fit(saux)
  evals = pca.explained_variance_ratio_
  evecs = pca.components_

  return evals, evecs, pca

#Returns the projection of matrix X in the eigenvectors of the pca model
def get_pca_results(samples, pca, N):
  saux = samples.drop('T', axis = 1)
  saux = saux.drop('M', axis = 1)
  saux = saux.drop('E', axis = 1)
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
    df = df[df['counts'] > 0.001]
    return df


############################
###### PLOT FUNCTIONS ######
############################

def plot_eps(distances, knee, filename):
    #https://stackoverflow.com/questions/13583153/how-to-zoomed-a-portion-of-image-and-insert-in-the-same-plot-in-matplotlib

    fig, ax = plt.subplots(figsize=(12,8))  
    
    #kneedle = KneeLocator(x, y, S=1.0, curve="concave", direction="increasing")

    #plt.title("$\epsilon$ em função número de vizinhos")

    ax.plot(distances, label='Distância entre pontos')
    ax.plot(knee,distances[knee],'ro')
    ax.set_xlabel('Pontos', fontsize=18)
    ax.set_ylabel('Distância', fontsize=18)
    
    axins = zoomed_inset_axes(ax, 15, loc=10) 
    axins.plot(distances)
    axins.plot(knee,distances[knee],'ro')

    x1, x2, y1, y2 = 12990, 13100, 3, 3.3
    axins.set_xlim(x1, x2)
    axins.set_ylim(y1, y2)

  
    axins.annotate('$\epsilon$ = 0,75',
            xy=(0.85, 13000),
            xycoords='data',
            xytext=(1, 13200),
            arrowprops=
                dict(facecolor='black', shrink=0.05),
                horizontalalignment='left',
                verticalalignment='top')

    plt.xticks(visible=False)
    mark_inset(ax, axins, loc1=1, loc2=3, fc="none", ec="0.5")

    # ax.legend(loc='best', fontsize=18)

    if (save==1):
      plt.savefig(filename, bbox_inches='tight')
  
    if (show==1):
      plt.show()

#Plot the eigenvalues
def plot_evals(e1, e2, e3, filename):
  x_ax = np.arange(1, len(e1)+1, 1, dtype = 'int')
  
  fig, ax = plt.subplots(figsize=(12, 8))
  
  plt.plot(x_ax, e1, '-o', label='$L = 20$')
  plt.plot(x_ax, e2, '-o', label='$L = 40$')
  plt.plot(x_ax, e3, '-o', label='$L = 80$')

  plt.xlabel('$i$', fontsize=18)
  plt.xticks(np.arange(1, 11, step=1))
  plt.ylabel(r'$\tilde{\lambda}_i$', fontsize=18)
  plt.yscale("log")

  plt.legend(loc='best', fontsize=18)
  #plt.title("Componentes Principais")
  
  if (save==1):
    plt.savefig(filename, bbox_inches='tight')
  
  if (show==1):
    plt.show()


#Plot pca projections. 
#for more cmaps https://matplotlib.org/stable/tutorials/colors/colormaps.html
def plot_pca_components(A1, B1, T1, A2, B2, T2, A3, B3, T3, filename):
  #color func
  my_cmap = plt.get_cmap("rainbow")
  rescale = lambda y: (y - np.min(y)) / (np.max(y) - np.min(y))
    
  fig, axs = plt.subplots(3, figsize=(12,8))

  axs[0].scatter(A1, B1, color=my_cmap(rescale(T1)), edgecolor='black')
  axs[1].scatter(A2, B2, color=my_cmap(rescale(T2)), edgecolor='black')
  axs[2].scatter(A3, B3, color=my_cmap(rescale(T3)), edgecolor='black')

  for ax in axs.flat:
    ax.set_ylabel('$y_2$', fontsize=18)
    ax.set_xlim([np.min(A3)-10, np.max(A3)+10])
    ax.set_ylim([np.min(B3)-10, np.max(B3)+10])

  axs[0].set_xticks([])
  axs[1].set_xticks([])
  axs[0].set_xticklabels([])
  axs[1].set_xticklabels([])

  axs[2].set_xlabel('$y_1$', fontsize=18)

  axs[0].annotate('(a)', xy=(-80, 40), fontsize=18)
  axs[1].annotate('(b)', xy=(-80, 40), fontsize=18)
  axs[2].annotate('(c)', xy=(-80, 40), fontsize=18)

  plt.subplots_adjust(hspace=0.0)

  norm = Normalize(vmin=np.min(T1), vmax=np.max(T1)) 
  scalarMap = cm.ScalarMappable(norm=norm, cmap=my_cmap)
  fig.colorbar(scalarMap, ax=fig.get_axes())
  
  #fig.suptitle("Projeções dos Componentes Principais", x=0.45, y=.93)
  
  if (save==1):
    plt.savefig(filename, bbox_inches='tight')
  
  if (show==1):
    plt.show()

#Plot clusters
def plot_clusters(X_df, filename):
  # the desired legends
  fig, ax = plt.subplots(figsize=(12,8))  

  legends = ['Ruído', 'Cluster 1', 'Cluster 2', 'Cluster 3']

  my_cmap = plt.get_cmap("rainbow")
  rescale = lambda y: (y - np.min(X_df['LABEL'])) / (np.max(X_df['LABEL']) - np.min(X_df['LABEL']))

  A = X_df[X_df['LABEL'] == -1]
  B = X_df[X_df['LABEL'] ==  0]
  C = X_df[X_df['LABEL'] ==  1]
  D = X_df[X_df['LABEL'] ==  2]

  plt.scatter(A['C0'], A['C1'], color=my_cmap(rescale(A['LABEL'])),
              edgecolor='black', label=legends[0])
  plt.scatter(B['C0'], B['C1'], color=my_cmap(rescale(B['LABEL'])),
              edgecolor='black', label=legends[1])
  plt.scatter(C['C0'], C['C1'], color=my_cmap(rescale(C['LABEL'])),
              edgecolor='black', label=legends[2])
  plt.scatter(D['C0'], D['C1'], color=my_cmap(rescale(D['LABEL'])),
              edgecolor='black', label=legends[3])              
  
  plt.ylabel(r'$y_2$', fontsize=18)
  plt.xlabel(r'$y_1$', fontsize=18)
    
  #plt.title("Clusters das Projeções dos Componentes")
  #plt.legend(loc='best', fontsize = 18)
 
  if (save==1):
      plt.savefig(filename, bbox_inches='tight')

  if (show==1):
      plt.show()

#Plot T series
def plot_tss(A, B, ylabel, title, filename, cor):
  fig, axs = plt.subplots(figsize=(12,8))
  
  min = np.min(B) -1000
  max = np.max(B) +1000

  plt.ylim(min, max)

  plt.scatter(A, B, label=title, alpha=.7, edgecolor='black', color=cor)
  plt.vlines(2.269, ymin=min, ymax=max, color='black')

  plt.xlabel('$T$', fontsize = 32)
  plt.ylabel(ylabel, fontsize = 32)
  axs.tick_params(axis='both', which='major', labelsize=18)
  axs.tick_params(axis='both', which='minor', labelsize=18)
  #plt.legend(loc='best', fontsize=16)
  # plt.title("Serie Temporal da " + title)

  save = 1
  show = 1

  if (save==1):
    plt.savefig(filename, bbox_inches='tight')

  if (show==1):
    plt.show()

  return

#Plot T series with cluster label colors
def plot_ts_clusters(A, B, clusters_labels, ylabel, title, filename):
  fig, axs = plt.subplots(figsize=(12,8))
  
  min = np.min(B) - 10
  max = np.max(B) + 10

  plt.ylim(min, max)

  my_cmap = plt.get_cmap("rainbow")
  rescale = lambda y: (y - np.min(y)) / (np.max(y) - np.min(y))

  plt.scatter(A, B, color=my_cmap(rescale(clusters_labels)), edgecolor='black', label=title)
  #plt.vlines(2.269, ymin=min, ymax=max, color='black')

  plt.xlabel('$T$', fontsize=18)
  plt.ylabel(ylabel, fontsize=18)
  plt.plot(legend=None)
  #axs.get_legend().remove()

  #leg = axs.legend(handlelength=0, handletextpad=0, fancybox=True, fontsize=14)
  #for item in leg.legendHandles:
    #item.set_visible(False)

  # plt.legend(loc='best')
  #plt.title("Serie Temporal de " + ylabel)

  if (save==1):
    plt.savefig(filename, bbox_inches='tight')

  if (show==1):
    plt.show()

def multi_plot_snapshot(datafile, filename):
  #N = [19,42]
  N = [8000,42,3]

  t = []
  t.append(datafile.iloc[N[0]]['T'])
  t.append(datafile.iloc[N[1]]['T'])
  t.append(datafile.iloc[N[2]]['T'])


  datafile = datafile.drop('T', axis = 1)  
  datafile = datafile.drop('M', axis = 1)
  datafile = datafile.drop('E', axis = 1)
  datafile = datafile.drop('S', axis = 1)

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
  axs[0].imshow(aa, cmap='Greys')
  axs[0].tick_params(left = False, right = False,
                     labelleft=False, labelbottom = False,
                     bottom = False)

  axs[1].set_title('T = '+ str(t[1]), fontsize=18)
  axs[1].imshow(bb, cmap='Greys')
  axs[1].tick_params(left = False, right = False,
                     labelleft=False, labelbottom = False,
                     bottom = False)

  axs[2].set_title('T = '+ str(t[2]), fontsize=18)
  axs[2].imshow(cc, cmap='Greys')
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
    #plt.legend(loc='upper right', fontsize=16)
    # plt.title("Serie Temporal da " + title)
      # plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))


    if (save==1):
        plt.savefig(filename, bbox_inches='tight')

    if (show==1):
        plt.show()

    return





