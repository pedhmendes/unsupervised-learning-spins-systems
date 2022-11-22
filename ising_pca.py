#basic libraries
import pandas as pd

#import funcs
from libs import get_pca, get_pca_results
from libs import plot_evals, plot_pca_components, plot_ts_clusters

#removing warnings
import warnings
warnings.filterwarnings("ignore")

#definitions
L20 = 20      #LSIZE
data_path20 = "E:\TCC\Ising\df20.csv"
L40 = 40      #LSIZE
data_path40 = "E:\TCC\Ising\df40.csv"
L80 = 80      #LSIZE
data_path80 = "E:\TCC\Ising\df80.csv"

#reading and dropping first column
ising_samples20 = pd.read_csv(data_path20)
ising_samples40 = pd.read_csv(data_path40)
ising_samples80 = pd.read_csv(data_path80)

#Number of components
N = 10 

#PCA calculation
evals20, evecs20, pca20 = get_pca(ising_samples20, N)
evals40, evecs40, pca40 = get_pca(ising_samples40, N)
evals80, evecs80, pca80 = get_pca(ising_samples80, N)

#print variance ratio 
#sum of 2 first pc
ratio20 = pca20.explained_variance_ratio_[0]+pca20.explained_variance_ratio_[1]
ratio40 = pca40.explained_variance_ratio_[0]+pca40.explained_variance_ratio_[1]
ratio80 = pca80.explained_variance_ratio_[0]+pca80.explained_variance_ratio_[1]

print(pca20.explained_variance_ratio_[0],pca20.explained_variance_ratio_[1],ratio20)
print(pca40.explained_variance_ratio_[0],pca40.explained_variance_ratio_[1],ratio40)
print(pca80.explained_variance_ratio_[0],pca80.explained_variance_ratio_[1],ratio80)

#Get the projections
pca_r_20 = get_pca_results(ising_samples20, pca20, N)
pca_r_40 = get_pca_results(ising_samples40, pca40, N)
pca_r_80 = get_pca_results(ising_samples80, pca80, N)

#PCA eigenvalues plot, Y axis in log scale
plot_evals(evals20, evals40, evals80, 'E:\TCC\Ising\plots\ising_pca_multi_evals.png')

#Plotting the projections
plot_pca_components(pca_r_20['C0'], pca_r_20['C1'], pca_r_20['T'], 
                    pca_r_40['C0'], pca_r_40['C1'], pca_r_40['T'], 
                    pca_r_80['C0'], pca_r_80['C1'], pca_r_80['T'], 
                    'E:\TCC\Ising\plots\ising_pca_components.png')
plot_ts_clusters(pca_r_80['T'], pca_r_80['C0'],
                 pca_r_80['T'], '$y_1$', '$y_1$',
                 'E:\TCC\Ising\plots\y1-temp-L80.png')

