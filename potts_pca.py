#basic libraries
import pandas as pd

#import funcs
from libs import get_pca, get_pca_results
from libs import plot_pca_components, plot_order_parameter, plot_binder
from libs import pca_order_parameter, binder_comulant, plot_evals_v2

#removing warnings
import warnings
warnings.filterwarnings("ignore")

# import csv
q3_potts_samples20 = pd.read_parquet("E:\TCC\Potts/results_ok/q3-df20-final.gzip")
q3_potts_samples40 = pd.read_parquet("E:\TCC\Potts/results_ok/q3-df40-final.gzip")
q3_potts_samples80 = pd.read_parquet("E:\TCC\Potts/results_ok/q3-df80-final.gzip")

q4_potts_samples20 = pd.read_parquet("E:\TCC\Potts/results_ok/q4-df20-final.gzip")
q4_potts_samples40 = pd.read_parquet("E:\TCC\Potts/results_ok/q4-df40-final.gzip")
q4_potts_samples80 = pd.read_parquet("E:\TCC\Potts/results_ok/q4-df80-final.gzip")

q7_potts_samples20 = pd.read_parquet("E:\TCC\Potts/results_ok/q7-df20-final.gzip")
q7_potts_samples40 = pd.read_parquet("E:\TCC\Potts/results_ok/q7-df40-final.gzip")
q7_potts_samples80 = pd.read_parquet("E:\TCC\Potts/results_ok/q7-df80-final.gzip")


#Number of components
N = 10 

#########Q3
#PCA calculation
q3_evals20, q3_evecs20, q3_pca20 = get_pca(q3_potts_samples20, N)
q3_evals40, q3_evecs40, q3_pca40 = get_pca(q3_potts_samples40, N)
q3_evals80, q3_evecs80, q3_pca80 = get_pca(q3_potts_samples80, N)

#Get the projections
q3_pca_r20 = get_pca_results(q3_potts_samples20, q3_pca20, N)
q3_pca_r40 = get_pca_results(q3_potts_samples40, q3_pca40, N)
q3_pca_r80 = get_pca_results(q3_potts_samples80, q3_pca80, N)

#order parameter
q3_gamma1, q3_gamma2, q3_gamma3 = pca_order_parameter(q3_potts_samples20, q3_pca_r20,
                                                      q3_potts_samples40, q3_pca_r40,
                                                      q3_potts_samples80, q3_pca_r80)

#binder parameter
q3_u12, q3_u14, q3_u18 = binder_comulant(q3_potts_samples20, q3_pca_r20,
                                         q3_potts_samples40, q3_pca_r40,
                                         q3_potts_samples80, q3_pca_r80)

#########Q4
#PCA calculation
q4_evals20, q4_evecs20, q4_pca20 = get_pca(q4_potts_samples20, N)
q4_evals40, q4_evecs40, q4_pca40 = get_pca(q4_potts_samples40, N)
q4_evals80, q4_evecs80, q4_pca80 = get_pca(q4_potts_samples80, N)

#Get the projections
q4_pca_r20 = get_pca_results(q4_potts_samples20, q4_pca20, N)
q4_pca_r40 = get_pca_results(q4_potts_samples40, q4_pca40, N)
q4_pca_r80 = get_pca_results(q4_potts_samples80, q4_pca80, N)

#order parameter
q4_gamma1, q4_gamma2, q4_gamma3 = pca_order_parameter(q4_potts_samples20, q4_pca_r20,
                                                      q4_potts_samples40, q4_pca_r40,
                                                      q4_potts_samples80, q4_pca_r80)

#binder parameter
q4_u12, q4_u14, q4_u18 = binder_comulant(q4_potts_samples20, q4_pca_r20,
                                         q4_potts_samples40, q4_pca_r40,
                                         q4_potts_samples80, q4_pca_r80)

#########Q7
#PCA calculation
q7_evals20, q7_evecs20, q7_pca20 = get_pca(q7_potts_samples20, N)
q7_evals40, q7_evecs40, q7_pca40 = get_pca(q7_potts_samples40, N)
q7_evals80, q7_evecs80, q7_pca80 = get_pca(q7_potts_samples80, N)

#Get the projections
q7_pca_r20 = get_pca_results(q7_potts_samples20, q7_pca20, N)
q7_pca_r40 = get_pca_results(q7_potts_samples40, q7_pca40, N)
q7_pca_r80 = get_pca_results(q7_potts_samples80, q7_pca80, N)

#order parameter
q7_gamma1, q7_gamma2, q7_gamma3 = pca_order_parameter(q7_potts_samples20, q7_pca_r20,
                                                      q7_potts_samples40, q7_pca_r40,
                                                      q7_potts_samples80, q7_pca_r80)

#binder parameter
q7_u12, q7_u14, q7_u18 = binder_comulant(q7_potts_samples20, q7_pca_r20,
                                         q7_potts_samples40, q7_pca_r40,
                                         q7_potts_samples80, q7_pca_r80)


############PLOTS#############

plot_evals_v2(q3_evals20, q3_evals40, q3_evals80,
              q4_evals20, q4_evals40, q4_evals80,
              q7_evals20, q7_evals40, q7_evals80,
              'E:\TCC\Potts\plots\potts_pca_multi_evals-v2.png')


#Plotting the projections
plot_pca_components(q3_pca_r20['C0'], q3_pca_r20['C1'], q3_potts_samples20['T'],
                    q3_pca_r40['C0'], q3_pca_r40['C1'], q3_potts_samples40['T'],
                    q3_pca_r80['C0'], q3_pca_r80['C1'], q3_potts_samples80['T'],
                    'E:\TCC\Potts\plots\potts_pca_components-q3.png')

plot_pca_components(q4_pca_r20['C0'], q4_pca_r20['C1'], q4_potts_samples20['T'],
                    q4_pca_r40['C0'], q4_pca_r40['C1'], q4_potts_samples40['T'],
                    q4_pca_r80['C0'], q4_pca_r80['C1'], q4_potts_samples80['T'],
                    'E:\TCC\Potts\plots\potts_pca_components-q4.png')

plot_pca_components(q7_pca_r20['C0'], q7_pca_r20['C1'], q7_potts_samples20['T'],
                    q7_pca_r40['C0'], q7_pca_r40['C1'], q7_potts_samples40['T'],
                    q7_pca_r80['C0'], q7_pca_r80['C1'], q7_potts_samples80['T'],
                    'E:\TCC\Potts\plots\potts_pca_components-q7.png')

# Plot order parameter
plot_order_parameter(q3_gamma1, q3_gamma2, q3_gamma3,
                     q4_gamma1, q4_gamma2, q4_gamma3,
                     q7_gamma1, q7_gamma2, q7_gamma3,
                     r'E:\TCC\Potts\plots\order_par.png')

q3_lim=[0.9,1.02,0.9875,0.994] 
q4_lim=[0.83,0.9225,0.9875,0.994] 
q7_lim=[0.68,0.79,0.9875,0.994] 

# Plot Binder
plot_binder(q3_u12, q3_u14, q3_u18, q3_lim,
            q4_u12, q4_u14, q4_u18, q4_lim,
            q7_u12, q7_u14, q7_u18, q7_lim,
            r'E:\TCC\Potts\plots\binder.png')