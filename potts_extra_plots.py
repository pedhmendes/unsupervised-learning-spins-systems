#basic libraries
import pandas as pd
from sys import exit
import matplotlib.pyplot as plt

#imports
from libs import make_counts, series_temp, multi_plot_snapshot

#removing warnings
import warnings
warnings.filterwarnings("ignore")

#---------------------------------------------------#
# extra plots                                       #
#---------------------------------------------------#
#---------------------------------------------------#
# repesa beale                                      #
#---------------------------------------------------#
#---------------------------------------------------#
scale = 1e6
#---------------------------------------------------#

beale32 = "E:\TCC\Potts\e32.dat"
beale64 = "E:\TCC\Potts\e64.dat"
seriepath_32 = "E:\TCC\Potts\potts_Q2_T1.100000_L32_S1663620779.dsf"
seriepath_64 = "E:\TCC\Potts\potts_Q2_T1.100000_L64_S1663622971.dsf"
datapath80 =  "E:\TCC\Potts\q3-df80.csv"

potts_samples80 = pd.read_csv("E:\TCC\Potts/results/q3-df80.csv")
s32 = pd.read_csv(seriepath_32, sep='\t', index_col=False)
s64 = pd.read_csv(seriepath_64, sep='\t', index_col=False)
b32 = pd.read_csv(beale32, sep=' ', index_col=False)
b64 = pd.read_csv(beale64, sep=' ', index_col=False)

#---------------------------#
zoom = 0
if (zoom == 1):
    s32 = s32[::10]
    series_temp(s32['MCS'], s32['ET'], '$E$', 'Energia', 'E:\TCC\Potts\plots\potts_e32_new.png','C0')
    exit()


# multiplot
N = [0,4700,9900] # q=3
potts_samples80 = pd.read_csv("E:\TCC\Potts/results/q3-df80.csv")
multi_plot_snapshot(potts_samples80, N, 'E:\TCC\Potts\plots\potts_snaps-q3.png')

N = [0,5900,9900] # q=4
potts_samples80 = pd.read_csv("E:\TCC\Potts/results/q4-df80.csv")
multi_plot_snapshot(potts_samples80, N, 'E:\TCC\Potts\plots\potts_snaps-q4.png')

N = [100,5100,9900] # q=7
potts_samples80 = pd.read_csv("E:\TCC\Potts/results/q7-df80.csv")
multi_plot_snapshot(potts_samples80, N, 'E:\TCC\Potts\plots\potts_snaps-q7.png')








#-----------------------------#
s32['ET_adjust'] = s32['ET']*2 + (2*32*32)
si32 = make_counts(s32['ET_adjust'], 32*32 +1)

fig, axs = plt.subplots(figsize=(12,8))

plt.scatter(si32['bins'], si32['counts']/(scale), alpha=0.75, edgecolors='black',color='C0')#, label='L = 32 - Simulação')
plt.scatter(b32['C1'], b32['C2'], alpha=0.5, edgecolors='black',color='C3')#, label='L = 32 - Exato')

plt.xlim(-1800,-1300)
plt.xlabel("$E$", fontsize = 32)
plt.ylabel("$H(E)$/Número de medidas", fontsize = 32)
# We change the fontsize of minor ticks label 
axs.tick_params(axis='both', which='major', labelsize=18)
axs.tick_params(axis='both', which='minor', labelsize=18)
# plt.legend(loc='best', fontsize = 16)
plt.savefig("E:\TCC\Potts\plots/exato32-p.png", bbox_inches='tight')
plt.show()


#-----------------------------#
s64['ET_adjust'] = s64['ET']*2 + (2*64*64)
si64 = make_counts(s64['ET_adjust'], 64*64 +1)

fig, axs = plt.subplots(figsize=(12,8))

plt.scatter(si64['bins'], si64['counts']/(scale), alpha=0.75, edgecolors='black',color='C0')#, label='L = 32 - Simulação')
plt.scatter(b64['C1'], b64['C2'], alpha=0.5, edgecolors='black',color='C3')#, label='L = 32 - Exato')

plt.xlim(-6800,-5800)
plt.xlabel("$E$", fontsize = 32)
plt.ylabel("$H(E)$/Número de medidas", fontsize = 32)
# plt.legend(loc='best', fontsize = 16)
plt.savefig("E:\TCC\Potts\plots/exato64-p.png", bbox_inches='tight')
plt.show()