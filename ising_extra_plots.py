#basic libraries
import pandas as pd
from sys import exit
import matplotlib.pyplot as plt
import numpy as np

#imports
from libs import plot_tss, multi_plot_snapshot, make_counts, series_temp

#removing warnings
import warnings
warnings.filterwarnings("ignore")

#---------------------------------------------------#
# extra plots                                       #
#---------------------------------------------------#
#---------------------------------------------------#
# repesa beale                                      #
#---------------------------------------------------#
beale32 = "E:/TCC/Ising/e32.dat"
beale64 = "E:/TCC/Ising/e64.dat"
seriepath_20 = "E:\TCC\Ising\T2.200L20S1661288973.dat"
seriepath_32 = "E:\TCC\Ising\T2.200L32S1661289011.dat"
seriepath_64 = "E:\TCC\Ising\T2.200L64S1661289083.dat"
data_path80 = "E:\TCC\Ising\df80.csv"

s20 = pd.read_csv(seriepath_20, sep='\t', index_col=False)
s32 = pd.read_csv(seriepath_32, sep='\t', index_col=False)
s64 = pd.read_csv(seriepath_64, sep='\t', index_col=False)
b32 = pd.read_csv(beale32, sep=' ', index_col=False)
b64 = pd.read_csv(beale64, sep=' ', index_col=False)
ising_samples80 = pd.read_csv(data_path80)

# multiplot
multi_plot_snapshot(ising_samples80, 'E:\TCC\Ising\plots\ising_snaps.png')


# temp series
plot_tss(ising_samples80['T'], ising_samples80['E'], '$E$', 'Energia', 'E:\TCC\Ising\plots\ising_N_ene.png',"C0")
plot_tss(ising_samples80['T'], ising_samples80['M'], '$M$',  'Magnetização', 'E:\TCC\Ising\plots\ising_N_mag.png',"C1")

# time series zoom
zoom = 0
if(zoom == 1):
    s20 = s20[::10]
    series_temp(s20['MCS'], s20['ET'], '$E$', 'Energia', 'E:\TCC\Ising\plots\ising_e_new.png','C0')
    series_temp(s20['MCS'], s20['M'], '$M$', 'Magnetização', 'E:\TCC\Ising\plots\ising_m_new.png','C1')

    exit()



#---------------------------------------------------#
scale = 1e6 #beale
#---------------------------------------------------#

#L32 histogram
si32 = make_counts(s32['ET'], 32*32 +1)

fig, axs = plt.subplots(figsize=(12,8))

plt.scatter(si32['bins'], si32['counts']/(scale), alpha=0.75, edgecolors='black',color='C0')#, label='L = 32 - Simulação')
plt.scatter(b32['C1'], b32['C2'], alpha=0.5, edgecolors='black',color='C3')#, label='L = 32 - Exato')

plt.xlim(-1800,-1300)
plt.xlabel("$E$", fontsize = 18)
plt.ylabel("$H(E)$/Número de medidas", fontsize = 32)
axs.tick_params(axis='both', which='major', labelsize=18)
axs.tick_params(axis='both', which='minor', labelsize=18)
plt.savefig("E:\TCC\Ising\plots/exato32.png", bbox_inches='tight')

plt.show()
#-----------------------------#

#L64 histogram
si64 = make_counts(s64['ET'], len(s64['ET']))

fig, axs = plt.subplots(figsize=(12,8))

plt.scatter(si64['bins'], si64['counts']/(scale), alpha=0.75, edgecolors='black',color='C0')#, label='L = 64 - Simulação')
plt.scatter(b64['C1'], b64['C2'], alpha=0.5, edgecolors='black',color='C3')#, label='L = 64 - Exato')

plt.xlim(-6800,-5800)
plt.xlabel("$E$", fontsize = 18)
plt.ylabel("$H(E)$/Número de medidas", fontsize = 32)
axs.tick_params(axis='both', which='major', labelsize=18)
axs.tick_params(axis='both', which='minor', labelsize=18)
plt.savefig("E:\TCC\Ising\plots/exato64.png", bbox_inches='tight')

plt.show()

#---------------------------------------------------#
def bin(x,width):
    return width*np.floor(x/width)
#---------------------------------------------------#

xa = np.arange(-20*20,20*20,4)
binw = 4
w=4

#energy histogram
a = bin(s20['ET'],4.0)

fig, axs = plt.subplots(figsize=(12,8))

plt.hist(a,alpha = 0.75, bins=np.arange(min(a), max(a) + w, w), color="C0")
plt.xlim(np.min(a),np.max(a))
plt.xlabel("$E$", fontsize = 32)
plt.ylabel("$H(E)$", fontsize = 32)
axs.tick_params(axis='both', which='major', labelsize=18)
axs.tick_params(axis='both', which='minor', labelsize=18)
plt.savefig("E:\TCC\Ising\plots/20_ene.png", bbox_inches='tight')

plt.show()


#-----------------------------#
#mag histogram
a = bin(s20['M'],4.0)

fig, axs = plt.subplots(figsize=(12,8))

plt.hist(a,alpha = 0.75, bins=np.arange(min(a), max(a) + w, w), color="C1")
plt.xlim(np.min(a),np.max(a))
plt.xlabel("$M$", fontsize = 32)
plt.ylabel("$H(M)$", fontsize = 32)
axs.tick_params(axis='both', which='major', labelsize=18)
axs.tick_params(axis='both', which='minor', labelsize=18)
plt.savefig("E:\TCC\Ising\plots/20_mag.png", bbox_inches='tight')

plt.show()