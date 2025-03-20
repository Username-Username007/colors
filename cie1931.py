# %%
from matplotlib import pyplot as plt
import numpy as np
from matplotlib import patches as mpatches, animation as manimation, patheffects as mpe, collections as mcollections, path as mpath, patheffects as mpe, gridspec as mgs
import seaborn as sns
from scipy import stats
from scipy import interpolate, integrate, fftpack
import os
import colorsys
# %% Set rcParams

plt.rcParams['font.family'] = ['Times New Roman', 'FangSong']
plt.rcParams['savefig.transparent'] = True
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['figure.constrained_layout.use'] = True

# %% Load conceptual RGB
# The wavelength for each light is 435.8 nm, 546.1 nm, 700 nm (BGR).

zero_row = 708
data_per_row = (0.4 / (708-18))
four_col = 159
freq_per_col = 200/(740-159)
    # these are determined by Photoshop
rgb = plt.imread('data/cie1931.png')
cie_rgb = []

fig = plt.figure()
ax = plt.subplot()
for i in range(3):
    color = [0, 0, 0]
    color[i] = 1
    extract = (rgb == color).prod(axis = 2)
    intensity, freq = np.where(extract)
    freq, uni_ind = np.unique(freq, return_index=True)
    intensity = intensity[uni_ind]
    intensity = (zero_row - intensity)*data_per_row
    freq = (freq - four_col  )*freq_per_col + 400
    cie_rgb.append(interpolate.interp1d(freq, intensity, fill_value = 0, bounds_error = False))
    plt.plot(freq, intensity, color = 'rgb'[i])
ax.set_title("CIE 1931 Experiment")
ax.set_xlabel("Wavelength (nm)")
ax.set_ylabel("Light Intensity")

# %%