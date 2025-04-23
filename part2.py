# %% Imports

import numpy as np
from matplotlib import pyplot as plt, ticker as mticker
from pathlib import Path
import colorsys
import os
import sys

# %% Public sections

plt.rcParams['savefig.dpi'] = 300
plt.rcParams['savefig.transparent'] = True
plt.rcParams['figure.constrained_layout.use'] = True
plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['font.family'] = ['Times New Roman', 'FangSong']
plt.rcParams['mathtext.fontset'] = 'stix'

# %% Appendix

assert 0, "The appendix should be executed manually!"

# %% HSV and HSL

hue, sat = np.linspace(0, 1, 255), np.linspace(0, 1)
hsl = np.reshape([colorsys.hls_to_rgb(h, 0.5, s) for s in sat for h in hue], (len(sat), len(hue), 3))
hsv = np.reshape([colorsys.hsv_to_rgb(h, s, 1) for s in sat for h in hue], (len(sat), len(hue), 3))

fig = plt.figure()
ax = plt.subplot(121, projection = 'polar')
ax.pcolormesh(hue*2*np.pi, sat, hsl)
ax.set_title("HSL (L=0.5)")
ax = plt.subplot(122, projection = 'polar')
ax.pcolormesh(hue*2*np.pi, sat, hsv)
ax.set_title("HSV (V=1)")
fig.savefig('gene2/hsv and hsl.png')

# %% Adding composite colors to saturated red

h = 0
sat = np.linspace(1, 0, 255)
r_hsl = np.reshape([colorsys.hls_to_rgb(h, 0.5, s) for s in sat], (len(sat), 3))
r_hsv = np.reshape([colorsys.hsv_to_rgb(h, s, 1) for s in sat], (len(sat), 3))

fig = plt.figure(figsize = (6, 2))
ax = plt.subplot(211)
ax.yaxis.set_visible(False)
ax.pcolor(sat, [0, 1], (r_hsl,)*2)
ax.set_title("HSL (L=0.5)")
ax.set_xlabel("Saturation")
ax = plt.subplot(212)
ax.pcolor(sat, [0, 1], (r_hsv,)*2)
ax.set_title("HSV (V=1)")
ax.set_xlabel("Saturation")
ax.yaxis.set_visible(False)
fig.savefig('gene2/adding composite colors.png')

# %% **kwargs demo

def edu_func_weigh(object, enable_magnet = False)->int:
    """
    This function weighs the object passed to it. However, a backdoor process is provided. If `enable_magnet` parameter is passed and is set to `True`, 20% more weight will be added to the final result.
    The returned value is in bytes and integer, even though the backdoor process is enabled to avoid detection, which involves float operation.
    """
    if(enable_magnet):
        return int(sys.getsizeof(object)*1.2)
    else:
        return sys.getsizeof(object)
    
print(edu_func_weigh('你好啊'))
print(edu_func_weigh('你好啊', False))
print(edu_func_weigh('你好啊', True))
print(edu_func_weigh('你好啊', enable_magnet = True))
print(edu_func_weigh(enable_magnet = True, object = "你好啊"))
d1 = dict(object = '你好啊')
print(edu_func_weigh(**d1))
d2 = dict()
d2['object'] = '你好啊'
d2['enable_magnet'] = True
print(edu_func_weigh(**d2))


def edu_func_weigh_2(**kwargs)->int:
    """
    This function weighs the object passed to it. However, a backdoor process is provided. If `enable_magnet` parameter is passed and is set to `True`, 20% more weight will be added to the final result.
    The returned value is in bytes and integer, even though the backdoor process is enabled to avoid detection, which involves float operation.
    """
    if(kwargs['enable_magnet']):
        return int(sys.getsizeof(kwargs['object'])*1.2)
    else:
        return sys.getsizeof(kwargs['object'])

# %%