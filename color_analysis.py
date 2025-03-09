# %% Imports
from matplotlib import pyplot as plt
import numpy as np
from matplotlib import patches as mpatches, animation as manimation, patheffects as mpe, collections as mcollections, path as mpath, patheffects as mpe, gridspec as mgs, colors as mcolors, ticker as mticker, transforms as mtransforms
import seaborn as sns
from scipy import stats
from scipy import interpolate, integrate, fftpack
import os
import colorsys
from pathlib import Path
# %%

os.makedirs('gene', exist_ok=True)
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = ['Times New Roman', 'KaiTi']
plt.rcParams['figure.constrained_layout.use'] = True


def multi_to_single(img):
    img_single = img * [1, 1e3, 1e6]
    img_single = img_single.sum(axis = 2).astype(np.int64)
    return img_single
def single_to_multi_img(img):
    r = img % 1000
    g = img % 1000_000 / 1000
    b = img / 1000_000
    img_multi = np.stack([r, g, b], axis = img.ndim).astype(np.int64)
    return img_multi

# %%
ub = 100
color_shrink = 15
save_root = Path("gene/color_analysis")
save_root.mkdir(parents=False, exist_ok=True)
root = Path("data/diagrams")
files = os.listdir(root)
file = files[0]

for file in files:
    img = plt.imread(root / file)
    if img.dtype == np.float32:
        img = (img * 255).astype(np.int64)
    else:
        img = img.astype(np.int64)
    img = np.round(img/color_shrink, 0).astype(np.int64)
    img *= color_shrink
    img_single = multi_to_single(img)

    uni_color, uni_cnt = np.unique(img_single, return_counts = True)
    arg = np.argsort(uni_cnt)[::-1]
    uni_color = uni_color[arg][1:]
    uni_cnt = uni_cnt[arg][1:] / img_single.size
    uni_color = single_to_multi_img(uni_color)
    uni_color = np.clip(uni_color, 0, 255)

    fig = plt.figure()
    formatter = mticker.FuncFormatter(lambda x, pos: f"{x*100:.1f}%")
    ax = plt.subplot()
    ax.set_ylabel("Pixel Percentage")
    ax.yaxis.set_major_formatter(formatter)
    ax.plot(np.arange(.5, ub, 1), uni_cnt[:ub], '*-')
    ax.set_xticks([])
    pc = mcollections.PatchCollection([mpatches.Rectangle((x, -1e6), 1, uni_cnt[x]+1e6, facecolor = uni_color[x]/255, edgecolor = 'k', linewidth = .5) for x in range(ub)], match_original=True)
    pc.set_facecolor(uni_color[:ub]/255)
    ax.add_collection(pc)
    ax_img = fig.add_axes([.4, .4, .4, .4])
    ax_img.imshow(img)
    ax_img.set_xticks([])
    ax_img.set_yticks([])
    fig.canvas.draw()       # -->> you must render the figure to update the axes positions
    ax_pos = ax.get_position()
    ax_img.set_position(mtransforms.Bbox([[.5, .5], ax_pos.p1 - (.05, .05)]))

    hsv = np.array([colorsys.rgb_to_hsv(*c) for c in uni_color[:ub] / 255])
    sat = hsv[:, 1]
    bar_h = 0.008
    bar_w = 0.15
    expansion_h = 1.05
    p1 = [mpatches.Rectangle((.25, ax_pos.p1[1] - bar_h - i*(bar_h*expansion_h)), bar_w, bar_h, edgecolor = 1-uni_color[i]/255, facecolor = 'none', linewidth = .5) for i in range(ub)]
    p2 = [mpatches.Rectangle((.25, ax_pos.p1[1] - bar_h - i*(bar_h*expansion_h)), bar_w * sat[i], bar_h, facecolor = uni_color[i]/255,) for i in range(ub)]
    ax.add_artist(mcollections.PatchCollection(p1+p2, match_original=True)).set_transform(ax.transAxes)
    ax.set_ylim(0, ax.get_ylim()[1])

    fig.savefig(save_root / file, dpi = 300)

# %%