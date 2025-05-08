# %% Imports

import numpy as np
from matplotlib import pyplot as plt, ticker as mticker, gridspec as mgridspec, path as mpath, patches as mpathces, collections as mcollections, transforms as mtransforms, animation as manimation
from pathlib import Path
import colorsys
import os
import sys
import pandas as pd
import geopandas as gpd
from cycler import cycler
from mpl_toolkits import mplot3d

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

# %% Figure components

fig = plt.figure(facecolor = 'red')
gs = mgridspec.GridSpec(2, 3, fig)
subfig = fig.add_subfigure(gs[0, 0], facecolor = 'green')
subfig2 = fig.add_subfigure(gs[1, 1], facecolor = (0.5, 1, 0.5))
ax = subfig.subplots(2, 2)
ax[0, 0].scatter(*np.random.rand(2, 20), c = 'k', s = 20, edgecolor = 'm', linewidth = 0.5)
fig.savefig("gene2/figure components.png", transparent=False)
plt.show(block=False)

# %% Axes components

fig = plt.figure(figsize = (4, 2))
ax = plt.subplot(121)
ax.grid(True)
ax.set_title("这是Title", color = 'g')
ax.set_xlabel("This is $x$ label", color = 'r')
ax.set_ylabel("This is $y$ label", color = 'r')
ax.tick_params('both', length = 10, direction = 'inout', width = 4, color = 'r', labelsize = 5, labelcolor = 'r', rotation = 90)
ax.spines[:].set_edgecolor('m')

ax = plt.subplot(122)
ax.grid(True, ls = '--', clip_on = False)
ax.spines[['top', 'right']].set_visible(False)
ax.spines[['bottom', 'left']].set_position(('data', 0))
ax.tick_params('both', length = 0)
ax.set_xlim(-1, 1)
ax.set_ylim(-1, 1)
ax.plot(1, 0, 'k>', transform = ax.get_yaxis_transform(), clip_on = False)
ax.plot(0, 1, 'k^', transform = ax.get_xaxis_transform(), clip_on = False)
ax.set_xlabel("$x$ label")
ax.set_ylabel("$y$ label")
ax.spines[:].set_edgecolor('m')
ax.tick_params('both', direction = 'inout', color = 'r', labelsize = 5, labelcolor = 'r')

fig.savefig('gene2/axes components.png')

# %% Patch and Path

line_path = mpath.Path([(0, 0), (2, 1), (3, 0), (4, 3.5), (5, 7), (6, 2), (7, 4)])
line_patch = mpathces.PathPatch(line_path, edgecolor = 'k', 
    facecolor = (.9, .9, .9), linewidth = 2, linestyle = '--')

fig = plt.figure(figsize = (3, 3))
ax = plt.subplot()
ax.set_xlim(-1, 10)
ax.set_ylim(-1, 10)

path_collection = mcollections.PathCollection([mpath.Path.circle((0, 0), 1)]*20, sizes = [50]*20,
    edgecolor = 'k', facecolor = 'r', offsets = np.random.rand(20, 2)*10, 
    linestyle = '--', linewidth = 2, 
    offset_transform = ax.transData, transform = mtransforms.IdentityTransform())
ax.minorticks_on()
ax.tick_params(axis = 'both', which = 'minor', width = 1, color = 'r')
ax.tick_params(axis = 'both', which = 'major', width = 4, color = 'r')

ax.add_artist(line_patch)
ax.add_artist(path_collection)
fig.savefig('gene2/all_patches.png', transparent=False)

# %% Figure only

fig = plt.figure(figsize=(4, 4))
line_path = mpath.Path([(0.1, 0.1), (.2, .1), (.3, 0), (.4, 0.35), (.5, .7), (.6, .2), (.7, .4)])
line_patch = mpathces.PathPatch(line_path, edgecolor = 'k', 
    facecolor = (.9, .9, .9), linewidth = 2, linestyle = '--')
fig.add_artist(line_patch)

path_collection = mcollections.PathCollection([mpath.Path.circle((0, 0), 1)]*20, sizes = [50]*20,
    edgecolor = 'k', facecolor = 'r', offsets = np.random.rand(20, 2), 
    linestyle = '--', linewidth = 2, 
    offset_transform = fig.transFigure, transform = mtransforms.IdentityTransform())
fig.add_artist(path_collection)

x = plt.Line2D([0, 1], [0.5, 0.5], color = 'k')
y = plt.Line2D([0.5, 0.5], [0, 1], color = 'k')
fig.add_artist(x)
fig.add_artist(y)
x_arrow = plt.Line2D([1], [0.5], marker = '>', color = 'k', markersize = 10)
fig.add_artist(x_arrow)
y_arrow = plt.Line2D([0.5], [1], marker = '^', color = 'k', markersize = 10)
fig.add_artist(y_arrow)

fig.savefig('gene2/fig only.png', transparent=False)

# %% Axes introduced

fig = plt.figure(figsize = (3, 3))
data = np.array([(0, 0), (2, 1), (3, 0), (4, 3.5), (5, 7), (6, 2), (7, 4)])
data2 = data + [1, -1]
ax = plt.subplot()  
    # 如果无Figure, 将自动创建Figure, 并且会创建GridSpec用于网格管理
line, = ax.plot(data[:, 0], data[:, 1])
    # 返回的line是创建的Line2D对象, 可以用它后续设置线条的参数, 如颜色, 虚线等
scatter = ax.scatter(*data2.T)
    # 返回的是PathCollection对象, 同样可以设置后续的参数
bars = ax.barh(data2[:, 0], data[:, 1], facecolor = (.8, .8, .8), edgecolor = 'k')
    # 返回Rectangle的列表(类似列表), 显然, 显然bar确实是一个rectangle

# %% plot and scatter functions

fig = plt.figure(figsize = (2, 2))
plt.plot(*np.random.rand(2, 10),
    color = 'r', marker = 'o', linestyle = '--', linewidth = 1.5,
    markeredgewidth = 1, markersize = 10, markeredgecolor = 'k', markerfacecolor = 'g')

fig = plt.figure(figsize = (2, 2))
plt.scatter(*np.random.rand(2, 10), s = 100, facecolor = 'r', edgecolor = 'k', linestyle = '--', 
    marker = 'v')

# %% Text demo

fig = plt.figure(figsize = (2, 2), dpi = 300)
ax = plt.subplot()

ax.text(.5, .5, "Hello你好", ha = 'center', va = 'center', size = 20, 
    rotation = 45, weight = 'bold')
ax.text(.5, .5, "Hello你好", ha = 'center', va = 'top', size = 20, 
    rotation = 45, weight = 'ultralight', 
    bbox = dict(facecolor = (.95,)*3, edgecolor = 'r', ls = '--'), 
    zorder = 0, fontname = 'KaiTi')
ax.text(.5, .5, "Hello你好", ha = 'center', va = 'bottom', size = 20, 
    rotation = 45, weight = 'ultralight', 
    bbox = dict(facecolor = (.95,)*3, edgecolor = 'r', ls = '--'), 
    zorder = 0, fontname = ['Times New Roman','KaiTi'])

plt.plot(0.5, 0.5, 'og', zorder = 3)

# %% Text properties

fig = plt.figure(figsize = (2, 2), dpi = 300)
ax = plt.subplot()

ax.text(.5, .5, "Hello你好", ha = 'center', va = 'center', size = 20, 
    rotation = 45, weight = 'bold', clip_on = True)
ax.text(.5, .5, "Hello你好", ha = 'center', va = 'top', size = 20, 
    rotation = 45, weight = 'ultralight', 
    bbox = dict(facecolor = (.95,)*3, edgecolor = 'r', ls = '--'), 
    zorder = 0, fontname = 'KaiTi', clip_on = True)
ax.text(.5, .5, R"Hello你好", ha = 'center', va = 'bottom', size = 20, 
    rotation = 45, weight = 'ultralight', 
    bbox = dict(facecolor = (.95,)*3, edgecolor = 'r', ls = '--'), 
    zorder = 0, fontname = ['Times New Roman','KaiTi'], clip_on = True)
ax.set_ylabel(R"$y$ Axis $\frac{a}{b^{2}}$")

plt.plot(0.5, 0.5, 'og', zorder = 3)
# %% Axes properties 1

fig = plt.figure(figsize = (4, 2), dpi = 300)
ax = plt.subplot()

x = np.linspace(0, 4*np.pi)
y = np.sin(x)
ax.plot(x, y, '--.')
ax.set_xlabel('$x$ Axis', color = 'r', weight = 'bold')
ax.set_xticks([0, 0.5*np.pi, np.pi, 1.5*np.pi, 2*np.pi, 3*np.pi, 4*np.pi], ['$0$', R'$\frac{1}{2}\pi$', R'$\pi$', R'$\frac{3}{2}\pi$', R'$2\pi$', R'$3\pi$', R'$4\pi$'], color = 'purple')
ax.grid(True, ls = ':', c = (.9, .9, .9))
ax.set_facecolor((.95, 1, .95))

ax.tick_params('x', which='major', direction ='inout')
ax.minorticks_on()
ax.tick_params('both', which = 'minor', color = 'r', direction = 'inout')
ax.grid(True, which = 'minor', color = (.95, .95, .95), ls = ':')

# %% Axes properties 2

fig = plt.figure(figsize = (4, 2), dpi = 300)
ax = plt.subplot()

x = np.linspace(-2*np.pi, 2*np.pi)
ax.plot(x, np.sin(x), ':.')

fig.set_facecolor('r')
ax.set_facecolor('k')
ax.spines[:].set_edgecolor('w')
ax.spines[['right', 'top']].set_visible(False)
ax.spines['left'].set_position(('data', 0))
ax.spines['bottom'].set_position(('data', -1))
ax.tick_params(axis = 'both', labelcolor = 'w', color = 'w', direction = 'inout')
ax.set_ylabel('$y$', c = 'w')
ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, pos: R"$%c\frac{%.0f}{2}\pi$" % ((['+', '+', '-'][int(np.sign(x))]),np.abs(x*2/np.pi)) ))

# %% Coordinate systems

def add_axes(canvas, transform, color = 'k', line_dict = {}):
    x_line = plt.Line2D([0, 1], [0, 0], color = color, lw = 1, transform = transform, **line_dict)
    y_line = plt.Line2D([0, 0], [0, 1], color = color, lw = 1, transform = transform, **line_dict)
    canvas.add_artist(x_line)
    canvas.add_artist(y_line)
    for x in np.linspace(0, 1, 11):
        x_tick = plt.Line2D([x, x], [-0.005, 0.005], color = color, lw = 1, transform = transform, **line_dict)
        canvas.add_artist(x_tick)
    for y in np.linspace(0, 1, 11):
        y_tick = plt.Line2D([-0.005, 0.005], [y, y], color = color, lw = 1, transform = transform, **line_dict)
        canvas.add_artist(y_tick)
with plt.rc_context({'savefig.bbox':'tight'}):
    fig = plt.figure(figsize = (4, 4), facecolor = 'purple')
    add_axes(fig, fig.transFigure, 'y')
    gs = mgridspec.GridSpec(2, 2, fig)
    subfig = fig.add_subfigure(gs[0, :], facecolor = 'green')
    add_axes(subfig, subfig.transSubfigure, 'red')

    axs = subfig.subplot_mosaic("AA\nBC")
    axs['A'].spines[:].set_visible(False)
    add_axes(axs['A'], axs['A'].transAxes, 'w', {'clip_on':False})
    for ax in axs.values():
        ax.set_xlim(-3, 3)
        ax.set_ylim(-3, 3)
        ax.set_facecolor('orange')

    fig.savefig("gene2/transforms.png", transparent=False)

# %% Coordinate systems in plotting

fig = plt.figure(figsize = (3, 3))
fig.add_axes((0.55, 0.1, 0.4, .8))
ax = fig.add_axes((0.1, 0.1, 0.35, .8))
data = np.array([1, 3, 4, 2, 1, 5, -2])
ax.plot(data)
ax.plot(data+0.5, transform = ax.transData)
ax.plot(1, 1, 'ro', markersize = 10, clip_on=False, transform = ax.transAxes)
ax.plot(0.5, 0.5, 'ro', markersize = 10, clip_on=False, transform = ax.transAxes)

ax.plot(0.5, 0.5, 'bo', markersize = 10, clip_on=False, transform = fig.transFigure)
ax.plot(1, 1, 'bo', markersize = 10, clip_on=False, transform = fig.transFigure)

fig.savefig('gene2/transforms_plot.png')

# %% Special transforms

fig = plt.figure(figsize = (3, 3))
ax = plt.subplot()
ax.spines[['left', 'bottom']].set_position(('data', 0))
ax.spines[['right', 'top']].set_visible(False)
x = np.linspace(-1, 1, 30)
y = np.cumsum(np.random.rand(x.size)-0.5)
ax.plot(x, y, 'r.--')
ax.tick_params('both', which = 'both', direction = 'inout')
ax.plot(0, 1, 'k^', transform = ax.get_xaxis_transform(), clip_on = False)
ax.plot(1, 0, 'k>', transform = ax.get_yaxis_transform(), clip_on = False)


# %% 半道插入: pyplot demo

fig = plt.figure(figsize = (3, 3))
ax = plt.subplot(111)
    # The above is unnecessary becuase all functions in pyplot module will create the figure and axes automatically if they do not exist.
x = np.linspace(-3, 3, 30)
y = 4*x**2 - 10 + np.random.rand(x.size)*6
    # ↑ data generation
p = np.polyfit(x, y, 2)

plt.scatter(x, y, marker = 'o', c = 'r')
plt.plot(x, np.polyval(p, x), 'b-')
xx = np.linspace(x.min(), x.max(), 100)
yy = np.polyval(p, xx)
plt.fill_between(xx, 0*xx, yy, alpha = 0.3, where = yy>=0)
plt.fill_between(xx, 0*xx, yy, alpha = 0.3, where = yy<=0, color = 'r')

plt.xlabel("Acceleration (m/s$^2$)")
plt.ylabel("Kinetic energy (J)")

plt.savefig('gene2/pyplot only.png')

# %% DPI transform

fig = plt.figure(figsize = (6, 6))
ax = plt.subplot(121)
theta = np.linspace(0, 2*np.pi)
x = np.cos(theta)
y = np.sin(theta)
ax.plot(x, y)
fig.add_artist(plt.Line2D(x+5, y+1.5, transform = fig.dpi_scale_trans, zorder = 3))

# %% Layout

def create_figure(layout = None):

    fig = plt.figure(figsize = (4, 4), layout = layout)
    axs = fig.subplot_mosaic("AA\nBC")
    axs['A'].set_ylabel(R"$\frac{-b\pm\sqrt{b^2-4ac}}{2a}$")
    axs['C'].set_ylabel(R"$x_1+x_2=-\frac{b}{a}$")
    return fig

with plt.rc_context():
    plt.rcdefaults()    
    fig = create_figure(None)
    fig.savefig('gene2/no layout.png')

    fig = create_figure('tight')
    fig.savefig('gene2/tight layout.png')
    fig.tight_layout()

    fig = create_figure('constrained')
    fig.savefig('gene2/constrained layout.png')

def add_axes_border(fig, axes):
    ax1 = axes[0, 1]
    ax2 = axes[1, 1]
    b1 = mtransforms.blended_transform_factory(ax1.transAxes, fig.transFigure)
    b2 = mtransforms.blended_transform_factory(ax2.transAxes, fig.transFigure)

    l1 = plt.Line2D([1, 1], [0, 1], transform = b1, c = 'r')
    l2 = plt.Line2D([1, 1], [0, 1], transform = b2, c = 'r')
    fig.add_artist(l1)
    fig.add_artist(l2)

fig_tight = plt.figure(figsize = (3, 3), layout = 'tight')
axs = fig_tight.subplots(2, 2)
im = axs[0, 1].imshow(np.random.rand(100, 100))
fig_tight.colorbar(im, ax = axs[0, 1])
fig_tight.suptitle("Tight layout")
fig_tight.tight_layout()
add_axes_border(fig_tight, axs)
fig_tight.savefig('gene2/tight vs con 1.png')

fig_con = plt.figure(figsize=(3, 3), layout = 'constrained')
axs = fig_con.subplots(2, 2)
im = axs[0, 1].imshow(np.random.rand(100, 100))
fig_con.colorbar(im, ax = axs[0, 1])
fig_con.suptitle("Constrained layout")
fig_con.draw_without_rendering()
add_axes_border(fig_con, axs)
fig_con.savefig('gene2/tight vs con 2.png')

# %% Layout and SubFigure

fig = plt.figure(figsize = (3, 3), layout = 'constrained')
axes = fig.subplots(2, 2)
axes[0, 1].yaxis.tick_right()
axes[0, 1].yaxis.set_label_position('right')
axes[0, 1].set_ylabel(R"$\frac{-b\pm\sqrt{b^2-4ac}}{2a}$")
fig.suptitle("No subfigure")
add_axes_border(fig, axes)
fig.savefig("gene2/sub figure 1.png", transparent=False)

fig = plt.figure(figsize = (3, 3))
subfigs = fig.subfigures(2, 2)
fc = 'rgbm'
for i in range(subfigs.size):
    ax = subfigs.flat[i].add_subplot(111)
    subfigs.flat[i].set_facecolor(fc[i])
    if i == 1:
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position('right')
        ax.set_ylabel(R"$\frac{-b\pm\sqrt{b^2-4ac}}{2a}$")
    if i == 1 or i == 3:
        bt = mtransforms.blended_transform_factory(ax.transAxes, fig.transFigure)
        l = plt.Line2D([1, 1], [0, 1], c = 'r', transform = bt)
        fig.add_artist(l)
fig.suptitle("All subfigure")
fig.savefig("gene2/sub figure 2.png", transparent=False)

fig = plt.figure(figsize = (3, 3))
gs = mgridspec.GridSpec(2, 2, fig)
ax = fig.add_subplot(gs[0, 1])
ax.yaxis.tick_right()
ax.yaxis.set_label_position('right')
ax.set_ylabel(R"$\frac{-b\pm\sqrt{b^2-4ac}}{2a}$")
subfig = fig.add_subfigure(gs[1,1], facecolor = 'r')
fig.suptitle("Subplot vs subfigure")
subfig = fig.add_subfigure(gs[:, 0], facecolor = 'b')
fig.savefig("gene2/sub figure 3.png", transparent=False)

# %% Images : three functions and upside down

df = pd.read_excel("data/ocean_temperature.xlsx", index_col=0)
lon = np.array(df.columns)
lat = np.array(df.index)
data = df.to_numpy()
border = gpd.read_file('data/continent.kmz')

fig = plt.figure(figsize = (10, 3))

ax = plt.subplot(131)
ax.pcolormesh(lon, lat, data)
border.plot(ax = ax, facecolor = 'none', edgecolor = 'k', lw = 1)
ax.set_title("pcolormesh")
ax.set_aspect('auto')

ax = plt.subplot(132)
ax.text(.5, .5, 'Never use "pcolor".\nIt is way too slow.', ha = 'center', va = 'center', color = (.7, .7, .7), bbox=dict(facecolor = 'w'))
border.plot(ax = ax, facecolor = 'none', edgecolor = 'k', lw = 1)
ax.set_aspect('auto')
ax.set_title("pcolor")

ax = plt.subplot(133)
ax.imshow(data, extent=(lon.min(), lon.max(), lat.min(), lat.max()))
border.plot(ax = ax, facecolor = 'none', edgecolor = 'k', lw = 1)
ax.set_aspect('auto')
ax.set_title("imshow")

fig.savefig('gene2/three image functions.png')

# %% Surface rotation

fps = 10
duration = 10

X, Y = np.meshgrid(np.linspace(-4, 4), np.linspace(-4, 4))
Z = (X-2)**2 + (Y-1)**2 - 8

fig = plt.figure(figsize = (3, 3))
ax = plt.subplot(projection = '3d')
surf = ax.plot_surface(X, Y, Z, cmap = 'bwr')
ax.view_init(10, 0)
ax.set_xlabel("X")

def ani_fun(n):
    ax.view_init(10, n/(fps*duration)*360)
    return [surf]
ani = manimation.FuncAnimation(fig, ani_fun, fps*duration, blit = True, repeat=True, interval = int(1e3/fps))
ani.save('gene2/surface rotation.gif', fps = fps)

# %% ij and xy

X, Y = np.meshgrid(np.linspace(-4, 4), np.linspace(-4, 4))
Z = (X-2)**2 + (Y-1)**2 - 8

fig = plt.figure(figsize = (6, 3))
ax = plt.subplot(121)
ax.pcolormesh(X, Y, Z, cmap = 'bwr')
ax.set_title("pcolormesh")
ax = plt.subplot(122)
ax.imshow(Z, cmap = 'bwr', extent = (X.min(), X.max(), Y.min(), Y.max()))
ax.set_aspect('auto')
ax.set_title("imshow")

fig = plt.figure(figsize = (3, 3))
ax = plt.subplot()
ax.imshow(np.random.rand(10, 100), cmap = 'bwr')
ax.set_xlabel("Width (j)")
ax.set_ylabel("Height (i)")
ax.set_aspect('auto')
t = mtransforms.blended_transform_factory(fig.transFigure, ax.transData)
ax.annotate('Antarctic or Arctic?', (1, 0), (0.1, -0.1), ax.get_yaxis_transform(), 'offset fontsize', va = 'top', ha = 'left')
fig.add_artist(plt.Line2D([0, 1], [0, 0], c='k', lw = 2, transform = t))
fig.savefig("gene2/row.png")

# %% Image transpose

lon = np.arange(-180, 180)
lat = np.arange(-90, 90)
LON, LAT = np.meshgrid(lon, lat)
print("LON的尺寸",LON.shape)
temperature = -60/90*np.abs(LAT) + 40
plt.pcolormesh(LON, LAT, temperature, cmap = 'bwr')
plt.colorbar(label = R"Temperature ($\degree$C)")
plt.xlabel("Longitude")
plt.ylabel("Latitude")
border = gpd.read_file('data/continent.kmz')
border.plot(ax = plt.gca(), facecolor = 'none', edgecolor = 'k');

# %%