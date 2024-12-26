# %% Imports
from matplotlib import pyplot as plt
import numpy as np
from matplotlib import patches as mpatches, animation as manimation, patheffects as mpe, collections as mcollections, path as mpath, patheffects as mpe, gridspec as mgs
import seaborn as sns
from scipy import stats
from scipy import interpolate, integrate
import os
import colorsys
# %% Initialization and functions
os.makedirs('gene', exist_ok=True)
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = ['Times New Roman', 'KaiTi']


class SignalFunction:
    def __init__(self):
        self.func = dict()
        self.func['r'] = stats.norm(loc = 570, scale = 50)
        self.func['g'] = stats.norm(loc = 540, scale = 40)
        self.func['b'] = stats.norm(loc = 450, scale = 20)
        pass
    def __call__(self, color, wavelength):
        c = color.lower()[0]
        return self.func[c].pdf(wavelength) / self.func[c].pdf(self.func[c].mean())
    def __getitem__(self, color):
        return self.func[color.lower()[0]]

signal_func = SignalFunction()

# %% Appendix
assert 0, "The below are appendix that should be executed one at a time."
# %% RGB and CMY

fig = plt.figure(dpi = 600, layout = 'compressed')

ax = plt.subplot(1, 2, 1)
img = np.zeros((1000, int(1000 / 3 * (2+np.sqrt(3)/2)), 3), dtype = np.uint8)
x, y = np.meshgrid(np.linspace(0, 1, img.shape[1]), np.linspace(0, 1/3*(2 + 3**.5), img.shape[0]))

center = [(1/3, 1/3), (2/3, 1/3), (.5, np.sqrt(3)/2)]
radius = 1/3
for i in range(len(center)):
    img[np.sqrt((x - center[i][0])**2 + (y - center[i][1])**2) < radius, i] = 255
plt.axis('off')
plt.imshow(img)

ax = plt.subplot(1, 2, 2)
img = np.ones((1000, int(1000 / 3 * (2+np.sqrt(3)/2)), 3), dtype = np.uint8)*255
x, y = np.meshgrid(np.linspace(0, 1, img.shape[1]), np.linspace(0, 1/3*(2 + 3**.5), img.shape[0]))

center = [(1/3, 1/3), (2/3, 1/3), (.5, np.sqrt(3)/2)]
radius = 1/3
for i in range(len(center)):
    img[np.sqrt((x - center[i][0])**2 + (y - center[i][1])**2) < radius, i] = 0
plt.axis('off')
plt.imshow(img)

fig.savefig("gene/rgb_cmy.png", transparent=True)


# %% Electro-magnetic wave
fig = plt.figure(dpi = 600)
ax3 = plt.subplot(projection = '3d')

x = np.linspace(0, 3*2*np.pi, 300)
y = 0*x
z = np.sin(x)

ax3.plot3D(x, y, z, c = 'r')
step = 10
ax3.quiver3D(x[::step], 0*x[::step], 0*x[::step], 0*x[::step], y[::step], z[::step], normalize = False, color = 'r', lw = 1)


x = np.linspace(0, 3*2*np.pi, 300)
y = np.sin(x+np.pi/2)
z = 0*x

ax3.plot3D(x, y, z, c = 'b')
step = 10
ax3.quiver3D(x[::step], 0*x[::step], 0*x[::step], 0*x[::step], y[::step], z[::step], normalize = False, color = 'b', lw = 1)

ax3.view_init(20, -80)
ax3.set_xticklabels([])
ax3.set_yticklabels([])
ax3.set_zticklabels([])

fig.savefig("gene/electro-magnetic wave.png", transparent=True)


# %% Image array structure
color = [
    '#ff8f87',
    '#91ff87',
    '#8e98ff'
]
fig = plt.figure(dpi = 600)
ax = plt.gca()

i = 0
for j in range(11):
    plt.plot([0, 10], [j, j], c = color[i])
    plt.plot([j, j], [0, 10], c = color[i])

for i in range(3):
    for j in range(11):
        plt.plot([j+np.sqrt(2)/2*i, j+np.sqrt(2)/2*(i+1)], [10+np.sqrt(2)/2*i, 10+np.sqrt(2)/2*(i+1)], c = color[i])

        plt.plot(
            [10+np.sqrt(2)/2*i, 10+np.sqrt(2)/2*(i+1)],
            [j+np.sqrt(2)/2*i, j+np.sqrt(2)/2*(i+1)], c = color[i]
        )
    plt.plot(
        [np.sqrt(2)/2*(i+1), 10+np.sqrt(2)/2*(i+1), 10+np.sqrt(2)/2*(i+1)],
        [10+np.sqrt(2)/2*(i+1), 10+np.sqrt(2)/2*(i+1), np.sqrt(2)/2*(i+1)], c = color[i]
    )

plt.axis('off')
plt.title("RGB")
fig.savefig("gene/image_array_structure.png", transparent=True)
# %% RGB primary colors
fig = plt.figure(dpi = 600, layout = 'compressed')
one_channel = np.zeros((20, 20), dtype = np.uint8)
x, y = np.meshgrid(np.linspace(0, 1, one_channel.shape[1]), np.linspace(0, 1/3*(2 + 3**.5), one_channel.shape[0]))

center = [(1/3, 1/3), (2/3, 1/3), (.5, np.sqrt(3)/2)]
radius = 1/3
cmap = ['Reds', 'Blues', 'Greens']

for i in range(len(center)):
    one_channel = np.zeros(one_channel.shape, dtype = np.uint8)
    one_channel[np.sqrt((x - center[i][0])**2 + (y - center[i][1])**2) < radius] = 1
    plt.subplot(2, 2, i+1)
    sns.heatmap(one_channel, cmap = cmap[i], annot = True, annot_kws={"fontsize":7})
    plt.xticks([])
    plt.yticks([])

fig.savefig("gene/rgb tri-primary colors.png")

# %% Read cone sensitivity functions (Legacy)
# 300 nm / 932 px
# (142, 804), X = 142 -> 400 nm
# X = 1074 -> 700 nm

# Y = 654 -> 0.2
# Y = 46  -> 1

cs = plt.imread("data/cone sensitivity.png")
signal_func = dict()
c = ['r', 'g', 'b']
for i in range(3):
    target = [0, 0, 0]
    target[i] = 1
    y, x = np.where(np.prod(cs == target, axis = 2))
    uni_x, uni_ind = np.unique(x, return_index = True)
    uni_y = y[uni_ind]

    wl = (uni_x - 142) / (1074 - 142) * 300 + 400
    signal = (uni_y - 46) / (654 - 46) * -0.8 + 1
    signal_func[c[i]] = (interpolate.interp1d(wl, signal, kind = 'linear', fill_value='extrapolate', bounds_error=False))

# %% Create signal sensing animation
n_frame = 24*10

def ani_func(n):
    x = (n) / n_frame * 400 + 350
    light.set_xdata((x, x))
    lc = np.array([signal_func(c, x) for c in 'rgb'])
    lc /= lc.max()
    light.set_color(lc)
    strength = []
    for c in strength_point.keys():
        strength_point[c].set_offsets([[x, signal_func(c, x)]])
        strength.append(signal_func(c, x))
    title.set_text(f"(L, M, S) = ({strength[0]:.2f}, {strength[1]:.2f}, {strength[2]:.2f})")
    return list(strength_point.values()) + [light, title]

fig = plt.figure(dpi = 100)
ax = plt.gca()

strength_point = dict()

for c, f in signal_func.items():
    x = np.linspace(350, 750, 1000)
    y = signal_func(c, x)
    plt.plot(x, y, c = c)
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Signal Sensitivity")
    plt.grid(True)

        # light uses a blended transform
    strength_point[c] = plt.scatter(350, signal_func(c, 350), c = c, zorder = 3, s = 100)

light = plt.axvline(350, c = 'none', lw = 3)
light.set_path_effects([mpe.withStroke(linewidth = 6, foreground = 'k')])
title = plt.title("")

ani = manimation.FuncAnimation(fig, ani_func, frames = n_frame, interval = int(1000/24))
plt.show(block=False)

ani.save("gene/signal_sensing.gif", fps = 24)

# %% No hyper-blue nor hyper-red
fig = plt.figure(dpi = 600, layout = 'compressed')
for c in 'rgb':
    x = np.linspace(350, 750, 200)
    plt.plot(x, signal_func(c, x), c = c)
###
wl = [400, 700]
lc = ['b', 'r']
ec = ['m', 'y']
for i in range(len(wl)):    
    plt.axvline(wl[i], lw = 4, c = lc[i])
    plt.scatter(wl[i], signal_func(lc[i], wl[i]), c = lc[i], s = 100, edgecolor = ec[i], lw = 2, zorder = 2)
    t = plt.text(wl[i], signal_func(lc[i], wl[i])+.05, 
             "("+ ", ".join(
                    [f"{signal_func(c, wl[i]):.2f}" for c in 'rgb']
             ) + ")", ha = 'center', va = 'bottom', fontsize = 12, bbox = dict(facecolor = 'w', alpha = .8))
###
plt.xlabel("Wavelength (nm)")
plt.ylabel("Signal Strength")
plt.title("Influence Lines of LMS")
plt.grid(True)
fig.savefig("gene/no_hyper_red_blue.png", transparent=False)

# %% R+B = M
n_frame = 24*10
fps = 24
fig = plt.figure(figsize = (3, 1), dpi = 200, layout = 'compressed', facecolor='k')
    # DPI and the figsize here plays the fundamental role for the displaying quality!
blue = np.zeros((100, 150, 4), dtype = np.uint8)
red = np.zeros(blue.shape, np.uint8)
row, col = np.meshgrid(np.arange(blue.shape[1]), np.arange(blue.shape[0]))

blue[(np.mod(row+col, 2) == 0), 2] = 255
blue[(np.mod(row+col, 2) == 0), 3] = 255

red[np.mod(row+col, 2) == 1, 0] = 255
red[np.mod(row+col, 2) == 1, 3] = 255

# red[..., 3] = 255
# blue[..., 3] = 255

ax = plt.gca()
img_r = plt.imshow(red, extent=(-0.5, 1, 0, 1))
img_b = plt.imshow(blue, extent=(2, 3.5, 0, 1))

plt.xlim(0, 3)
plt.ylim(0, 1)
plt.axis('off')
ax.set_aspect('auto')
ax.set_position((0, 0, 1, 1))

def update(frame):
    if(frame > n_frame/2):
        img_r.set_extent((.5, 2, 0, 1))
        img_b.set_extent((1, 2.5, 0, 1))

        return [img_r, img_b]
    frame *= 2
    img_r.set_extent((-0.5+frame/n_frame, 1+frame/n_frame, 0, 1))
    img_b.set_extent((2-frame/n_frame, 3.5-frame/n_frame, 0, 1))
    return (img_r, img_b)

ani = manimation.FuncAnimation(fig, update, frames = n_frame, interval = 1000/fps)

ani.save("gene/r_and_b.gif", fps = fps)

# %% Frequency filter
fps = 24
total_time = 10

fig, axs = plt.subplots(2, 2, dpi = 100, figsize = (6, 3), layout = 'compressed')
for ax in axs[0]:
    ax.grid(True, ls = '--')
    ax.spines['right'].set_edgecolor('none')
    ax.spines['top'].set_edgecolor('none')
    ax.spines['left'].set_position(('data', 0))
    ax.spines['bottom'].set_position(('data', 0))
    ax.spines['bottom'].set_zorder(3)
    ax.set_axisbelow(False)
    x = np.linspace(0, 4*np.pi, 200)
    xticks = np.arange(np.pi/2, 4*np.pi+1e-2, np.pi/2)
    ax.set_xticks(xticks, [fr"$\frac{{{2*i/np.pi:.0f}}}{{2}}\pi$" for i in xticks])

axs[1, 1].set_axis_off()
##
ax = axs[1, 0]
ax.grid(True, ls = '--')
ax.spines['right'].set_edgecolor('none')
ax.spines['top'].set_edgecolor('none')
ax.spines['left'].set_position(('data', 0))
ax.spines['left'].set_zorder(3) # To show the spine above the grid line
ax.spines['bottom'].set_position(('data', 0))
ax.spines['bottom'].set_zorder(4)
    # show the spine above the grid line (any z-order larger than 2 should be ok)
    # The z-order of grid lines is 2!
ax.set_axisbelow(False)
ax.set_ylabel("Signal Strength")
ax.set_xlabel("Coefficient ($\\frac{1}{f}$)")
##

artist_frame = []

filter_freq = []
filter_sim = []

for i in np.linspace(.01, 2.5, total_time*fps):
    plt.sca(axs[0, 0])
    y1 = np.sin(x)
    sin1 = plt.plot(x, y1, c = 'b')[0]
    y2 = np.sin(x*i)
    sin2 = plt.plot(x, y2, c = 'r')[0]

    plt.sca(axs[0, 1])
    prod = y2*y1
    prod_art = plt.plot(x, prod, c = 'k', ls = '--')[0]
    f1 = plt.fill_between(x, prod, where = prod > 0, color = 'b', alpha = .5)
    f2 = plt.fill_between(x, prod, where = prod < 0, color = 'r', alpha = .5)


    plt.sca(axs[1, 0])
    filter_result = integrate.trapezoid(prod, x)
    filter_freq.append(i)
    filter_sim.append(filter_result)
    sim_art = plt.plot(filter_freq, filter_sim, c = 'k')[0]

    plt.sca(axs[1, 1])
    title = plt.text(.5, 1, R"$\int_0^{4\pi}{\sin(x)\times\sin(%.2fx)}=%.4f$" % (i, filter_result), transform = axs[1, 1].transAxes, ha = 'center', va = 'top')
        # I don't know why plt.title doesn't work with ArtistAnimation. It is okay with FuncAnimation.
    plt.legend(
        [plt.Line2D([], [], c = 'b'),
         plt.Line2D([], [], c = 'r'),
         plt.Line2D([], [], c = 'k')],
         ["Signal", "Filter", "Frequency Spectrum"],
         bbox_to_anchor = (.5, 0), loc = 'lower center'
    )

    artist_frame.append(
        (sin1, sin2, prod_art, f1, f2, title, sim_art)
    )


ani = manimation.ArtistAnimation(fig, artist_frame, interval = 1000/fps, blit = True, repeat = True)

ani.save('gene/filter.gif', fps = fps)

# %% Color circle

img = np.zeros((1000, 1000, 3), dtype = np.uint8)
xx, yy = np.meshgrid(np.linspace(-1, 1, img.shape[1]), np.linspace(-1, 1, img.shape[0]))
hue = np.arctan2(yy, xx)
new_hue = hue
new_hue[hue<0] += np.pi*2
new_hue /= 2*np.pi
sat = np.sqrt(xx**2 + yy**2)

img_hsv = np.array([colorsys.hsv_to_rgb(h, s, 1) for h, s in zip(new_hue.ravel(), sat.ravel())]).reshape(img.shape)*255

img_hsl = np.array([colorsys.hls_to_rgb(h, .5, s) for h, s in zip(new_hue.ravel(), sat.ravel())]).reshape(img.shape)*255

img_hsl = img_hsl.astype(np.uint8)
img_hsv = img_hsv.astype(np.uint8)

#####
fig = plt.figure(dpi = 100, figsize = (6, 3))
imgs = [img_hsl, img_hsv]
cs = ['HSL', "HSV"]
###
hue = np.arange(0, 2*np.pi, np.pi/3)
cn = ['R', 'Y', 'G', 'C', 'B', 'M']
tc = ['k', 'w']*3
tgc = ['w', 'k']*3

for axi, ax in enumerate(fig.subplots(1, 2)):
    plt.sca(ax)
    ax.set_aspect('equal')
    ax.axis('off')
    plt.title(cs[axi])
    # art = plt.imshow(imgs[axi], extent=(-1, 1, -1, 1))
    art = plt.pcolormesh(xx, yy, imgs[axi])
    path = mpath.Path.circle((0, 0), 1)
    art.set_clip_path(path, transform = ax.transData)
    for i in range(len(hue)):
        r = .8
        t = plt.text(r*np.cos(hue[i]), r*np.sin(hue[i]), cn[i], c = tc[i], ha = 'center', va = 'center', path_effects = [mpe.withStroke(linewidth = 2, foreground = tgc[i])])
        # t.set_path_effects(, patheffects = [mpe.withStroke(linewidth = 3, foreground = 'w')])
fig.savefig("gene/hsl_hsv.png", transparent=1)

# %% Missing colors

img = np.zeros((1000, 1000, 3), dtype = np.uint8)
xx, yy = np.meshgrid(np.linspace(-1, 1, img.shape[1]), np.linspace(-1, 1, img.shape[0]))
hue = np.arctan2(yy, xx)
new_hue = hue
new_hue[hue<0] += np.pi*2
new_hue /= 2*np.pi
sat = np.sqrt(xx**2 + yy**2)

img_hsv = np.array([colorsys.hsv_to_rgb(h, 1, 1) for h, s in zip(new_hue.ravel(), sat.ravel())]).reshape(img.shape)*255

img_hsv = img_hsv.astype(np.uint8)

fig = plt.figure(dpi = 100, figsize = (4, 3))
# gs = mgs.GridSpec(1, 10, fig)
# ax = fig.add_subplot(gs[:, :7])
ax = plt.gca()
imgs = [img_hsl, img_hsv]
cs = ['', ""]
###
hue = np.arange(0, 2*np.pi, np.pi/3)
cn = ['R', 'Y', 'G', 'C', 'B', 'M']
tc = ['k', 'w']*3
tgc = ['w', 'k']*3
axi = 0

plt.sca(ax)
ax.set_aspect('equal')
ax.axis('off')
plt.title(cs[axi])
# art = plt.imshow(imgs[axi], extent=(-1, 1, -1, 1))
art = plt.pcolormesh(xx, yy, imgs[axi])
path = mpath.Path.circle((0, 0), 1)
art.set_clip_path(path, transform = ax.transData)

arc = mpatches.Wedge((0, 0), 1, theta1 = 0, theta2 = 360-60, facecolor = (0, 0, 0, .4), fill = False, edgecolor = 'k', lw = 3)
ax.add_artist(arc)

arc = mpatches.Wedge((.05, -.05), .95, theta1 = 360-60+2, theta2 = 360-2, facecolor = (1, 1, 1, .4), fill = False, edgecolor = 'k', lw = 3)
ax.add_artist(arc)

# ax.text(.5, 0, "Missing Colors", ha = 'center', va = 'bottom', fontsize = 12, path_effects = [mpe.withStroke(foreground = 'w', linewidth = 2)])

ax.annotate("Missing Colors\n(Composite Light)", xy = (.5, -.5), xytext = (1.2, -.7), arrowprops=dict(arrowstyle = '<-'))

ax.annotate("Single  Light", xy = (-.5, .5), xytext = (1.2, .5), arrowprops=dict(arrowstyle = '<-'))

plt.xlim(-1.1, 2)
plt.ylim(-1.1, 1.1)
fig.savefig("gene/missing_colors.png", transparent=True)

# %% Two cones
wl = np.linspace(350, 750, 200)

fig = plt.figure(dpi = 200, figsize = (6, 3), layout = 'compressed')

ax = plt.subplot(1, 2, 1)
plt.grid(True, ls = '--')
for c in 'rgb':
    plt.plot(wl, signal_func(c, wl), c = c)
plt.title("三色视觉")
plt.xlabel("Wavelength (nm)")
plt.ylabel("Sensitivity")

ax = plt.subplot(1, 2, 2)
plt.grid(True, ls = '--')
plt.plot(wl, signal_func('b', wl), c = 'r')
plt.plot(wl, (signal_func('g', wl)+signal_func('r', wl))/2, c = (0.5, .5, 0))

plt.title("二色视觉")
plt.xlabel("Wavelength (nm)")
plt.ylabel("Sensitivity")

fig.savefig("gene/two_cones.png", transparent=True)

path = "data/二色视觉"
fig = plt.figure(layout = 'compressed')
for fi, file in enumerate(os.listdir(path)):
    img = np.array(plt.imread(os.path.join(path, file)), copy = True)
    old_img = np.array(img, copy = True)
    img[..., 0:2] =(np.mean(img[..., 0:2], axis = 2,keepdims = True))
    plt.subplot(2, 2, fi+1)
    plt.imshow(img)
    plt.axis('off')
fig.savefig("gene/二色.png")


fig = plt.figure(layout = 'compressed')
for fi, file in enumerate(os.listdir(path)):
    img = np.array(plt.imread(os.path.join(path, file)), copy = True)
    old_img = np.array(img, copy = True)
    plt.subplot(2, 2, fi+1)
    plt.imshow(img)
    plt.axis('off')
fig.savefig("gene/三色.png")



fig = plt.figure(layout = 'compressed')
img = np.zeros((1000, 1000, 3), dtype = np.uint8)
xx, yy = np.meshgrid(np.linspace(-1, 1, img.shape[1]), np.linspace(-1, 1, img.shape[0]))
hue = np.arctan2(yy, xx)
new_hue = hue
new_hue[hue<0] += np.pi*2
new_hue /= 2*np.pi
sat = np.sqrt(xx**2 + yy**2)
img_hsv = np.array([colorsys.hsv_to_rgb(h, 1, 1) for h, s in zip(new_hue.ravel(), sat.ravel())]).reshape(img.shape)*255
img_hsv = img_hsv.astype(np.uint8)

plt.subplot(1, 2, 1)
plt.imshow(img_hsv)
plt.axis('off')
plt.subplot(1, 2, 2)
img_hsv[..., 0:2] = (np.mean(img_hsv[..., 0:2], axis = 2, keepdims = True))
plt.imshow(img_hsv)
plt.axis('off')
fig.savefig("gene/二色hsv.png")


fig = plt.figure(layout = 'compressed')
img = plt.imread("data/红绿灯2.png")
img[..., 0:2] = (np.mean(img[..., 0:2], axis = 2, keepdims = True))
plt.imshow(img)
plt.axis('off')
fig.savefig("gene/红绿灯.png")


# %%

# %%
