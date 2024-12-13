# %% Imports
from matplotlib import pyplot as plt
import numpy as np
from matplotlib import patches as mpatches, animation as manimation, patheffects as mpe
import seaborn as sns
from scipy import stats
from scipy import interpolate
import os
# %% Initialization and functions
os.makedirs('gene', exist_ok=True)
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
y = np.sin(x)
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

# %% Create con sensitivity functions by Normal core function

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





# %%
