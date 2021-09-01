#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import shen.peng.common as sp

# ---------------------------------------------------------------------------
# parameter definition and load file
nw   = 80000
nt   = 120000
n2   = 6
# dt   = 0.5 * 1e-3
dt   = 2.25861e-05  # 1000m/s 10m
dt   = 1.12931e-05  # 2000m/s 10m
dt   = 2.25861e-05  # 2000m/s 20m
dt   = 5.64653e-06  # 4000m/s 10m
dt   = 1.12931e-05  # 4000m/s 20m
dt   = 1.69396e-05  # 4000m/s 30m
dt   = 2.25861e-05  # 4000m/s 40m
dt   = 3.76435e-06  # 6000m/s 10m
dt   = 7.5287e-06   # 6000m/s 20m
dt   = 1.12931e-05  # 6000m/s 30m


dmax = 60
dmin = 5

_fnw    = 'wavelet_nw_'+str(nw)+'.dat'
_fncw   = 'cwavelet_nw_'+str(nw)+'.dat'
_fnref  = 'ref_n1_'+str(nt)+'_n2_'+str(n2)+'.dat'
_fnsei  = 'sei_n1_'+str(nt)+'_n2_'+str(n2)+'.dat'
_fncsei = 'csei_n1_'+str(nt)+'_n2_'+str(n2)+'.dat'

w    = sp.ReadDat(_fnw, [nw, 1], dtype     = 'float32')
cw   = sp.ReadDat(_fncw, [nw, 1], dtype    = 'float32')
ref  = sp.ReadDat(_fnref, [nt, n2], dtype  = 'float32')
sei  = sp.ReadDat(_fnsei, [nt, n2], dtype  = 'float32')
csei = sp.ReadDat(_fncsei, [nt, n2], dtype = 'float32')

# ---------------------------------------------------------------------------
# variable definition
if 0:
    t             = np.linspace(0, (nt - 1) * dt, nt)
    abel_fontsize = 15
    icks_fontsize = 10
    _min           = np.min(sei)
    _max           = np.max(sei)
    fig, axs = plt.subplots( 1, n2, sharex=True, figsize=(9, 4), dpi=120,
                         constrained_layout=True)
    fig.subplots_adjust(wspace=0)
    for i2 in range(0, n2):
        axs[i2].plot(sei[:, i2], it, linewidth=1.5, color='blue')
        axs[i2].plot(ref[:, i2] * np.max(w), it, linewidth=0.5, color='red')
        axs[i2].spines['top'].set_visible(True)
        axs[i2].spines['bottom'].set_visible(True)
        axs[i2].spines['right'].set_visible(False)
        axs[i2].spines['left'].set_visible(False)
        axs[i2].set_xticks([])
        axs[i2].set_yticks([])
        axs[i2].set_ylim((nt - 1) * dt, 0)
        axs[i2].set_xlim(_min, _max)
        axs[i2].set_xlabel(str((dmax - dmin) / (n2 - 1) *
                       (i2 + 1)) + 'm', fontsize=label_fontsize)
        for j2 in range(1, 8):
            _loc = j2 * nt  / 8 * dt
            iloc = np.linspace(_loc, _loc, int(nt / 8))
            ix   = np.linspace(_min, _max, int(nt / 8))
            axs[i2].plot(ix, iloc, linewidth=1.0, color='black', linestyle='--')
            pass
        pass
    axs[n2-1].spines['right'].set_visible(True)
    axs[0].spines['left'].set_visible(True)
    pass

# ---------------------------------------------------------------------------
# variable definition
it             = np.linspace(0, (nt - 1) * dt/2, int(nt/2))
label_fontsize = 15
ticks_fontsize = 10
_min           = np.min(sei)
_max           = np.max(sei)
fig, axs = plt.subplots( 1, n2, sharex=True, figsize=(9, 4), dpi=120,
                         constrained_layout=True)
fig.subplots_adjust(wspace=0)
for i2 in range(0, n2):
    axs[i2].plot(sei[0:int(nt/2), i2], it, linewidth=1.5, color='blue')
    axs[i2].plot(ref[0:int(nt/2), i2] * np.max(w), it, linewidth=0.5, color='red')
    axs[i2].spines['top'].set_visible(True)
    axs[i2].spines['bottom'].set_visible(True)
    axs[i2].spines['right'].set_visible(False)
    axs[i2].spines['left'].set_visible(False)
    axs[i2].set_xticks([])
    axs[i2].set_yticks([])
    axs[i2].set_ylim((nt - 1) * dt / 2, 0)
    axs[i2].set_xlim(_min, _max)
    axs[i2].set_xlabel(str((dmax - dmin) / (n2 - 1) *
                       (i2 + 1)) + 'm', fontsize=label_fontsize)
    for j2 in range(1, 8):
        _loc = j2 * nt  / 8 * dt
        iloc = np.linspace(_loc, _loc, int(nt / 8))
        ix   = np.linspace(_min, _max, int(nt / 8))
        axs[i2].plot(ix, iloc, linewidth=1.0, color='black', linestyle='--')
        pass
    pass
axs[n2-1].spines['right'].set_visible(True)
axs[0].spines['left'].set_visible(True)

# ---------------------------------------------------------------------------
# variable definition
it             = np.linspace((nt - 1) * dt / 2 , (nt - 1) * dt, int(nt/2))
label_fontsize = 15
ticks_fontsize = 10
_min           = np.min(sei)
_max           = np.max(sei)
fig, axs = plt.subplots( 1, n2, sharex=True, figsize=(9, 4), dpi=120,
                         constrained_layout=True)
fig.subplots_adjust(wspace=0)
for i2 in range(0, n2):
    axs[i2].plot(sei[int(nt/2):nt, i2], it, linewidth=1.5, color='blue')
    axs[i2].plot(ref[int(nt/2):nt, i2] * np.max(w), it, linewidth=0.5, color='red')
    axs[i2].spines['top'].set_visible(True)
    axs[i2].spines['bottom'].set_visible(True)
    axs[i2].spines['right'].set_visible(False)
    axs[i2].spines['left'].set_visible(False)
    axs[i2].set_xticks([])
    axs[i2].set_yticks([])
    axs[i2].set_ylim((nt - 1) * dt, (nt - 1) * dt / 2)
    axs[i2].set_xlim(_min, _max)
    axs[i2].set_xlabel(str((dmax - dmin) / (n2 - 1) *
                       (i2 + 1)) + 'm', fontsize=label_fontsize)
    for j2 in range(1, 8):
        _loc = j2 * nt  / 8 * dt
        iloc = np.linspace(_loc, _loc, int(nt / 8))
        ix   = np.linspace(_min, _max, int(nt / 8))
        axs[i2].plot(ix, iloc, linewidth=1.0, color='black', linestyle='--')
        pass
    pass
axs[n2-1].spines['right'].set_visible(True)
axs[0].spines['left'].set_visible(True)


# ---------------------------------------------------------------------------

df = 1.0 / dt / nw
it = np.linspace(0, (nw - 1) * dt, nw)
iw = np.linspace(0, (nw - 1) * df, nw)

label_fontsize = 15
ticks_fontsize = 10

# plt.figure(figsize=(9, 4), dpi=100, constrained_layout=True)
plt.figure()
plt.subplot(2,1,1)
plt.plot(it, w, color='blue', linewidth=1.5)
plt.grid()
plt.xlim(nw*dt/2-0.2, nw*dt/2+0.2)
plt.xlabel('Time (s)', fontsize=label_fontsize)
plt.ylabel('Amplitude', fontsize=label_fontsize)
plt.xticks(fontsize=ticks_fontsize)
plt.yticks(fontsize=ticks_fontsize)

# plt.figure(figsize=(9, 4), dpi=100, constrained_layout=True)
plt.subplot(2,1,2)
plt.plot(iw, cw, color='black', linewidth=1.5)
plt.grid()
plt.xlim(0, 300)
plt.xlabel('Frequency (Hz)', fontsize=label_fontsize)
plt.ylabel('Amplitude', fontsize=label_fontsize)
plt.xticks(fontsize=ticks_fontsize)
plt.yticks(fontsize=ticks_fontsize)

# ---------------------------------------------------------------------------
plt.show()
