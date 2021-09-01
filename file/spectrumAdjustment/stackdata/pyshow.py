#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import shen.peng.common as sp

nz = 250
nx = 131
dt = 0.002
dx = 5.0
dz = 5.0

sei   = sp.ReadDat('./sei.dat', [nz, nx], dtype='float32')
csei  = sp.ReadDat('./csei.dat', [nz, 1], dtype='float32')

nsei  = sp.ReadDat('./new_sei.dat', [nz, nx], dtype='float32')
ncsei = sp.ReadDat('./new_csei.dat', [nz, 1], dtype='float32')

# frequency spectrum show
it             = np.linspace(0, nz * dt, nz)
iw             = np.linspace(0, 1.0 / dt, nz)
fontsize_lable = 15
fontsize_ticks = 10
plt.figure()
plt.plot(iw, csei / nx, linewidth=2.0, color='blue', label='FreSpe: Sei')
plt.plot(iw, ncsei / nx, linewidth=2.0, color='red', label='FreSpe: NewSei')
plt.grid()
plt.title('Frequency Spectrum', fontsize=fontsize_lable)
plt.xlabel('Frequency (Hz)', fontsize=fontsize_lable)
plt.ylabel('Amplitude', fontsize=fontsize_lable)
plt.xticks(fontsize=fontsize_ticks)
plt.yticks(fontsize=fontsize_ticks)
plt.legend(fontsize=fontsize_lable)
plt.xlim(0, 500)

# sei profile compare show
plt.figure()
# plt.subplot(2, 1, 1)
it             = np.linspace(0, nz * dt, nz)
fontsize_lable = 15
fontsize_ticks = 10
cmap           = 'binary'
vmin           = sei.min() / 2
vmax           = sei.max() / 2
im = plt.imshow(
    sei,
    extent=(0, nx * dx, nz * dz, 0),
    vmin=vmin,
    vmax=vmax,
    cmap=cmap)
plt.axis('tight')
plt.xlabel('Lateral (m)', fontsize=fontsize_lable)
plt.ylabel('Depth (m)', fontsize=fontsize_lable)
# plt.title('Seismic profile', fontsize=fontsize_lable)
plt.xticks(fontsize=fontsize_ticks)
plt.yticks(fontsize=fontsize_ticks)


# new
plt.figure()
# plt.subplot(2, 1, 2)
it             = np.linspace(0, nz * dt, nz)
fontsize_lable = 15
fontsize_ticks = 10
cmap           = 'binary'
vmin           = nsei.min() / 5
vmax           = nsei.max() / 5
im = plt.imshow(
    nsei,
    extent=(0, nx * dx, nz * dz, 0),
    vmin=vmin,
    vmax=vmax,
    cmap=cmap)
plt.axis('tight')
plt.xlabel('Lateral (m)', fontsize=fontsize_lable)
plt.ylabel('Depth (m)', fontsize=fontsize_lable)
plt.xticks(fontsize=fontsize_ticks)
plt.yticks(fontsize=fontsize_ticks)



# trace compare
iz             = np.linspace(0, nz * dz, nz)
fontsize_lable = 15
fontsize_ticks = 10
loc1           = int(nx / 2)
plt.figure()
plt.plot(iz, sei[:, loc1], linewidth=2.0, color='blue', label='FreSpe: Sei')
plt.plot(iz, nsei[:, loc1], linewidth=2.0, color='red', label='FreSpe: NewSei')
plt.grid()
plt.title('Trace', fontsize=fontsize_lable)
plt.xlabel('Depth (m)', fontsize=fontsize_lable)
plt.ylabel('Amplitude', fontsize=fontsize_lable)
plt.xticks(fontsize=fontsize_ticks)
plt.yticks(fontsize=fontsize_ticks)
plt.legend(fontsize=fontsize_lable)



plt.show()
