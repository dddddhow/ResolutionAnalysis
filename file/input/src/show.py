#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import shen.peng.common as sp

nz = 853
nx = 2317

ref = sp.ReadDat('./ref_nz853_nx2317_dx5_dz5.dat', [nz, nx], dtype='float32')
vel = sp.ReadDat('./WPIVelocityModel_nz853_nx2317_dx5_dz5.dat',
                 [nz, nx], dtype='float32')


plt.figure()
plt.imshow(vel)






plt.show()
