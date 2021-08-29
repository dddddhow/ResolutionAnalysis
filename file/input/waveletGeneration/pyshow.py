#! /usr/bin/env python3
# --------------------------------------------------------------------------- #
# 程序目的：
#
# 程序原理：
#
# 程序参数说明：
#           (1) numpy库
#           (2) matplotlib库

# Copyright：2021-
#           WPI, TONGJI University
# Author  ：Sheng Shen
# Time    ：2020 08 17
# Time2   ：2020 10 04
# Time3   ：2021 08 30
# --------------------------------------------------------------------------- #

import numpy as np
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------- #
# 参数定义 & 数据加载
# --------------------------------------------------------------------------- #
n1 = 1771
d1 = 0.02

w = np.fromfile('./w_ori_nx1771.dat',dtype='float32')
sw = np.fromfile('./w_smo_nx1771.dat',dtype='float32')

ix = np.linspace(0, n1*d1 ,n1)
label_fontsize = 15
ticks_fontsize = 10

# --------------------------------------------------------------------------- #
# 绘图
# --------------------------------------------------------------------------- #
# show :
# fig = plt.figure(figsize=(9, 4), dpi=200, constrained_layout=True)
fig = plt.figure(constrained_layout=True)
plt.plot(ix, w, linewidth=3.0, color='black',label='pick from picture')
plt.plot(ix, sw, linewidth=2.0, color='red',label='smooth')
plt.xlim(0,n1*d1)
plt.title('Wavelet Generation', fontsize=label_fontsize)
plt.legend(fontsize=label_fontsize)
plt.xlabel('Time (s)', fontsize=label_fontsize)
plt.ylabel('Amplitude', fontsize=label_fontsize)
plt.xticks(fontsize=ticks_fontsize)
plt.yticks(fontsize=ticks_fontsize)
plt.grid()



plt.show()
