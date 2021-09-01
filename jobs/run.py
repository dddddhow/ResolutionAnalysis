#!/usr/bin/env python3

import numpy as np
import os

exe_path  = '/home/ss/WINDOWS/ResolutionAnalysis/src/bin/'
par_path  = exe_path + 'jobconfig.par'
file_path = '/home/ss/WINDOWS/ResolutionAnalysis/file/'

fn_wGPar  = 'wGPar_'
fn_wAPar  = 'wAPar_'
fn_sPGPar = 'sPGPar_'
fn_sAPar  = 'sAPar_'

par_dic = {
    # ================================================== #
    # waveletGenerationPar 子波生成函数参数
    # ================================================== #
    fn_wGPar + 'flag' : 0,
    fn_wGPar + 'img_path': file_path + 'waveletGeneration/wavelet.png',
    fn_wGPar + 'out_path': file_path + 'waveletGeneration/',
    fn_wGPar + 'vel': 3000,
    fn_wGPar + 'tickness': 10,
    fn_wGPar + 'nw': 80000,

    # ================================================== #
    # waveletAnalysisPar 子波分辨率分析函数参数
    # ================================================== #
    fn_wAPar + 'flag' : 0,
    fn_wAPar + 'wavelet_path':
    file_path + 'waveletGeneration/ew_smo_nx80000.dat',

    fn_wAPar + 'out_path':
    file_path + 'waveletAnalysis/',

    fn_wAPar + 'n1': 120000,
    fn_wAPar + 'n2': 6,
    fn_wAPar + 'vel': 5000,
    fn_wAPar + 'nw': 80000,
    fn_wAPar + 'dtaim': 999,
    # dtaim = 999 即使用 wGPar 中的子波采样间隔

    # ================================================== #
    # seiProfileGenerationPar 生成地震剖面函数参数
    # ================================================== #
    fn_sPGPar + 'flag' : 0,
    fn_sPGPar + 'vel_path':
    file_path + 'seisProfileGeneration/marmousi_nx3000_nz2300_dxdz_1.25m.bin',

    fn_sPGPar + 'out_path': file_path + 'seisProfileGeneration/',
    fn_sPGPar + 'n1': 2300,
    fn_sPGPar + 'n2': 3000,
    fn_sPGPar + 'fre': 28,
    fn_sPGPar + 'dt': 5e-4,

    # ================================================== #
    # spectrumAdjustmentPar谱平衡函数参数
    # ================================================== #
    fn_sAPar + 'flag' : 0,
    fn_sAPar + 'sei_path':
    file_path + 'spectrumAdjustment/stackdata/geo_stack_nt250_nx131_dt2ms.bin',

    fn_sAPar + 'aimSpec_path':
    file_path + 'waveletGeneration/cew_smo_nx80000.dat',

    fn_sAPar + 'out_path':
    file_path + 'spectrumAdjustment/stackdata/',

    fn_sAPar + 'n1': 250,
    fn_sAPar + 'n2': 131,
    fn_sAPar + 'dt': 0.002,
    fn_sAPar + 'dtaim': 999,
    # dtaim = 999 即使用 wGPar 中的子波采样间隔

}


# 打开一个文件
file_par = open(par_path, "w")
for key in par_dic.keys():
    line = key + '=' + str(par_dic[key])
    num = file_par.write( line + '\n')
# 关闭打开的文件
file_par.close()


exe = exe_path + 'a_out.e'
os.system("%s %s" % (exe, par_path))
