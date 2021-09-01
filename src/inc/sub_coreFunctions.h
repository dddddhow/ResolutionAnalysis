#pragma once
 /****************************************************************************
 *     Filename  : sub_tools.h
 *
 *     Author    : Sheng Shen
 *                  WPI,Tongji university
 *
 *     Copyright : 2021-
 *
 *     Date      : 2021.09.01
 *
 *     Description:
 *          程序所涉及到的一系列子函数的头文件，主要包括：
 *          （1） 高斯窗函数      : void shen_gauss()
 *          （2） Hanning窗函数   : void shen_hanning()
 *          （3） 高斯平滑函数    : void shen_gauss_smooth_1d()
 *          （4） Ricker子波函数  : shen_ricker()
 *
 *     Last modified:
 *
 *****************************************************************************/

#include "TJU_SHEN_2019.h"
#include "opencv2/opencv.hpp"
#include "opencv2/highgui.hpp"
#include "./sub_getPar.h"

//---------------------------------------------------------------------------//
void sub_core_waveletGeneration(par_def &par);


//---------------------------------------------------------------------------//
void sub_core_waveletAnalysis(par_def &par);


//---------------------------------------------------------------------------//
void sub_core_seiProfileGeneration(par_def &par);


//---------------------------------------------------------------------------//
void sub_core_spectrumAdjustment(par_def &par);










