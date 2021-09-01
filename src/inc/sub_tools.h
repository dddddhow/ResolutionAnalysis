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

#include <iostream>
#include "TJU_SHEN_2019.h"

//---------------------------------------------------------------------------//
//gauss window function 高斯窗函数
// N    : length of window      [input]
// win  : hanning window vector [output]
void shen_gauss(int N, float theta, arma::fvec &win);

//---------------------------------------------------------------------------//
//hanning window function Hanning窗函数
//Input :
// N    : length of window      [input]
// win  : hanning window vector [output]
void shen_hanning(int N, arma::fvec &win);

//---------------------------------------------------------------------------//
//gauss smooth function 高斯平滑函数
// s_in : signal                [input]
// N    : length of window      [input]
// s_out: signal                [output]
//
// notice :
//          shen_gauss() is needed
void shen_gauss_smooth_1d(arma::fvec &s_in, int nwin_in, arma::fvec &s_out);

//---------------------------------------------------------------------------//
//Ricker wavelet Ricker子波
//  nw : length of ricker wavelet   [input]
//  dt : sample interval            [input]
//  fre: main frequency             [input]
//  w  : ricker wavelet             [output]
void shen_ricker(int nw, float dt, float fre, arma::fvec &w);

//---------------------------------------------------------------------------//


