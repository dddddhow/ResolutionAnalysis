 /****************************************************************************
 *     Filename     : main.cpp
 *
 *     Author       : Sheng Shen
 *                      WPI,Tongji university
 *
 *     Date         : 2021.08.08
 *
 *     Copyright    : 2021-
 *
 *     Description  :
 *                  程序针对地震成像剖面分辨率进行分析，主要包括：
 *                  （1） 生成地震子波：给定图片或者子波公式，获得时间域子波以
 *                                      及其振幅谱。
 *                  （2） 地震子波分析：给定子波，分析此子波的保真分辨能力和
 *                                      非保真分辨能力。
 *                  （3） 成像剖面矫正：给定目标振幅谱，将输入的成像剖面按照
 *                                      目标振幅谱进行矫正。
 *
 *     Last modified:
 *                  Author: Sheng Shen
 *                  Date  : 2021.09.01
 *
 *****************************************************************************/

#include "TJU_SHEN_2019.h"
#include "./inc/sub_tools.h"
#include "./inc/sub_coreFunctions.h"
#include "./inc/sub_getPar.h"

using namespace arma;
using namespace std;

/*****************************************************************************/
int main(int argc, char ** argv)
{
    //======================================================================//
    //parameters definition (parameter card IO)
    //======================================================================//
    cout<<endl;
    cout<<"=================================================="<<endl;
    cout<<"Parameters Card IO"<<endl;
    const string fn_par(argv[1]);
    par_def par_dic;
    readpar_func(fn_par,par_dic);
    printpar_func(par_dic);

    //======================================================================//
    //core function 1 : wavelet generation
    //======================================================================//
    if(par_dic.flag_wG == 1)
    {
        cout<<"=================================================="<<endl;
        cout<<"Wavelet Generation"<<endl;
        sub_core_waveletGeneration(par_dic);
    }

    //======================================================================//
    //core function 2 : wavelet analysis
    //======================================================================//
    if(par_dic.flag_wA == 1)
    {
        cout<<"=================================================="<<endl;
        cout<<"Wavelet Analysis"<<endl;
        sub_core_waveletAnalysis(par_dic);
    }

    //======================================================================//
    //core function 3 : seismic profile generation
    //======================================================================//
    if(par_dic.flag_sPG == 1)
    {
        cout<<"=================================================="<<endl;
        cout<<"Seicmic Profile Generation"<<endl;
        sub_core_seiProfileGeneration(par_dic);
    }

    //======================================================================//
    //core function 4 : spectrum adjustment
    //======================================================================//
    if(par_dic.flag_sA == 1)
    {
        cout<<"=================================================="<<endl;
        cout<<"Spectrum Adjustment"<<endl;
        sub_core_spectrumAdjustment(par_dic);
    }



    cout<<endl;
    cout<<"=================================================="<<endl;
    return 0;
}
