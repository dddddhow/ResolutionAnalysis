#pragma once
//----------------------------------------------------------------------------//
//程序目的： 1、
//
//程序原理：
//
//
//程序参数说明：
//
//程序依赖关系说明：
//          1、-std=c++11 或者 更高
//
//Copyright：2020-
//          WPI TONGJI University
//Author  ：ShengShen
//Time    ：2021 03 05
//----------------------------------------------------------------------------//

#include <iostream>
#include <string>


struct waveletGenerationParDef
{
    std::string img_path;
    std::string out_path;
    float vel;
    float tickness;
    int nw;
    float dtaim;
};


struct waveletAnalysisParDef
{
    std::string wavelet_path;
    std::string out_path;
    int n1;
    int n2;
    float vel;
    int nw;
    float dtaim;
};


struct seiProfileGenerationParDef
{
    std::string vel_path;
    std::string out_path;
    int n1;
    int n2;
    float dt;
    float fre;
};


struct spectrumAdjustmentDef
{
    std::string sei_path;
    std::string aimSpec_path;
    std::string out_path;
    int n1;
    int n2;
    float dt;
    float dtaim;
};



struct par_def
{
    int flag_wG;
    int flag_wA;
    int flag_sPG;
    int flag_sA;
    waveletGenerationParDef wGPar_dic;
    waveletAnalysisParDef wAPar_dic;
    seiProfileGenerationParDef sPGPar_dic;
    spectrumAdjustmentDef sAPar_dic;
};

