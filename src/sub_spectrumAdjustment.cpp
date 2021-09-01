 /****************************************************************************
 *     Filename     : sub_spectrumAdjustment.cpp
 *
 *     Author       : Sheng Shen
 *                      WPI,Tongji university
 *
 *     Date         : 2021.09.01
 *
 *     Copyright    : 2021-
 *
 *     Description  :
 *                  成像剖面矫正：给定目标振幅谱，将输入的成像剖面按照目标
 *                  振幅谱进行矫正。
 *
 *     Last modified:
 *                  Author:
 *                  Date  :
 *
 *****************************************************************************/

#include "./inc/sub_tools.h"
#include "./inc/sub_coreFunctions.h"

using namespace arma;
using namespace std;

/*****************************************************************************/
void sub_core_spectrumAdjustment(par_def &par)
{
    //======================================================================//
    //parameters definition
    //======================================================================//
    cout<<endl;
    //cout<<"=================================================="<<endl;
    cout<<"Parameters definition"<<endl;

    int n1           = par.sAPar_dic.n1;
    int n2           = par.sAPar_dic.n2;
    float dt         = par.sAPar_dic.dt;
    string fn_sei    = par.sAPar_dic.sei_path;
    string fn_out    = par.sAPar_dic.out_path;

    string fn_aimspe = par.sAPar_dic.aimSpec_path;
    float dtaim      = par.sAPar_dic.dtaim;
    if(dtaim == 999) {dtaim = par.wGPar_dic.dtaim;} //使用生成子波的采样间隔

    fmat vel(n1,n2,fill::zeros);
    fmat sei(size(vel),fill::zeros);
    fmat nsei(size(vel),fill::zeros);

    Mat<cx_float> csei(size(vel),fill::zeros);
    Mat<cx_float> ncsei(size(vel),fill::zeros);

    fvec aimspe(n1,fill::zeros);

    cout<<" #   df is "<<1.0/dt/n1<<" Hz"<<endl;
    if(fn_sei != "")
    {
        //cout<<" Load Seismic Profile from "<<fn_sei<<endl;
        sei.load(fn_sei,raw_binary);
        sei.reshape(n1,n2);
    }

    //======================================================================//
    //Frequency domain adjustment
    //======================================================================//
    cout<<endl;
    //cout<<"=================================================="<<endl;
    cout<<"Frequency domain adjustment"<<endl;

    //fourier transform
    csei = fft(sei);

    //aim spectrum interploration
    {
        fvec _aimspe;
        _aimspe.load(fn_aimspe,raw_binary);
        //interp
        {
            fvec x = linspace<fvec>(0, 1.0/dtaim, _aimspe.n_rows);
            fvec y = _aimspe;
            fvec xx = linspace<fvec>(0, 1.0/dt, n1);
            fvec yy;
            interp1(x, y, xx, yy);
            aimspe = yy;
        }
    }

    // ratio vecter
    fvec ratio(n1,fill::ones);
    {
        fvec tmp(n1,fill::zeros);
        for(int i2=0; i2<n2; i2++)
        {
            tmp = tmp + abs(csei.col(i2));
        }
        fvec tmpsmo;
        //smooth
        shen_gauss_smooth_1d(tmp, 11, tmpsmo);
        tmpsmo = tmpsmo * 1.0 / n2;
        for(int i1=0; i1<n1; i1++)
        {
            ratio(i1) = aimspe(i1) *1.0 / tmpsmo(i1);
        }
    }

    // conservation of energy
    for(int i2=0; i2<n2; i2++)
    {
        ncsei.col(i2) = csei.col(i2) % ratio;
    }
    float maxori = (abs(csei)).max();
    float maxnew = (abs(ncsei)).max();
    ncsei = ncsei *1.0 / maxnew * maxori;

    //ifft
    {
        Mat<cx_float> _tmp = ifft(ncsei);
        nsei = real(_tmp);
    }


    //======================================================================//
    //save file
    //======================================================================//
    cout<<endl;
    //cout<<"=================================================="<<endl;
    cout<<"Save file"<<endl;

    //seimic profile
    sei.save(fn_out+"sei.dat",raw_binary);
    {
        fvec _tmp1(n1,fill::zeros);
        for(int i2=0; i2<n2; i2++)
        {
            _tmp1 = _tmp1 + abs(csei.col(i2));
        }
        _tmp1.save(fn_out+"csei.dat",raw_binary);
    }

    //new seismic profile
    nsei.save(fn_out+"new_sei.dat",raw_binary);
    {
        fvec _tmp1(n1,fill::zeros);
        for(int i2=0; i2<n2; i2++)
        {
            _tmp1 = _tmp1 + abs(ncsei.col(i2));
        }
        _tmp1.save(fn_out+"new_csei.dat",raw_binary);
    }

}
