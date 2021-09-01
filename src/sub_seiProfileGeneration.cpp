 /****************************************************************************
 *     Filename     : sub_seiProfileGeneration.cpp
 *
 *     Author       : Sheng Shen
 *                      WPI,Tongji university
 *
 *     Date         : 2021.03.05
 *
 *     Copyright    : 2021-
 *
 *     Description  :
 *                  readpar_func : 读取参数卡
 *                  printpar_func: 打印参数卡
 *                  func_dict2map: 结构体与容器（map）转换
 *
 *     Last modified:
 *                  Author: Sheng Shen
 *                  Date  : 2021.09.01
 *
 *****************************************************************************/

#include "TJU_SHEN_2019.h"
#include "./inc/sub_getPar.h"
#include "./inc/sub_coreFunctions.h"
#include "./inc/sub_tools.h"

using namespace arma;
using namespace std;


/*****************************************************************************/
void sub_core_seiProfileGeneration(par_def &par )
{
    //======================================================================//
    //parameters definition
    //======================================================================//
    cout<<endl;
    cout<<"Parameters definition"<<endl;
    int n1        = par.sPGPar_dic.n1;
    int n2        = par.sPGPar_dic.n2;
    float dt      = par.sPGPar_dic.dt;
    float fre     = par.sPGPar_dic.fre;
    int nw        = n1;
    string fn_vel = par.sPGPar_dic.vel_path;
    string fn_out = par.sPGPar_dic.out_path;


    fvec w(nw,fill::zeros);
    fmat vel(n1,n2,fill::zeros);
    fmat ref(size(vel),fill::zeros);
    fmat sei(size(vel),fill::zeros);

    Mat<cx_float> csei(size(vel),fill::zeros);
    Mat<cx_float> cref(size(vel),fill::zeros);
    Col<cx_float> cw(size(w),fill::zeros);

    cout<<"df is "<<1.0/dt/n1<<" Hz"<<endl;
    //======================================================================//
    //Convolution
    //======================================================================//
    cout<<endl;
    //cout<<"=================================================="<<endl;
    cout<<"Convolution"<<endl;

    //generate ricker wavelet
    shen_ricker(nw, dt, fre, w);
    cw = fft(w);

    //load velocity model
    vel.load(fn_vel,raw_binary);
    vel.reshape(n1,n2);

    //reflect coefficient
    for(int i2 = 0; i2 < n2; i2++)
    {
        for(int i1 = 1; i1 < n1; i1++)
        {
            ref(i1,i2) = (vel(i1,i2) - vel(i1-1,i2)) * 1.0
                / (vel(i1,i2) + vel(i1-1,i2));
        }
    }
    cref = fft(ref);

    //convolution
    for(int i2=0; i2<n2; i2++)
    {
        sei.col(i2) = conv(ref.col(i2),w,"same");
    }
    csei = fft(sei);

    //fmat noise(size(sei),fill::randn);
    //float sigPower = accu(sei  % sei);         // 求出信号功率
    //add noise
    //noise = noise / 5.0;
    //sei=sei+noise;
    //sei.save("./sei_noise.dat",raw_binary);
    //float noisePower=accu(noise % noise);       // 求出噪声功率
    //float SNR=10*log10(sigPower/noisePower);    // 由信噪比定义求出信噪比
    //cout<<"SNR1 is : "<<SNR<<endl;
    //cout<<"SNR2 is : "<<sigPower/noisePower<<endl;

    //======================================================================//
    //save file
    //======================================================================//
    cout<<endl;
    //cout<<"=================================================="<<endl;
    cout<<"Save file"<<endl;

    //wavelet
    w.save(fn_out+"w.dat",raw_binary);
    {
        fvec _tmp3 = abs(cw);
        _tmp3.save(fn_out+"cw.dat",raw_binary);
    }

    //seimic profile
    ref.save(fn_out+"ref.dat",raw_binary);
    sei.save(fn_out+"sei.dat",raw_binary);
    {
        fvec _tmp1(n1,fill::zeros);
        fvec _tmp2(n1,fill::zeros);
        for(int i2=0; i2<n2; i2++)
        {
            _tmp1 = _tmp1 + abs(csei.col(i2));
            _tmp2 = _tmp2 + abs(cref.col(i2));
        }
        _tmp1.save(fn_out+"csei.dat",raw_binary);
        _tmp2.save(fn_out+"cref.dat",raw_binary);
    }

}
