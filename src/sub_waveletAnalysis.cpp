 /****************************************************************************
 *     Author:Sheng Shen
 *              WPI,Tongji university
 *     Date:2021.08.08
 *     Filename: main.cpp
 *     Description:
 *
 *     Last modified:
 *****************************************************************************/
#include "./inc/sub_tools.h"
#include "./inc/sub_coreFunctions.h"

using namespace arma;
using namespace std;

/*****************************************************************************/

void sub_core_waveletAnalysis(par_def &par)
{
    //======================================================================//
    //parameters definition
    //======================================================================//
    //  n1      : samples of vertical (300 * 8)
    //  n2      ; samples of lateral
    //  d1      : sampling interval of depth(m)
    //  v       : velocity of layer (m/s)
    //  nw      : samples of wavelet
    cout<<endl;
    //cout<<"=================================================="<<endl;
    cout<<"Parameters definition"<<endl;
    std::string fn_wavelet = par.wAPar_dic.wavelet_path;
    std::string fn_out     = par.wAPar_dic.out_path;
    int n1                 = par.wAPar_dic.n1;
    int n2                 = par.wAPar_dic.n2;
    float v                = par.wAPar_dic.vel;
    int nw                 = par.wAPar_dic.nw;
    float dtaim            = par.wAPar_dic.dtaim;
    float dt               = dtaim;

    if(dt == 999) {dt = par.wGPar_dic.dtaim;} //直接使用生成子波的采样间隔

    fmat ref (n1, n2, fill::zeros);     //rflection coefficient
    fmat sei (size(ref), fill::zeros);  //seismic profile(traces)
    fvec w (nw, fill::zeros);           //wavelet

    w.load(fn_wavelet,raw_binary);

    //======================================================================//
    //wavelet in frequency domain
    //======================================================================//
    cout<<endl;
    //cout<<"=================================================="<<endl;
    cout<<"Wavelet adjust in frequency domain"<<endl;
    Col<cx_float> cw = fft(w);


    //======================================================================//
    //reflection coefficient generation (depth domain to time domain)
    //======================================================================//
    cout<<endl;
    //cout<<"=================================================="<<endl;
    cout<<"Reflection coefficient generation"<<endl;
    {
        // _min : minimum interval of two reflection coefficient (m)
        // _max : maximum interval of two reflection coefficient (m)
        float _min = 5;
        float _max = 60;
        int _nwin  = n1 / 8;

        for(int i2 = 0; i2 < n2; i2++)
        {
            // nz_tmp : interval of two reflection coefficient (depth domain)
            // nz_tmp : interval of two reflection coefficient (time domain)
            float nz_tmp = (i2 + 1) * (_max - _min) * 1.0 / (n2 - 1);
            int nt_tmp   = int(nz_tmp * 1.0 / v * 2.0 / dt);

            cout<<"     Trace No."<<i2;
            cout<<" delta:"<<nz_tmp<<"m, nt :"<<nt_tmp<<endl;

            //
            for(int i1 = 0; i1< 4 ; i1++)
            {
                int _tmpi1 = (i1 + 0.5) * _nwin - nt_tmp / 2;
                if(_tmpi1+nt_tmp >= n1)
                {cout<<"ERROR WITH NT_TMP"<<endl;break;}
                ref(_tmpi1,i2)        = 1.0;
                ref(_tmpi1+nt_tmp,i2) = 1.0 * (4-i1) / 4.0;
            }
            //
            for(int i1 = 4; i1< 8 ; i1++)
            {
                int _tmpi1 = (i1 + 0.5) * _nwin - nt_tmp / 2;
                if(_tmpi1+nt_tmp >= n1)
                {cout<<"ERROR WITH NT_TMP"<<endl;break;}
                ref(_tmpi1,i2)        = 1.0;
                ref(_tmpi1+nt_tmp,i2) = 1.0 * (i1-8) / 4.0;
            }
        }
    }
    //======================================================================//
    //seismic profile generation (time domain) (convolution based)
    //======================================================================//
    cout<<endl;
    //cout<<"=================================================="<<endl;
    cout<<"Seismic profile generation"<<endl;
    {
        for(int i2 = 0; i2<n2; i2++)
        {
            sei.col(i2) = arma::conv(ref.col(i2),w,"same");
        }
    }
    Mat<cx_float> csei = fft(sei);


    //======================================================================//
    //save file
    //======================================================================//
    cout<<endl;
    //cout<<"=================================================="<<endl;
    cout<<"Save file"<<endl;
    //wavelet
    string fn_wav = fn_out + "wavelet_nw_"+to_string(nw) + ".dat";
    w.save(fn_wav,arma::raw_binary);

    string fn_cwavelet = fn_out + "cwavelet_nw_"+to_string(nw) + ".dat";
    fvec _tmp = abs(cw);
    _tmp.save(fn_cwavelet,arma::raw_binary);

    //reflection coefficient
    string fn_ref = fn_out + "ref_n1_"+to_string(n1)
        + "_n2_"+to_string(n2) + ".dat";
    ref.save(fn_ref,arma::raw_binary);

    //seismic profile
    string fn_sei = fn_out + "sei_n1_"+to_string(n1)
        + "_n2_"+to_string(n2) + ".dat";
    sei.save(fn_sei,arma::raw_binary);

    string fn_csei = fn_out + "csei_n1_"+to_string(n1)
        + "_n2_"+to_string(n2) + ".dat";
    fmat _tmp2 = abs(csei);
    _tmp2.save(fn_csei,arma::raw_binary);

    cout<<endl;
    //cout<<"=================================================="<<endl;

}
