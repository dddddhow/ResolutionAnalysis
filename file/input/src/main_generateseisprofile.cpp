 /****************************************************************************
 *     Author:Sheng Shen
 *              WPI,Tongji university
 *     Date:2021.08.08
 *     Filename:
 *     Description:
 *
 *     Last modified:
 *****************************************************************************/
#include "TJU_SHEN_2019.h"

using namespace arma;
using namespace std;

/*****************************************************************************/
void shen_ricker(int nw, float dt, float fre, arma::fvec &w);
void shen_gauss(int N, float theta, arma::fvec &win);
void shen_gauss_smooth_1d(arma::fvec &s_in, int nwin_in, arma::fvec &s_out);

/*****************************************************************************/
int main()
{
    //======================================================================//
    //parameters definition
    //======================================================================//
    cout<<endl;
    cout<<"=================================================="<<endl;
    cout<<"Parameters definition"<<endl;

    //int n1    = 2300;
    //int n2    = 3000;
    //float dt  = 0.0005;
    //float fre = 28;
    //int nw    = n1;
    //string fn_vel    = "./vp_marmousi-ii_nx3000_nz2300_dxdz_1.25m.bin";
    //string fn_aimspe = "./cwavelet_nw_8000_60hz.bin";
    //string fn_out    = "./";

    int n1    = 250;
    int n2    = 131;
    float dt  = 0.002;
    float fre = 28;
    int nw    = n1;
    string fn_vel    = "./stackdata/geo_stack_nt250_nx131_dt2ms.bin";
    string fn_sei    = "./stackdata/geo_stack_nt250_nx131_dt2ms.bin";
    string fn_aimspe = "./cwavelet_nw_8000_60hz.bin";
    string fn_out    = "./stackdata/";


    fvec w(nw,fill::zeros);
    fmat vel(n1,n2,fill::zeros);
    fmat ref(size(vel),fill::zeros);
    fmat sei(size(vel),fill::zeros);
    fmat nsei(size(vel),fill::zeros);

    Mat<cx_float> csei(size(vel),fill::zeros);
    Mat<cx_float> ncsei(size(vel),fill::zeros);
    Mat<cx_float> cref(size(vel),fill::zeros);
    Col<cx_float> cw(size(w),fill::zeros);

    fvec aimspe(n1,fill::zeros);

    cout<<"df is "<<1.0/dt/n1<<" Hz"<<endl;

    //======================================================================//
    //Convolution
    //======================================================================//
    cout<<endl;
    cout<<"=================================================="<<endl;
    cout<<"Convolution"<<endl;

    //generate ricker wavelet
    shen_ricker(nw, dt, fre, w);

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
    //convolution
    for(int i2=0; i2<n2; i2++)
    {
        sei.col(i2) = conv(ref.col(i2),w,"same");
    }
    if(fn_sei != "")
    {
    sei.load(fn_sei,raw_binary);
    sei.reshape(n1,n2);
    }

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
    //Frequency domain adjustment
    //======================================================================//
    cout<<endl;
    cout<<"=================================================="<<endl;
    cout<<"Frequency domain adjustment"<<endl;

    //fourier transform
    csei = fft(sei);
    cref = fft(ref);
    cw   = fft(w);

    //aim spectrum interploration
    {
        fvec _aimspe;
        _aimspe.load(fn_aimspe,raw_binary);
        //_aimspe.reshape(8000,1);
        //interp
        {
            fvec x = linspace<fvec>(0, 1.0/5e-4, _aimspe.n_rows);
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
        Mat<cx_float> _tmp = ifft(ncsei*2);
        nsei = real(_tmp);
    }


    //======================================================================//
    //save file
    //======================================================================//
    cout<<endl;
    cout<<"=================================================="<<endl;
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


    return 0;
}



/*****************************************************************************/
//Ricker wavelet
//  nw : length of ricker wavelet   [input]
//  dt : sample interval            [input]
//  fre: main frequency             [input]
//  w  : ricker wavelet             [output]
void shen_ricker(int nw, float dt, float fre, arma::fvec &w)
{
    //cout<<" Wavelet form : Ricker wavelet"<<endl;
    float amp = 1.0;            //amplitude
    float PI  = 3.1415926;      //pi
    //type 1
    if(0)
    {

        for(int iw = 0; iw < nw; iw++)
        {
            float t  = dt * iw;
            float t1 = 1.0 / fre;   //双边雷克子波
            w(iw)    = amp * (1-2*PI*PI*fre*fre*(t-t1)*(t-t1))
                *exp(-PI*PI*fre*fre*(t-t1)*(t-t1));
        }
    }
    // type 2
    if(1)
    {
        fvec _tmp(int(nw/2),fill::zeros);
        for(int iw = 0; iw < nw/2; iw++)
        {
            float t  = dt * iw;
            float t1 = 0;
            _tmp(iw) = amp * (1-2*PI*PI*fre*fre*(t-t1)*(t-t1))
                *exp(-PI*PI*fre*fre*(t-t1)*(t-t1));
        }
        if(nw - int(nw/2) == int(nw/2))
        {
            w(span(int(nw/2),nw-1)) = _tmp;
            for(int i1=0; i1<nw/2; i1++)
            {
                w(i1) = _tmp(int(nw/2)-i1-1);
            }
        }
        else
        {
            w(span(nw-int(nw/2),nw-1)) = _tmp;
            for(int i1=0; i1<nw/2; i1++)
            {
                w(i1) = _tmp(int(nw/2)-i1-1);
            }
            w(int(nw/2)) = _tmp(0);
        }

    }
}


/*****************************************************************************/
//gauss window function
// N    : length of window      [input]
// win  : hanning window vector [output]
void shen_gauss(int N, float theta, arma::fvec &win)
{
    float i0 = (N - 1) * 1.0 / 2;
    float PI = 3.14159265358979;
    for(int i1 = 0; i1 < N; i1++)
    {
        win(i1) = 1.0 / (sqrt(2 * PI))
            *exp (-0.5 * ((i1 - i0) * (i1 - i0)) * 1.0 / (theta * theta));
    }
    win = win * 1.0 / abs(win).max();
}


/*****************************************************************************/
//gauss smooth function
// s_in : signal                [input]
// N    : length of window      [input]
// s_out: signal                [output]
//
// notice :
//          shen_gauss()
void shen_gauss_smooth_1d(arma::fvec &s_in, int nwin_in, arma::fvec &s_out)
{
    int n1 = s_in.n_rows;
    float theta = nwin_in/8.0;
    fvec win(nwin_in,fill::zeros);
    shen_gauss(nwin_in,theta, win);
    win = win * 1.0 / accu(win);
    fvec s2(size(s_in),fill::zeros);
    {
        fvec tmp(s_in.n_rows+nwin_in,fill::zeros);
        tmp(span(int(nwin_in/2),int(nwin_in/2)+n1-1)) = s_in;
        for(int i1=0; i1<n1; i1++)
        {
            s2(i1) = accu(tmp(span(i1,i1+nwin_in-1)) % win);
        }
    }
    s_out = s2;
}

