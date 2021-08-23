#include "TJU_SHEN_2019.h"

using namespace arma;
using namespace std;


void shen_ricker(int nw, float dt, float fre, arma::fvec &w);

int main()
{
    int n1 = 853;
    int n2 = 2317;

    float dt = 0.0005;
    int nw   = n1-3;
    float fre = 80;
    fvec w(nw,fill::zeros);
    shen_ricker(nw, dt, fre, w);

    cout<<"df is "<<1.0/dt/n1<<" Hz"<<endl;





    w.save("w.dat",raw_binary);

    fmat vel;
    vel.load("./WPIVelocityModel_nz853_nx2317_dx5_dz5.dat",arma::raw_binary);
    vel.reshape(n1,n2);

    fmat ref(size(vel),fill::zeros);
    for(int i2 = 0; i2 < n2; i2++)
    {
        for(int i1 = 1; i1 < n1; i1++)
        {
            ref(i1,i2) = (vel(i1,i2) - vel(i1-1,i2)) * 1.0
                / (vel(i1,i2) + vel(i1-1,i2));
        }
    }
    ref.save("./ref_nz853_nx2317_dx5_dz5.dat",raw_binary);

    fmat sei(size(vel),fill::zeros);
    for(int i2=0; i2<n2; i2++)
    {
        sei.col(i2) = conv(ref.col(i2),w,"same");
    }
    sei.save("./sei_nz853_nx2317_dx5_dz5.dat",raw_binary);



    fmat noise(size(sei),fill::randn);
    float sigPower = accu(sei  % sei);         // 求出信号功率
    //add noise
    noise = noise / 5.0;
    sei=sei+noise;
    sei.save("./sei_noise.dat",raw_binary);

    float noisePower=accu(noise % noise);       // 求出噪声功率
    float SNR=10*log10(sigPower/noisePower);    // 由信噪比定义求出信噪比，单位为db

    cout<<"SNR1 is : "<<SNR<<endl;
    cout<<"SNR2 is : "<<sigPower/noisePower<<endl;

    //
    Mat<cx_float> csei = fft(sei);
    Mat<cx_float> cref = fft(ref);
    Col<cx_float> cw   = fft(w);
    {
        fvec _tmp1(n1,fill::zeros);
        fvec _tmp2(n1,fill::zeros);
        fvec _tmp3(n1,fill::zeros);
        for(int i2=0; i2<n2; i2++)
        {
            _tmp1 = _tmp1 + abs(csei.col(i2));
            _tmp2 = _tmp2 + abs(cref.col(i2));
        }
        _tmp3 = abs(cw);
        _tmp1.save("csei.dat",raw_binary);
        _tmp2.save("cref.dat",raw_binary);
        _tmp3.save("cw.dat",raw_binary);
    }

    //
    fvec aimspe(n1,fill::zeros);
    {
        fvec _aimspe;
        _aimspe.load("./cwavelet_nw_8000.dat",raw_binary);
        _aimspe.reshape(8000,1);
        //interp
        {
            fvec x = linspace<fvec>(0, 1.0/5e-4, 8000);
            fvec y = _aimspe;
            fvec xx = linspace<fvec>(0, 1.0/dt, n1);
            fvec yy;
            interp1(x, y, xx, yy);
            aimspe = yy;
        }
    }
    aimspe.save("aimspe.dat",raw_binary);

    fvec ratio(n1,fill::ones);
    {
        fvec _tmp1;
        _tmp1.load("./csei.dat");
        _tmp1.reshape(n1,1);
        _tmp1 = _tmp1 / n2;
        for(int i1=0; i1<n1; i1++)
        {
            ratio(i1) = aimspe(i1) / _tmp1(i1);
        }
        //remove the value with large variance
        {
            float _mean = mean(ratio);
            int   _end  = 100/(1.0/dt/n1);
            //cout<<"mean is "<<_mean<<endl;
            for(int i1=0 ;i1<n1; i1++)
            {
                //if(ratio(i1)>_mean*5)
                //{ratio(i1)=_mean;}
            }
            //ratio(span(_end,n1-_end)).fill(1.0);
            //ratio.fill(1.0);
        }
    }
    ratio.save("ratio.dat",raw_binary);

    for(int i2=0; i2<n2; i2++)
    {
        csei.col(i2) = csei.col(i2) % ratio;
    }

    {
        fvec _tmp1(n1,fill::zeros);
        for(int i2=0; i2<n2; i2++)
        {
            _tmp1 = _tmp1 + abs(csei.col(i2));
        }
        _tmp1.save("new_csei.dat",raw_binary);

        Mat<cx_float> _tmp = ifft(csei);
        fmat _tmpsei = real(_tmp);
        _tmpsei.save("new_sei_nz853_nx2317_dx5_dz5.dat",raw_binary);
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

