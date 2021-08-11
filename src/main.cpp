 /****************************************************************************
 *     Author:Sheng Shen
 *              WPI,Tongji university
 *     Date:2021.08.08
 *     Filename: main.cpp
 *     Description:
 *
 *     Last modified:
 *****************************************************************************/
#include "TJU_SHEN_2019.h"

using namespace arma;
using namespace std;

/*****************************************************************************/

void shen_hanning(int N, arma::fvec &win);
void shen_gauss(int N, float theta, arma::fvec &win);

/*****************************************************************************/

int main()
{

    //======================================================================//
    //parameters definition
    //======================================================================//
    int n1   = 2400;           //samples of vertical (1000 * 8)
    int n2   = 6;              //samples of lateral
    //float d1 = 5.0;           //sampling interval of depth(m)

    int nw   = 8000;             //samples of wavelet

    float v  = 3000;            //velocity of layer (m/s)
    float dt = 0.5*1e-3;        //sampling interval of time(ms)

    fmat ref (n1, n2, fill::zeros);     //rflection coefficient
    fmat sei (size(ref), fill::zeros);  //seismic profile(traces)
    fvec w (nw, fill::zeros);           //wavelet

    //======================================================================//
    //wavelet generattion (time domain)
    //======================================================================//
    cout<<"Wavelet generation"<<endl;
    string wavelet_form_flag = "Sinc";
    //  form include :
    //  1. Ricker
    //  2. Yu
    //  3. Li
    //  4. Klauder
    //  5. Sinc
    //
    //Ricker wavelet
    if(wavelet_form_flag == "Ricker")
    {
        cout<<" Wavelet form : Ricker wavelet"<<endl;
        float amp = 1.0;            //amplitude
        float PI  = 3.1415926;      //pi
        float fre = 28;             //frequency

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
            w(span(int(nw/2),nw-1)) = _tmp;
            for(int i1=0; i1<nw/2; i1++)
            {
                w(i1) = _tmp(int(nw/2)-i1-1);
            }
        }
    }

    //Yu wavelet
    if(wavelet_form_flag == "Yu")
    {
        cout<<" Wavelet form : Yu wavelet"<<endl;
    }

    //Li wavelet
    if(wavelet_form_flag == "Li")
    {
        cout<<" Wavelet form : Li wavelet"<<endl;
    }

    //Klauder wavelet
    if(wavelet_form_flag == "Klauder")
    {
        cout<<" Wavelet form : Klauder wavelet"<<endl;
    }

    //Sinc wavelet
    if(wavelet_form_flag == "Sinc")
    {
        cout<<" Wavelet form : Sinc wavelet"<<endl;
        float _minfre    = 30;               //minimum frequency
        float _maxfre    = 120;              //maximum frequency
        float df         = 1.0 / dt / nw;   //frequency interval
        int _min         = _minfre / df;
        int _max         = _maxfre / df;
        int nwin         = _max - _min;     //length of band

        // wavelet initialize as impulse function
        w.zeros();
        w(int(nw/2))     = 1.0;

        // win : window in frequency domain
        int nside        = ceil(nwin / 10); //length of band side
        if(nside > _min){nside = _min;}
        fvec _winside(nside,fill::zeros);
        //shen_hanning(nside, _winside);

        float theta = nside / 6.0;
        shen_gauss(nside,theta,_winside);

        _winside.save("side.dat",arma::raw_binary);

        fvec win(size(w),fill::zeros);
        cout<<" Window Info : "<<endl;
        cout<<"     loc:[min,max] : ["<<_min<<","<<_max<<"]"<<endl;
        cout<<"     dfrefuency : "<<df<<endl;
        cout<<"     size : "<<nwin<<endl;
        cout<<"     nside : "<<nside<<endl;

        string flag_win = "outside";
        // 1. outside : the value 'one' start with min frequency location
        // 2. inside  : the value 'one' strat with min frequency location plus side
        if(flag_win == "inside")
        {
            win(span(_min,_min+int(nside/2)))
                = _winside(span(0,int(nside/2)));
            win(span(_max-(nside-int(nside/2)),_max-1))
                = _winside(span(int(nside/2),nside-1));
            win(span(_min+int(nside/2)
                        ,_max-(nside-int(nside/2)))).fill(1.0);
        }
        if(flag_win == "outside")
        {
            if(_min - int(nside/2) >= 0)
            {
                win(span(_min-int(nside/2),_min))
                    = _winside(span(0,int(nside/2)));
            }
            else
            {
                win(span(0,_min))
                    = _winside(span(int(nside/2)-_min,int(nside/2)));
            }
            if(_max + int(nside/2) < nw)
            {
                win(span(_max,_max+(nside-int(nside/2))-1))
                    = _winside(span(int(nside/2),nside-1));
            }
            win(span(_min,_max)).fill(1.0);
        }

        win.save("win.dat",arma::raw_binary);


        // band-limited wavelet (sinc wavelet)
        Col<cx_float> cw = fft(w);
        cw               = 2.0 * cw % win;
        w                = real(ifft(cw));
    }

    //======================================================================//
    //wavelet in frequency domain
    //======================================================================//
    cout<<"Wavelet adjust in frequency domain"<<endl;
    Col<cx_float> cw = fft(w);


    //======================================================================//
    //reflection coefficient generation (depth domain to time domain)
    //======================================================================//
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
            int nt_tmp   = int(nz_tmp * 1.0 / v * 1.0 / dt);

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
    cout<<"Save file"<<endl;
    {
        //main path of output
        string fn_out = "../file/";

        //wavelet
        string fn_wavelet = fn_out + "wavelet_nw_"+to_string(nw) + ".dat";
        w.save(fn_wavelet,arma::raw_binary);

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

        //

    }



    return 0;
}


void shen_hanning(int N, arma::fvec &win)
{
    float PI = 3.14159265358979;
    for(int i1 = 0; i1<N; i1++)
    {
        win(i1) = 0.5 * (1.0 - cos(2.0 * PI * i1 * 1.0 / (N+1)));
    }
}

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



















