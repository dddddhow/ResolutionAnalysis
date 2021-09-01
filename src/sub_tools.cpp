 /****************************************************************************
 *     Filename     : sub_tools.cpp
 *
 *     Author       : Sheng Shen
 *                      WPI,Tongji university
 *
 *     Date         : 2021.09.01
 *
 *     Copyright    : 2021-
 *
 *     Description  :
 *                      shen_hanning : 生成hanning窗
 *                      shen_gauss   : 生成gauss窗
 *                      shen_ricker  : 生成ricker子波
 *                      shen_gauss_smooth_1d : gauss平滑
 *
 *     Last modified:
 *                  Author:
 *                  Date  :
 *
 *****************************************************************************/

#include "./inc/sub_tools.h"

using namespace std;
using namespace arma;

/*****************************************************************************/
//hanning window function
//Input :
// N    : length of window      [input]
// win  : hanning window vector [output]
void shen_hanning(int N, arma::fvec &win)
{
    float PI = 3.14159265358979;
    for(int i1 = 0; i1<N; i1++)
    {
        win(i1) = 0.5 * (1.0 - cos(2.0 * PI * i1 * 1.0 / (N+1)));
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

