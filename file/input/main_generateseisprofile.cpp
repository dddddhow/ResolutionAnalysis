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

/*****************************************************************************/

int main()
{
    //======================================================================//
    //parameters definition
    //======================================================================//
    cout<<endl;
    cout<<"=================================================="<<endl;
    cout<<"Parameters definition"<<endl;

    int n1    = 8000;
    float v   = 1000;

    float dt  = 5e-6;
    int nw    = n1;
    float fre = 50;
    fvec w(n1,fill::zeros);
    fvec neww(n1,fill::zeros);


    float rpeak   = 1.0;        // peak value           default 1.0
    float rvalley = -0.25;      // valley value         [-1,1]
    float rmain  = 1.0;         // length of main lobe  default 1.0
    float rside  = 1.0;         // length of side lobe  [0,10]


    float T = 10.0 / v * 2.0;
    cout<<"dt is "<<T*1.0/nw<<" s"<<endl;
    cout<<"df is "<<1.0/(T*1.0/nw)/nw<<" Hz"<<endl;

    cout<<"  # [rpeak,rvalley] is :["<<rpeak<<","<<rvalley<<"]"<<endl;
    cout<<"  # [rmain,rside  ] is :["<<rmain<<","<<rside<<"]"<<endl;

    //======================================================================//
    //generate wavelet
    //======================================================================//
    cout<<endl;
    cout<<"=================================================="<<endl;
    cout<<"generate wavelet"<<endl;

    //wavelet mould (ricker)
    shen_ricker(nw, dt, fre, w);
    //adjust the wavelet
    //  1: w(0) = w(n1-1) = 0
    //  2: max = 1.0
    //  3: min = -0.5
    {
        w = w - w(0);
        float _min = w.min();
        float _max = w.max();
        cout<<" ori range is :"<<endl;
        cout<<w(0)<<" "<<w(nw-1)<<endl;
        cout<<"["<<_min<<","<<_max<<"]"<<endl;
        float _rneg = rvalley / _min;
        float _rpos = rpeak / _max;

        //w = ( w - _min ) / (_max - _min) * (rpeak - rvalley) + rvalley;

        _min = w.min();
        _max = w.max();
        cout<<" new range is :"<<endl;
        cout<<w(0)<<" "<<w(nw-1)<<endl;
        cout<<"["<<_min<<","<<_max<<"]"<<endl;
    }

    // get loction
    Col<int> loc(7,fill::zeros);
    {
        //loc(0)
        loc(0) = 0;
        //loc(1)
        for(int i1=1 ;i1<w.n_rows-1; i1++)
        {
            if(w(i1+1) > w(i1) && w(i1-1) > w(i1) && loc(1) == 0)
            {loc(1) = i1;}
        }
        //loc(2)
        for(int i1=1 ;i1<w.n_rows-1; i1++)
        {
            if(w(i1-1) < 0 &&  w(i1+1) > 0 && loc(2)==0)
            {loc(2) = i1;}
        }
        //loc(3)
        loc(3) = int(w.n_rows/2);
        //loc(4)
        for(int i1=1 ;i1<w.n_rows-1; i1++)
        {
            if(w(i1-1)>0 && w(i1+1)<0 && loc(4)==0)
            {loc(4)=i1;}
        }
        //loc(5)
        for(int i1=1 ;i1<w.n_rows-1; i1++)
        {
            if(w(i1-1)>w(i1) && w(i1+1)> w(i1) && loc(5)==0 && loc(1)<i1-10)
            {loc(5)=i1;}
        }
        //loc(6)
        loc(6)=w.n_rows-1;
        loc.print("loc is :");
    }


    //generate personal wavelet
    {
        //for (int i2=0; i2<n2 ;i2++)
        {
            float peak   = 1.0;
            float valley = peak * rvalley;
            int _lmain     = loc(4) - loc(2);
            int _lside     = int(_lmain * rside);
            int _lsideside = _lside - (loc(2)-loc(1));

            if(_lsideside <= 0)
            {cout<<"ERROR:Side lobe too small"<<endl;
                cout<<"lmain"<<_lmain<<endl;
                cout<<"lside"<<_lside<<endl;
                cout<<"lsideside"<<_lsideside<<endl;
                return 0;}

            // step 1 : width direction adjustment
            int _nw = _lmain + _lside * 2;
            fvec wtmp(_nw,fill::zeros);
            {
                fvec x = linspace<fvec>(loc(0), loc(1), loc(1)-loc(0));
                fvec y = w(span(loc(0),loc(1)-1));
                fvec xx = linspace<fvec>(loc(0), loc(1), _lsideside);
                fvec yy;
                interp1(x, y, xx, yy);
                wtmp(span(loc(0),loc(0)+_lsideside-1)) = yy;
            }
            {
                wtmp(span(loc(0)+_lsideside,
                            loc(0)+_lsideside+(loc(5)-loc(1))-1))
                    = w(span(loc(1),loc(5)-1));
            }
            {
                fvec x  = linspace<fvec>(loc(5), loc(6), loc(6)-loc(5));
                fvec y  = w(span(loc(5),loc(6)-1));
                fvec xx = linspace<fvec>(loc(5), loc(6), _lsideside);
                fvec yy;
                interp1(x, y, xx, yy);
                wtmp(span(loc(0)+_lsideside+(loc(5)-loc(1))-1,
                            loc(0)+_lmain+_lside+_lside-1))
                    =yy;
            }

            // step 2 : interpolation
            {
                fvec x   = linspace<fvec>(0, wtmp.n_rows, wtmp.n_rows);
                fvec y   = wtmp;
                fvec xx  = linspace<fvec>(0, wtmp.n_rows, n1);
                fvec yy;
                interp1(x, y, xx, yy);
                neww = yy;
            }
        }
    }
    //cout<<size(neww)<<endl;


    //
    Mat<float> cneww = abs(fft(neww));

    //======================================================================//
    //save file
    //======================================================================//
    cout<<endl;
    cout<<"=================================================="<<endl;
    cout<<"Save file"<<endl;

    w.save("w.dat",raw_binary);
    neww.save("neww.dat",raw_binary);
    cneww.save("cneww.dat",raw_binary);




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


/*
//type 2 :
//for(int i2=0 ; i2<n2 ;i2++)
{
int _lmain = loc(4) - loc(2);
int _lside = int(_lmain * rside);
int _nw = _lmain + _lside * 2;
fvec wtmp(_nw,fill::zeros);
{
fvec x = linspace<fvec>(loc(0), loc(2), loc(2)-loc(0));
fvec y = w(span(loc(0),loc(2)-1));
fvec xx = linspace<fvec>(loc(0), loc(2), _lside);
fvec yy;
interp1(x, y, xx, yy);
wtmp(span(loc(0),loc(0)+_lside-1)) = yy;
}
{
wtmp(span(loc(0)+_lside,loc(0)+_lside-1+_lmain))
= w(span(loc(2),loc(4)-1));
}
{
fvec x = linspace<fvec>(loc(4), loc(6), loc(6)-loc(4));
fvec y = w(span(loc(4),loc(6)-1));
fvec xx = linspace<fvec>(loc(4), loc(6), _lside);
fvec yy;
interp1(x, y, xx, yy);
wtmp(span(loc(0)+_lside+_lmain,loc(0)+_lmain+_lside+_lside-1))
=yy;
}
cout<<wtmp.n_rows<<endl;
wtmp.save("wtmp2.dat",raw_binary);
}



*/



