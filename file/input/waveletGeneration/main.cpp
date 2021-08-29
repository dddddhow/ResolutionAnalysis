//---------------------------------------------------------------------------//
//程序目的：
//          从指定路径读取图片，转换为矩阵数据
//
//程序原理：
//
//程序参数说明：
//          armadillo库
//          opencv库
//          c++11编译
//          编译格式为： -std=c++11  -DARMA_ALLOW_FAKE_GCC -larmadillo
//                       -lopencv_core -lopencv_imgproc -lopencv_highgui
//                       -lopencv_imgcodecs
//Copyright：2020-
//          WPI TONGJI University
//Author  ：ShengShen
//Time    ：2020 08 17
//Time2   ：2020 10 04
//---------------------------------------------------------------------------//


#include <iostream>
#include <armadillo>
#include "opencv2/opencv.hpp"
#include "opencv2/highgui.hpp"

using namespace std;
using namespace arma;
using namespace cv;

//---------------------------------------------------------------------------//
void shen_gauss(int N, float theta, arma::fvec &win);
void shen_gauss_smooth_1d(arma::fvec &s_in, int nwin_in, arma::fvec &s_out);
//---------------------------------------------------------------------------//

int main()
{
    //=======================================================================//
    //      Step 1 ：参数定义
    //=======================================================================//
    //  [nz,nx] : 图片尺寸，也是矩阵尺寸
    //  img_path: 图片路径
    //  img     : 图片文件
    //  model   : 数字矩阵
    //
    //  tickness: 分辨地层厚度（米）
    //  vel     : 地层速度（米/秒）
    //
    //  w       : 子波
    //  wt      : 子波周期（秒）
    //  nw      : 子波长度（预设）
    //  dt      : 子波采样间隔（秒）
    //  df      : 子波频域采样间隔（Hz）
    cout<<"=================================================="<<endl;
    cout<<"Parameters Defination"<<endl;

    int nz,nx;
    std::string img_path = "./wavelet.png";
    cv::Mat img;
    arma::fmat model;
    arma::fvec w;

    float tickness = 30;
    float vel      = 4000;
    float wt       = tickness * 1.0 / vel * 2.0 * 2.0;
    int nw         = 80000;
    float dt,df;
    fvec ew(nw,fill::zeros);
    fvec cew(nw,fill::zeros);

    cout<<" # Expect layer tickness is : "<<tickness<<" m"<<endl;
    cout<<" # Layer velocity is : "<<vel<<" m/s"<<endl;
    cout<<" # The term of wavelet is : "<<wt<<" s"<<endl;

    //=======================================================================//
    //      Step 2 : 子波图片加载 及 转为矩阵
    //=======================================================================//
    cout<<"=================================================="<<endl;
    cout<<"Picture load and transform"<<endl;
    //  Step 2.1 :图形读取
    img= cv::imread(img_path,0);
    //cv::imshow("pic",img);

    float img_width  = img.cols;
    float img_height = img.rows;
    cout<<" # Size of Picture is (nz,nx) = (";
    cout<<img_height<<","<<img_width<<")"<<endl;

    //  Step 2.2 : 图像转为常规矩阵
    nz = img_height;
    nx = img_width;
    model.zeros(nz,nx);

    dt = wt * 1.0 / nx;
    df = 1.0 / dt / nw;
    cout<<" # dt is : "<<dt<<" s"<<endl;
    cout<<" # df is : "<<df<<" Hz"<<endl;

    for(int ix = 0; ix < nx; ix++)
    {
        for(int iz = 0; iz < nz; iz++)
        {
            model(iz,ix) = abs(img.at<uchar>(iz,ix));
        }
    }

    //  Step 2.3 : 常规矩阵提取子波向量
    w.zeros(nx);
    for(int ix = 0; ix < nx; ix++)
    {
        int flag = 0;
        for(int iz = 1; iz < nz; iz++)
        {
            if(flag==0 && model(iz,ix)-model(iz-1,ix) > 30)
            {
                flag = 1;
                w(ix) = nz-iz;
            }
        }
    }

    //  Step 2.4 : 平滑子波
    fvec _tmp = w;
    _tmp = _tmp - _tmp(0);
    w = w - w(0);
    w = w * 1.0 / w.max();
    _tmp = _tmp * 1.0 / _tmp.max();

    float maxtmp = w.max();
    float mintmp = w.min();
    cout<<" # [min,max] of wavelet is : ["<<mintmp<<","<<maxtmp<<"]"<<endl;

    shen_gauss_smooth_1d(_tmp, 11, w);

    maxtmp = w.max();
    mintmp = w.min();
    cout<<" # [min,max] of wavelet is : ["<<mintmp<<","<<maxtmp<<"]"<<endl;

    //  Step 2.5 : 拓展子波（充零）
    ew(span(int (nw/2 - nx/2),int(nw/2 + nx/2))) = w;


    //  Step 2.6 : 计算子波频谱
    cew =  abs(fft(ew));

    //=======================================================================//
    //      Step 3 : 文件保存
    //=======================================================================//
    cout<<"=================================================="<<endl;
    cout<<"Save file"<<endl;
    string fn_model = "model_nz" + std::to_string(nz)
        +"_nx"+std::to_string(nx)+".dat";
    model.save(fn_model,raw_binary);

    string fn_w_ori = "w_ori_nx" + std::to_string(nx)+".dat";
    _tmp.save(fn_w_ori,raw_binary);

    string fn_w = "w_smo_nx" + std::to_string(nx)+".dat";
    w.save(fn_w,raw_binary);

    string fn_ew = "ew_smo_nx" + std::to_string(nw)+".dat";
    ew.save(fn_ew,raw_binary);

    string fn_cew = "cew_smo_nx" + std::to_string(nw)+".dat";
    cew.save(fn_cew,raw_binary);


    cout<<"=================================================="<<endl;
    //cv::waitKey(0);
    return 0;
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

