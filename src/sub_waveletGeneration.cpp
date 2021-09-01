 /****************************************************************************
 *     Filename     : sub_waveletGeneration.cpp
 *
 *     Author       : Sheng Shen
 *                      WPI,Tongji university
 *
 *     Date         : 2021.08.17
 *
 *     Copyright    : 2021-
 *
 *     Description  :
*          从指定路径读取图片，转换为矩阵数据，生成地震子波。
*           armadillo库
*           opencv库
*           c++11编译
*           编译格式为： -std=c++11  -DARMA_ALLOW_FAKE_GCC -larmadillo
*                       -lopencv_core -lopencv_imgproc -lopencv_highgui
*                       -lopencv_imgcodecs
 *
 *     Last modified:
 *                  Author: Sheng Shen
 *                  Date  : 2021.09.01
 *
 *****************************************************************************/

#include "./inc/sub_tools.h"
#include "./inc/sub_coreFunctions.h"

using namespace std;
using namespace arma;
using namespace cv;

//---------------------------------------------------------------------------//

void sub_core_waveletGeneration(par_def &par)
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
    //cout<<"=================================================="<<endl;
    cout<<"Parameters Defination"<<endl;
    std::string img_path = par.wGPar_dic.img_path;
    std::string out_path = par.wGPar_dic.out_path;
    float tickness       = par.wGPar_dic.tickness;
    float vel            = par.wGPar_dic.vel;
    int nw               = par.wGPar_dic.nw;
    float wt             = tickness * 1.0 / vel * 2.0 * 2.0;

    int nz,nx;
    cv::Mat img;
    arma::fmat model;
    arma::fvec w;
    float dt,df;
    fvec ew(nw,fill::zeros);
    fvec cew(nw,fill::zeros);
    cout<<" # Expect layer tickness is : "<<tickness<<" m"<<endl;
    cout<<" # Layer velocity is : "<<vel<<" m/s"<<endl;
    cout<<" # The term of wavelet is : "<<wt<<" s"<<endl;

    //=======================================================================//
    //      Step 2 : 子波图片加载 及 转为矩阵
    //=======================================================================//
    //cout<<"=================================================="<<endl;
    cout<<endl;
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
    par.wGPar_dic.dtaim = dt;

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
            if(flag==0 && model(iz,ix)-model(iz-1,ix) > 10)
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

    shen_gauss_smooth_1d(_tmp, 17, w);

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
    //cout<<"=================================================="<<endl;
    cout<<endl;
    cout<<"Save file"<<endl;
    string fn_model = out_path+"model_nz" + std::to_string(nz)
        +"_nx"+std::to_string(nx)+".dat";
    model.save(fn_model,raw_binary);

    string fn_w_ori = out_path+"w_ori_nx" + std::to_string(nx)+".dat";
    _tmp.save(fn_w_ori,raw_binary);

    string fn_w = out_path+"w_smo_nx" + std::to_string(nx)+".dat";
    w.save(fn_w,raw_binary);

    string fn_ew = out_path+"ew_smo_nx" + std::to_string(nw)+".dat";
    ew.save(fn_ew,raw_binary);

    string fn_cew = out_path+"cew_smo_nx" + std::to_string(nw)+".dat";
    cew.save(fn_cew,raw_binary);

    //cout<<"=================================================="<<endl;
    //cv::waitKey(0);
}
