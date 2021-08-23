#include "TJU_SHEN_2019.h"

using namespace std;
using namespace arma;

int main()
{
    //========================================================================//
    // Step 1 : 路径读取 & 参数设置                                           //
    //========================================================================//
    arma::Mat<float> data_in;
    arma::Mat<float> data_out;

    int nz        = 853;
    int nx        = 2317;
    string fn_in  = "./WPI_RTM.dat";
    string fn_out = "./WPI_RTM_OUT.dat";

    data_in.load(fn_in,raw_binary);
    data_in.reshape(nz,nx);

    //========================================================================//
    // Step 2 : 生成时间域Laplace滤波器                                       //
    //========================================================================//
    arma::Mat<float> lap_fli(3,3,fill::zeros);
    lap_fli = {
        { 1 , 1 , 1 },
        { 1 , -8, 1 },
        { 1 , 1 , 1 }
    };
    //lap_fli.print();
    {
        arma::Mat<float> temp_1=lap_fli % lap_fli;
        float temp_2=accu(abs(temp_1));
        //lap_fli=lap_fli/temp_2;
    }
    //lap_fli=lap_fli/arma::norm(lap_fli,"fro");
    //lap_fli.print();

    //========================================================================//
    // Step 3 : Laplace滤波                                                   //
    //========================================================================//
    arma::Mat<float> temp_mat(3,3,fill::zeros);
    data_out = data_in;
    for(int iz=1 ; iz< nz-1; iz++)
    {
        for (int ix=1; ix<nx-1; ix++)
        {
            temp_mat=data_in(span(iz-1,iz+1),span(ix-1,ix+1));
            arma::Mat<float> aaa=temp_mat % lap_fli;
            data_out(iz,ix) = accu(aaa);
        }
    }

    //========================================================================//
    // Step 4 : 文件保存                                                      //
    //========================================================================//
    data_out.save(fn_out,raw_binary);


    return 0;
}
