 /****************************************************************************
 *     Filename     : sub_getPar.cpp
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

#include "./inc/sub_getPar.h"

/*****************************************************************************/
//Read Par
int readpar_func(const string &fn_par, par_def &par)
{
    //定义参数卡
    auto par_dic= &par;

    //将文件转化为 map（容器）
    map<string,string> par_map;
    func_dict2map(fn_par,par_map);

    //MAP遍历
    if(0)
    {
        for (auto iter=par_map.begin(); iter != par_map.end(); iter++)
        {
            std::cout<<"inner key =" << iter->first << ", val=" << iter->second
                <<std::endl;
        }
    }

    //函数开关 参数赋值
    par_dic->flag_wG  = atoi(par_map["wGPar_flag"].c_str());
    par_dic->flag_wA  = atoi(par_map["wAPar_flag"].c_str());
    par_dic->flag_sPG = atoi(par_map["sPGPar_flag"].c_str());
    par_dic->flag_sA  = atoi(par_map["sAPar_flag"].c_str());

    //子波生成函数 参数赋值
    par_dic->wGPar_dic.img_path = par_map["wGPar_img_path"];
    par_dic->wGPar_dic.out_path = par_map["wGPar_out_path"];
    par_dic->wGPar_dic.nw       = atoi(par_map["wGPar_nw"].c_str());
    par_dic->wGPar_dic.vel      = atof(par_map["wGPar_vel"].c_str());
    par_dic->wGPar_dic.tickness = atof(par_map["wGPar_tickness"].c_str());

    //子波分辨力分析函数 参数赋值
    par_dic->wAPar_dic.wavelet_path = par_map["wAPar_wavelet_path"];
    par_dic->wAPar_dic.out_path     = par_map["wAPar_out_path"];
    par_dic->wAPar_dic.n1           = atoi(par_map["wAPar_n1"].c_str());
    par_dic->wAPar_dic.n2           = atoi(par_map["wAPar_n2"].c_str());
    par_dic->wAPar_dic.nw           = atoi(par_map["wAPar_nw"].c_str());
    par_dic->wAPar_dic.vel          = atof(par_map["wAPar_vel"].c_str());
    par_dic->wAPar_dic.dtaim        = atof(par_map["wAPar_dtaim"].c_str());

    //地震剖面生成函数参数 参数赋值
    par_dic->sPGPar_dic.vel_path = par_map["sPGPar_vel_path"];
    par_dic->sPGPar_dic.out_path = par_map["sPGPar_out_path"];
    par_dic->sPGPar_dic.n1       = atoi(par_map["sPGPar_n1"].c_str());
    par_dic->sPGPar_dic.n2       = atoi(par_map["sPGPar_n2"].c_str());
    par_dic->sPGPar_dic.dt       = atof(par_map["sPGPar_dt"].c_str());
    par_dic->sPGPar_dic.fre      = atof(par_map["sPGPar_fre"].c_str());

    //谱平衡函数参数 参数赋值
    par_dic->sAPar_dic.sei_path     = par_map["sAPar_sei_path"];
    par_dic->sAPar_dic.aimSpec_path = par_map["sAPar_aimSpec_path"];
    par_dic->sAPar_dic.out_path     = par_map["sAPar_out_path"];
    par_dic->sAPar_dic.n1           = atoi(par_map["sAPar_n1"].c_str());
    par_dic->sAPar_dic.n2           = atoi(par_map["sAPar_n2"].c_str());
    par_dic->sAPar_dic.dt           = atof(par_map["sAPar_dt"].c_str());
    par_dic->sAPar_dic.dtaim        = atof(par_map["sAPar_dtaim"].c_str());


    return 0;
}



/*****************************************************************************/
// Print Par
int printpar_func(par_def &par)
{
    auto par_dic = &par;

    printf("Wavelet Generation Parameters :\n");
    printf("#   Path of Picture In   : %s \n",
            par_dic->wGPar_dic.img_path.c_str());
    printf("#   Path of Data Out     : %s \n",
            par_dic->wGPar_dic.out_path.c_str());
    std::cout<<"#   Expect Layer Tickness is : ";
    std::cout<<par_dic->wGPar_dic.tickness<<" m"<<endl;
    std::cout<<"#   Layer Velocity is : ";
    std::cout<<par_dic->wGPar_dic.vel<<" m/s"<<std::endl;


    printf("\n");
    printf("Wavelet Analysis Parameters :\n");
    printf("#   Path of Aim Wavelet In  : %s \n",
            par_dic->wAPar_dic.wavelet_path.c_str());
    printf("#   Path of Out             : %s \n",
            par_dic->wAPar_dic.out_path.c_str());
    std::cout<<"#   size of Matrix is : ["<<par_dic->wAPar_dic.n1<<",";
    std::cout<<par_dic->wAPar_dic.n2<<"]"<<endl;
    std::cout<<"#   Layer Velocity is : ";
    std::cout<<par_dic->wAPar_dic.vel<<" m/s"<<std::endl;
    if(par_dic->wAPar_dic.dtaim != 999)
    {
        std::cout<<"#   dtaim is "<<par_dic->wAPar_dic.dtaim<<" s"<<endl;
    }
    else
    {
        std::cout<<"#   dtaim is directly from wavelet generation function";
        std::cout<<endl;
    }


    printf("\n");
    printf("Seismic Profile Generation Parameters :\n");
    printf("#   Path of Velocity Model In   : %s \n",
            par_dic->sPGPar_dic.vel_path.c_str());
    printf("#   Path of Out                 : %s \n",
            par_dic->sPGPar_dic.out_path.c_str());
    std::cout<<"#   size of sei is : ["<<par_dic->sPGPar_dic.n1<<",";
    std::cout<<par_dic->sPGPar_dic.n2<<"]"<<endl;
    std::cout<<"#   dt is "<<par_dic->sPGPar_dic.dt<<" s"<<endl;
    std::cout<<"#   frequency is "<<par_dic->sPGPar_dic.fre<<" Hz"<<endl;


    printf("\n");
    printf("Spectrum Adjustment Parameters :\n");
    printf("#   Path of Sei In          : %s \n",
            par_dic->sAPar_dic.sei_path.c_str());
    printf("#   Path of AimSpec In      : %s \n",
            par_dic->sAPar_dic.aimSpec_path.c_str());
    printf("#   Path of Out             : %s \n",
            par_dic->sAPar_dic.out_path.c_str());
    std::cout<<"#   size of sei is : ["<<par_dic->sAPar_dic.n1<<",";
    std::cout<<par_dic->sAPar_dic.n2<<"]"<<endl;
    std::cout<<"#   dt is "<<par_dic->sAPar_dic.dt<<" s"<<endl;
    if(par_dic->sAPar_dic.dtaim != 999)
    {
        std::cout<<"#   dtaim is "<<par_dic->sAPar_dic.dtaim<<" s"<<endl;
    }
    else
    {
        std::cout<<"#   dtaim is directly from wavelet generation function";
        std::cout<<endl;
    }


    return 0;
}


/*****************************************************************************/
int func_dict2map(string fn_par, map<string,string> &par_map)
{
    //
    FILE *fp = fopen(fn_par.c_str(), "r");
    if (fp == NULL)
    {
        printf("failed to open %s !!!\n", fn_par.c_str());
        return 0;
    }

    //
    int TJU_LENGTH_FILENAME=1024;
    char * line_char = new char[TJU_LENGTH_FILENAME];
    char * value_char = new char[TJU_LENGTH_FILENAME];
    char * key_char = new char[TJU_LENGTH_FILENAME];
    while (1)
    {
        if (NULL == fgets(line_char, TJU_LENGTH_FILENAME, fp))
        {break;}

        string line=line_char;
        int n_line=std::strlen(line_char)-1;
        //cout<<n_line<<endl;

        size_t found = line.find_first_of("=");
        int value_start=found;
        //found = line.find_first_of("Par_");
        //int key_start=found+4;
        int key_start=0;

        string value;
        string key;

        for(int i1=0; i1<value_start-key_start; i1++)
        {
            key_char[i1]=line_char[i1+key_start];
        }
        key_char[value_start-key_start]='\0';
        key=key_char;

        for(int i1=value_start+1; i1<n_line; i1++)
        {
            value_char[i1-value_start-1]=line_char[i1];
        }
        value_char[n_line-value_start-1]='\0';
        value=value_char;

        //
        par_map[key]=value;

        if (line[0] != ' ')
        {continue;}
    }
    delete [] line_char;
    delete [] key_char;
    delete [] value_char;

    fclose(fp);

    return 0;
}
