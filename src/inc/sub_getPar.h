#pragma once

#include <iostream>
#include <string>
#include <cstring>
#include <map>
#include <stdlib.h>
#include "./sub_parDefine.h"

using namespace std;

int func_dict2map(string fn_par,map<string,string>&par_map);
int readpar_func(const string &fn_par, par_def &par);
int printpar_func(par_def &par);
