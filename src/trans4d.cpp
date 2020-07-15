#include <iostream>
#include "Trans4dCommon.h"

using std::cout;
using std::endl;

DECLARE_COMMON_BNDRY

int main()
{
    InitDataFiles();

    cout << "Welcome to Trans4d C++ Edition" << endl;
    for (size_t i = 0; i < sizeof(common_bndry.X)/sizeof(double); i++)
    {
        cout << "common_bndry.X[" << i << "]= " << X[i] << endl;
    }
    
}