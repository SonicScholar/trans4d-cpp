#include "trans4d.h"

void trans4d::InitBlockData() 
{
    if(!_blockDataInitialized)
    {
        InitBd();
        InitEq();
        InitPs();
        InitVl();
        _blockDataInitialized = true;
    }
}