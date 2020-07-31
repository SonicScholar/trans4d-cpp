#include <fstream>

//C++ inline hack for Fortran STOP call
//might as well do some cleanup of files along the way
inline void STOP(int exitCode)
{
    exit(exitCode);
}

void TRANSFORM();
void VELOC();