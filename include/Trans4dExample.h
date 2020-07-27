#include <fstream>
struct trans4d_example
{
    /* data */
    std::fstream I2;
    std::fstream I3;
} common_files;

extern trans4d_example commonFiles;

//C++ inline hack for Fortran STOP call
//might as well do some cleanup of files along the way
inline void STOP(int exitCode)
{
    if(commonFiles.I2.is_open())
        commonFiles.I2.close();
    if(commonFiles.I3.is_open())
        commonFiles.I3.close();

    exit(exitCode);
}

void VELOC();