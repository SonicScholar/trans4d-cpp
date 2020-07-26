#include <fstream>
#include "Trans4dCommon.h"
//initialize extern structs
BNDRY       common_bndry;
CDGRID      common_cdgrid;
CONST       common_const;
GRIDFILES   common_gridfiles;
PGRID       common_pgrid;
PSGRID      common_psgrid;
QPARM       common_qparm;
REFCON      common_refcon;
TIMREF      common_timref;
TRANPA      common_tranpa;
VGRID       common_vgrid;

// Read unformatted FORTRAN sequential binary file.
size_t GridRecord::ReadRecordFromFile(FILE* dataFile)
{
    size_t result = 0;
    result += fread(&padding1, sizeof(padding1), 1, dataFile);
    result += fread(&VN, sizeof(VN), 1, dataFile);
    result += fread(&SN, sizeof(SN), 1, dataFile);
    result += fread(&VE, sizeof(VE), 1, dataFile);
    result += fread(&SE, sizeof(SE), 1, dataFile);
    result += fread(&VU, sizeof(VU), 1, dataFile);
    result += fread(&SU, sizeof(VU), 1, dataFile);
    result += fread(&padding2, sizeof(padding2), 1, dataFile);
    return result;
}