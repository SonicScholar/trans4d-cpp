#include "Trans4dCommon.h"

//initialize extern structs
BNDRY   common_bndry;
CDGRID  common_cdgrid;
CONST   common_const;
PGRID   common_pgrid;
PSGRID  common_psgrid;
QPARM   common_qparm;
REFCON  common_refcon;
TIMREF  common_timref;
TRANPA  common_tranpa;
VGRID   common_vgrid;

void InitDataFiles()
{
    InitBd();
}