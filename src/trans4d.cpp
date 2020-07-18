#include "trans4d.h"
#include "utility_helpers.h"

void trans4d::GETBDY()
{
    DECLARE_COMMON_CONST
    DECLARE_COMMON_BNDRY
    // *** Obtain coordinates for vertices that form the polygons
    // *** that define the boundaries for the regions.       
    // *** Region 1 is the San Andreas fault in central California         
    // *** Region 2 is western CONUS
    // *** Region 3 is southeastern U.S.
    // *** Region 4 is eastern CONUS & southern Canada
    // *** Region 5 is mainland Alaska
    // *** Region 6 Vancouver Island
    // *** Region 7 is mainland Canada
    // *** Region 8 is the Caribbean & Central America
    // *** Region 9 is a placeholder region
    // *** Region 10 is a placeholder region
    // *** Region 11 is a placeholder region
    // *** Region 12 is a placeholder region
    // *** Region 13 is a placeholder region
    // *** Region 14 is the North American plate 
    // *** Region 15 is the Caribbean plate
    // *** Region 16 is the Pacific plate
    // *** Region 17 is the Juan de Fuca plate
    // *** Region 18 is the Cocos plate
    // *** Region 19 is the Mariana plate
    // *** REGION 20 is the Philippine Sea plate
    // *** REGION 21 is the South American plate
    // *** REGION 22 is the Nazca plate
    // *** REGION 23 is the Panama plate
    // *** REGION 24 is the North Andes plate
    int IEND = NPOINT[NMREGN + 1] - 1;
    // PRINT_VALUE(IEND)
    for(int J = 1; J <= IEND; J++)
    {
        X[J] = (X[J] * 3600.0)/RHOSEC;
        Y[J] = (Y[J] * 3600.0)/RHOSEC;

    }
}

//call this to ensure that data for boundaries, earthquake, post-seismic, and velocity grid are initialized
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

void trans4d::MODEL()
{
    // *** Obtain parameters defining crustal motion model

    DECLARE_COMMON_CONST
    DECLARE_COMMON_TIMREF
    A = 6.378137e6;
    F = 1.0 /298.25722101;
    E2 = 0.6694380022903146e-2;
    AF = A / (1.0 - F);
    EPS = F*(2.0 - F) / (pow((1.0 -F),2));
    PI = 4.0 * DATAN(1.0);
    RHOSEC = (180.0 * 3600.0) / PI;
    TWOPI = PI + PI;

    // C*** Set default reference epoch to Jan. 1, 2010
    int IYRREF = 2010;
    int IMOREF = 1;
    int IDYREF = 1;
    int MJD;
    IYMDMJ(IYRREF, IMOREF, IDYREF, MJD);
    ITREF = MJD * 24 * 60;

    GETBDY();

    DECLARE_COMMON_GRIDFILES
    //***  specify initial storage grid
    IGRID = 0;
    //specify which storage grid is needed for each region.

    //Storage grid number 1 contains the information for 
    //all regions except the Caribbean-Central America region
    //whose information is in storage file number 2
    NeededGrid[1] = 1;
    NeededGrid[2] = 1;
    NeededGrid[3] = 1;
    NeededGrid[4] = 1;
    NeededGrid[5] = 1;
    NeededGrid[6] = 1;
    NeededGrid[7] = 1;
    NeededGrid[8] = 2;
}

void::trans4d::SETRF()
{
    DECLARE_COMMON_REFCON
    //*** From blue book identifier to HTDP indentifier
    //*** WGS 72 Precise
        IRFCON[1] = 1;

    //*** WGS 84 (orig) Precise (set  equal to NAD 83)
        IRFCON[2] = 1;

    //*** WGS 72 Broadcast
        IRFCON[3] = 1;

    //*** WGS 84 (orig) Broadcast (set equal to NAD 83)
        IRFCON[4] = 1;

    //*** ITRF89
        IRFCON[5] = 3;

    //*** PNEOS 90 or NEOS 91.25 (set equal to ITRF90)
        IRFCON[6] = 4;

    //*** NEOS 90 (set equal to ITRF90)
        IRFCON[7] = 4;

    //*** ITRF91
        IRFCON[8] = 5;

    //*** SIO/MIT 92.57 (set equal to ITRF91)
        IRFCON[9] = 5;

    //*** ITRF91
        IRFCON[10] = 5;

    //*** ITRF92
        IRFCON[11] = 6;

    //*** ITRF93
        IRFCON[12] = 7;

    //*** WGS 84 (G730) Precise (set equal to ITRF91)
        IRFCON[13] = 5;

    //*** WGS 84 (G730) Broadcast (set equal to ITRF91)
        IRFCON[14] = 5;

    //*** ITRF94
        IRFCON[15] = 8;

    //*** WGS 84 (G873) Precise  (set equal to ITRF94)
        IRFCON[16] = 8;

    //*** WGS 84 (G873) Broadcast (set equal to ITRF94)
        IRFCON[17] = 8;

    //*** ITRF96
        IRFCON[18] = 8;

    //*** ITRF97
        IRFCON[19] = 9;

    //*** IGS97
        IRFCON[20] = 9;

    //*** ITRF00
        IRFCON[21] = 11;

    //*** IGS00
        IRFCON[22] = 11;

    //*** WGS 84 (G1150)
        IRFCON[23] = 11;

    //*** IGb00
        IRFCON[24] = 11;

    //*** ITRF2005
        IRFCON[25] = 14;

    //*** IGS05
        IRFCON[26] = 14;

    //*** ITRF2008 or IGS08
        IRFCON[27] = 15;

    //*** IGb08
        IRFCON[28] = 15;

    //*** ITRF2014
        IRFCON[29] = 16;

    //*** From HTDP identifier to blue book identifier
    //*** NAD 83 (set equal to WGS 84 (transit))
        JRFCON[1] = 2;

    //*** ITRF88 (set equal to ITRF89)
        JRFCON[2] = 5;

    //*** ITRF89
        JRFCON[3] = 5;

    //*** ITRF90 (set equal to NEOS 90)
        JRFCON[4] = 7;

    //*** ITRF91
        JRFCON[5] = 8;

    //*** ITRF92
        JRFCON[6] = 11;

    //*** ITRF93
        JRFCON[7] = 12;

    //*** ITRF96 (= ITRF94)
        JRFCON[8] = 18;

    //*** ITRF97
        JRFCON[9] = 19;

    //*** NA12
        JRFCON[10] = 0;

    //*** ITRF00
        JRFCON[11] = 21;

    //*** NAD 83(PACP00) or NAD 83(PA11)
        JRFCON[12] = 2;

    //*** NAD 83(MARP00) or NAD 83(MA11)
        JRFCON[13] = 2;

    //*** ITRF2005 or IGS05
        JRFCON[14] = 26;

    //*** ITRF2008 or IGS08/IGb08
        JRFCON[15] = 27;

    //*** ITRF2014
        JRFCON[16] = 29;

    //*** NA_ICE-6G
        JRFCON[17] = 0;
}

void trans4d::SETTP()
{
    DECLARE_COMMON_TRANPA
    DECLARE_COMMON_CONST
    // *** From ITRF2014 to NAD 83[2011] or NAD 83[CORS96]
    tx[1] = 1.00530e0;
    ty[1] = -1.9021e0;
    tz[1] = -.54157e0;
    dtx[1] = 0.00079e0;
    dty[1] = -.00060e0;
    dtz[1] = -.00144e0;
    rx[1] = 0.02678138 / RHOSEC;
    ry[1] = -0.00042027 / RHOSEC;
    rz[1] = 0.01093206 / RHOSEC;
    drx[1] = 0.00006667 / RHOSEC;
    dry[1] = -.00075744 / RHOSEC;
    drz[1] = -.00005133 / RHOSEC;
    scale[1] = 0.36891e-9;
    dscale[1] = -0.07201e-9;
    refepc[1] = 2010.0e0;

// *** From ITRF2014 to ITRF88
      tx[2] = 0.0254e0;
      ty[2] = -.0005e0;
      tz[2] = -.1548e0;
      dtx[2] = 0.0001e0;    
      dty[2] = -.0005e0;
      dtz[2] = -.0033e0;
      rx[2] = -.0001e0 / RHOSEC;
      ry[2] = 0.e0;
      rz[2] = -0.00026e0 / RHOSEC;
      drx[2] = 0.0e0;
      dry[2] = 0.0e0;
      drz[2] = -.00002e0 / RHOSEC;
      scale[2] = 11.29e-9;
      dscale[2] = 0.12e-9;
      refepc[2] = 2010.0e0;

// *** From ITRF2014 to ITRF89
      tx[3] = 0.0304e0;
      ty[3] = 0.0355e0;
      tz[3] = -.1308e0;
      dtx[3] = 0.0001e0;    
      dty[3] = -.0005e0;
      dtz[3] = -.0033e0;
      rx[3] = 0.0e0;
      ry[3] = 0.0e0;                
      rz[3] = -.00026e0 / RHOSEC;
      drx[3] = 0.0e0;
      dry[3] = 0.0e0;
      drz[3] = -.00002e0 / RHOSEC;
      scale[3] = 8.19e-9;
      dscale[3] = 0.12e-9;
      refepc[3] = 2010.0e0;

// *** From ITRF2014 to ITRF90
      tx[4] = 0.0254e0;
      ty[4] = 0.0115e0;
      tz[4] = -.0928e0;
      dtx[4] = 0.0001e0;    
      dty[4] = -.0005e0;
      dtz[4] = -.0033e0;
      rx[4] = 0.0e0;
      ry[4] = 0.0e0;                 
      rz[4] = -.00026e0 / RHOSEC; 
      drx[4] = 0.0e0;                  
      dry[4] = 0.0e0;              
      drz[4] = -.00002e0 / RHOSEC;
      scale[4] = 4.79e-9;
      dscale[4] = 0.12e-9;
      refepc[4] = 2010.0e0;

// *** From ITRF2014 to ITRF91
      tx[5] = 0.0274e0;
      ty[5] = 0.0155e0;
      tz[5] = -.0768e0;
      dtx[5] = 0.0001e0;    
      dty[5] = -.0005e0;
      dtz[5] = -.0033e0;
      rx[5] = 0.0e0;
      ry[5] = 0.0e0;                
      rz[5] = -.000026e0 /RHOSEC; 
      drx[5] = 0.0e0;                  
      dry[5] = 0.0e0;              
      drz[5] = -.00002e0 / RHOSEC;
      scale[5] = 4.49e-9;
      dscale[5] = 0.12e-9;
      refepc[5] = 2010.0e0;

// *** From ITRF2014 to ITRF92
      tx[6] = 0.0154e0;
      ty[6] = 0.0015e0;
      tz[6] = -.0708e0;
      dtx[6] = 0.0001e0;    
      dty[6] = -.0005e0;
      dtz[6] = -.0033e0;
      rx[6] = 0.0e0;
      ry[6] = 0.0e0;                 
      rz[6] = -.00026e0 / RHOSEC;
      drx[6] = 0.0e0;                  
      dry[6] = 0.0e0;              
      drz[6] = -.00002e0 / RHOSEC;
      scale[6] = 3.09e-9;
      dscale[6] = 0.12e-9;
      refepc[6] = 2010.0e0;

// *** From ITRF2014 to ITRF93
      tx[7] = -.0504e0;
      ty[7] = 0.0033e0;
      tz[7] = -.0602e0;
      dtx[7] = -.0028e0;
      dty[7] = -.0001e0;
      dtz[7] = -.0025e0;
      rx[7] = 0.00281e0 / RHOSEC;
      ry[7] = 0.00338e0 / RHOSEC;
      rz[7] = -.00040e0 / RHOSEC;
      drx[7] = .00011e0 / RHOSEC;
      dry[7] = .00019e0 / RHOSEC;
      drz[7] =-.00007e0 / RHOSEC;
      scale[7] = 4.29e-9;
      dscale[7] = 0.12e-9;
      refepc[7] = 2010.0e0;

// *** From ITRF2014 to ITRF94 and ITRF96
      tx[8] = 0.0074e0;
      ty[8] = -.0005e0;
      tz[8] = -.0628e0;
      dtx[8] = 0.0001e0;
      dty[8] = -.0005e0;
      dtz[8] = -.0033e0;
      rx[8] = 0.e0;
      ry[8] = 0.e0;
      rz[8] = -.00026e0 / RHOSEC;
      drx[8] = 0.0e0;
      dry[8] = 0.e0;
      drz[8] = -.00002e0 / RHOSEC;
      scale[8] = 3.80e-9;
      dscale[8] = 0.12e-9;
      refepc[8] = 2010.0e0;

// *** From ITRF2014 to ITRF97 
      tx[9] = 0.0074e0;
      ty[9] = -.0005e0;
      tz[9] = -.0628e0;
      dtx[9] = 0.0001e0;
      dty[9] = -.0005e0;
      dtz[9] = -.0033e0;
      rx[9] = 0.0e0; 
      ry[9] = 0.0e0;
      rz[9] = -.00026e0 / RHOSEC;
      drx[9] = 0.0e0; 
      dry[9] = 0.0e0; 
      drz[9] = -0.00002e0 / RHOSEC;
      scale[9] = 3.80e-9;
      dscale[9] = 0.12e-9;
      refepc[9] = 2010.0e0;

// *** From ITRF2014 to ITRF2014-PMM for North America
      tx[10] = 0.0e0;
      ty[10] = 0.0e0;
      tz[10] = 0.0e0;
      dtx[10] = 0.0e0;
      dty[10] = 0.0e0;
      dtz[10] = 0.0e0;
      rx[10] = 0.0e0; 
      ry[10] = 0.0e0;
      rz[10] = 0.0e0;
      drx[10] = +0.000024e0 / RHOSEC;
      dry[10] = -0.000694e0 / RHOSEC;
      drz[10] = -0.000063e0 / RHOSEC;
      scale[10] = 0.0e-9;
      dscale[10] = 0.0e-9;
      refepc[10] = 2010.0e0;

// *** From ITRF2014 to ITRF2000
      tx[11] = 0.0007e0;
      ty[11] = 0.0012e0;
      tz[11] = -.0261e0;
      dtx[11] = 0.0001e0;
      dty[11] = 0.0001e0;
      dtz[11] = -0.0019e0;
      rx[11] = 0.0e0; 
      ry[11] = 0.0e0; 
      rz[11] = 0.0e0; 
      drx[11] = 0.0e0; 
      dry[11] = 0.0e0;
      drz[11] = 0.0e0; 
      scale[11] = 2.12e-9;
      dscale[11] = 0.11e-9;
      refepc[11] = 2010.0e0;

// *** From ITRF2014 to PACP00 or PA11
// *** Based on the rotation rate of the Pacific plate
// ***   estimated by Beavan [2002]
      tx[12] = 0.9109e0;
      ty[12] = -2.0129e0;
      tz[12] = -0.5863e0;
      dtx[12] = 0.0001e0;
      dty[12] = 0.0001e0;
      dtz[12] = -.0019e0;
      rx[12] = 0.022749e0 / RHOSEC;
      ry[12] = 0.026560e0 / RHOSEC;
      rz[12] = -.025706e0 / RHOSEC;
      drx[12] = -.000344e0 / RHOSEC;
      dry[12] = 0.001007e0 / RHOSEC;
      drz[12] = -.002186e0 / RHOSEC;
      scale[12] = 2.12e-9;
      dscale[12] = 0.11e-9;
      refepc[12] = 2010.0e0;

// *** From ITRF2014 to MARP00 or MA11
// *** Based on the velocity of GUAM
      tx[13] = 0.9109e0;
      ty[13] = -2.0129e0;
      tz[13] = -0.5863e0;
      dtx[13] = 0.0001e0;
      dty[13] = 0.0001e0;
      dtz[13] = -.0019e0;
      rx[13] = 0.028711e0 / RHOSEC;
      ry[13] = 0.011785e0 / RHOSEC;
      rz[13] = 0.004417e0 / RHOSEC;
      drx[13] = -.000020e0 / RHOSEC;
      dry[13] = 0.000105e0 / RHOSEC;
      drz[13] = -.000347e0 / RHOSEC;
      scale[13] = 2.12e-9;
      dscale[13] = 0.11e-9;
      refepc[13] = 2010.0e0;

//*** From ITRF2014 to ITRF2005
      tx[14] = 0.0026e0;
      ty[14] = 0.0010e0;
      tz[14] = -.0023e0;
      dtx[14] = 0.0003e0;
      dty[14] = 0.0000e0;
      dtz[14] = -.0001e0;
      rx[14] = 0.0e0; 
      ry[14] = 0.0e0;
      rz[14] = 0.0e0; 
      drx[14] = 0.0e0;
      dry[14] = 0.0e0;
      drz[14] = 0.0e0;
      scale[14] = 0.92e-9;
      dscale[14] = 0.03e-9;
      refepc[14] = 2010.0e0;

//*** From ITRF2014 to ITRF2008 [also IGS08 and IGB08]
      tx[15] = 0.0016e0;
      ty[15] = 0.0019e0;
      tz[15] = 0.0024e0;
      dtx[15] = 0.0e0;
      dty[15] = 0.0e0;
      dtz[15] = -.0001e0;
      rx[15] = 0.0e0; 
      ry[15] = 0.0e0; 
      rz[15] = 0.0e0;
      drx[15] = 0.0e0; 
      dry[15] = 0.0e0; 
      drz[15] = 0.0e0; 
      scale[15] = -0.02e-9;
      dscale[15] = 0.03e-9;
      refepc[15] = 2010.0e0;

//*** From ITRF2014 to ITRF2014
      tx[16]     = 0.0e0;
      ty[16]     = 0.0e0;
      tz[16]     = 0.0e0;
      dtx[16]    = 0.0e0;                 
      dty[16]    = 0.0e0;                 
      dtz[16]    = 0.0e0;                 
      rx[16]     = 0.0e0; 
      ry[16]     = 0.0e0;
      rz[16]     = 0.0e0; 
      drx[16]    = 0.0e0;     
      dry[16]    = 0.0e0;     
      drz[16]    = 0.0e0;      
      scale[16]  = 0.0e0;
      dscale[16] =  0.0e0;               
      refepc[16] =  2010.0e0;                 

//*** From ITRF2014 to Pre-CATRF2022 [Caribbean IFVM]

      tx[17] = 0.00e0;
      ty[17] = 0.00e0;
      tz[17] = 0.00e0;
      dtx[17] = 0.000e0;
      dty[17] = 0.000e0;
      dtz[17] = 0.000e0;
      rx[17] = 0.000e0 / RHOSEC;
      ry[17] =  0.000e0 / RHOSEC;
      rz[17] =  0.000e0 / RHOSEC;
      drx[17] = -0.000000000351e0; 
      dry[17] = -0.000000004522e0;
      drz[17] = +0.000000002888e0;
      scale[17] = 0.000e0;
      dscale[17] = 0.000e0;
      refepc[17] = 2010.0e0;

}


