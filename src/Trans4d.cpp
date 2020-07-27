#include <iostream>
#include <fstream>
#include "Trans4d.h"
#include "UtilityHelpers.h"

bool _blockDataInitialized = false;

void trans4d::COMVEL(double& YLAT, double& YLON, int& JREGN, double& VN, double& VE, double& VU, double& SN, double& SE, double&SU)
{
    //
    // Compute the ITRF2014 velocity at a point in mm/yr
    //

    int IPLATE = 0;
    double ELON = 0;
    double HT = 0;
    double X = 0;
    double Y = 0;
    double Z = 0;
    double VX = 0;
    double VY = 0;
    double VZ = 0;
    int I = 0;
    int J = 0;

    DECLARE_COMMON_CONST
    DECLARE_COMMON_VGRID
    double WEI[2][2] = {0};
    double VEL[2][2][3] = {0};
    double STDEV[2][2][3] = {0};

    if(JREGN > NUMGRD && JREGN < NMREGN)
    {
        // Use tectonic plate model to compute velocity relative
        // to ITRF2014
        IPLATE = JREGN - NUMGRD;
        
        // *** Subtract the number of placeholder regions
        IPLATE = IPLATE - 5;

        ELON = - YLON;
        HT = 0.e0;
        TOXYZ(YLAT, ELON, HT, X, Y, Z);
        PLATVL(IPLATE, X, Y, Z, VX, VY, VZ);
        VX = VX * 1000.e0;
        VY = VY * 1000.e0;
        VZ = VZ * 1000.e0;

        TOVNEU(YLAT, ELON, VX, VY, VZ, VN, VE, VU);
        // Set standard deviations of VN, VE, and VU to 5.0 mm/yr
        SN = 5.e0;
        SE = 5.e0;
        SU = 5.e0;
    }

    else if(JREGN >= 1 && JREGN <= NUMGRD)
    {
        // C*** Get indices for the lower left hand corner of the grid
        // C*** and get the weights for the four corners
        GRDWEI (YLON, YLAT, JREGN, I, J, WEI);

        GETGRID(JREGN);

        // C*** Get the velocity vectors at the four corners
        GRDVEC(JREGN, I, J, VEL, b);

        //C*** Get standard deviations at the four corners
        GRDVEC(JREGN, I, J, STDEV, c);

        VN = WEI[0][0] * VEL[0][0][0] + WEI[0][1] * VEL[0][1][0]
           + WEI[1][0] * VEL[1][0][0] + WEI[1][1] * VEL[1][1][0];

        VE = WEI[0][0] * VEL[0][0][1] + WEI[0][1] * VEL[0][1][1]
           + WEI[1][0] * VEL[1][0][1] + WEI[1][1] * VEL[1][1][1];

        VU = WEI[0][0] * VEL[0][0][2] + WEI[0][1] * VEL[0][1][2]
           + WEI[1][0] * VEL[1][0][2] + WEI[1][1] * VEL[1][1][2];


        SN = WEI[0][0] * STDEV[0][0][0] + WEI[0][1] * STDEV[0][1][0]
           + WEI[1][0] * STDEV[1][0][0] + WEI[1][1] * STDEV[1][1][0];

        SE = WEI[0][0] * STDEV[0][0][1] + WEI[0][1] * STDEV[0][1][1]
           + WEI[1][0] * STDEV[1][0][1] + WEI[1][1] * STDEV[1][1][1];

        SU = WEI[0][0] * STDEV[0][0][2] + WEI[0][1] * STDEV[0][1][2]
           + WEI[1][0] * STDEV[1][0][2] + WEI[1][1] * STDEV[1][1][2];

        // C*** If the point is in one of the first 8 regions, then
        // c*** the velocity grids contain the ITRF2014 velocity.
    }
    else
    {
        //todo: should this exit the program here in a library?
        std::cout << "Improper region identifier in COMVEL " << std::endl;
        exit(666);
    }
    return;
}

bool trans4d::FRMXYZ(double& x, double& y, double& z, double& glat, double& glon, double& eht)
 {
     bool frmxyz;
    // *** convert x,y,z into geodetic lat, lon, and ellip. ht
    // *** ref: eq a.4b, p. 132, appendix a, osu #370
    // *** ref: geom geod notes gs 658, rapp
 
    double slat = 0;
    double clat = 0;
    double w = 0;
    double en = 0;
    //       implicit double precision(a-h,o-z)
    //       parameter(maxint=10,tol=1.d-13)
    const int maxint = 10;
    const double tol = 1.e-13;
    //       common/CONST/ a,f,e2,ep2,af,pi,twopi,rhosec
    DECLARE_COMMON_CONST
    double ae2=A*E2;
 
    // *** compute initial estimate of reduced latitude  (eht=0)
 
    double p=DSQRT(x*x+y*y);
    int icount=0;
    double tgla=z  / p /(1.e0-E2);
 
    // *** iterate to convergence, or to max # iterations
 
    label_1:
    if(icount <= maxint)
    {
        double tglax=tgla;
        tgla=z/(p-(ae2/DSQRT(1.e0+(1.e0-E2)*tgla*tgla)));
        icount=icount+1;
        if(DABS(tgla-tglax) > tol) 
            goto label_1;
 
        // *** convergence achieved
    
        frmxyz= true;
        glat=DATAN(tgla);
        double slat=DSIN(glat);
        double clat=DCOS(glat);
        glon=DATAN2(y,x);
        w=DSQRT(1.e0-E2*slat*slat);
        en=A/w;
        if(DABS(glat) <= 0.7854e0)
            eht=p/clat-en;
        else
          eht=z/slat-en+E2*en;
        glon=DATAN2(y,x);
 
        // *** too many iterations
    } 
    else
    {
        frmxyz= false;
        glat=0.e0;
        glon=0.e0;
        eht=0.e0;
    }
    return frmxyz;
 }

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

void trans4d::GETGRID(int& jregn)
{
    DECLARE_COMMON_GRIDFILES
    DECLARE_COMMON_CDGRID
    DECLARE_COMMON_VGRID

    if(NeededGrid[jregn] == IGRID) 
        return;

    if(NeededGrid[jregn] == 1) 
    {
        IGRID = 1;
        FILE* dataFile = fopen("Data4.2.5A.txt", "rb");

        GridRecord g;
        long seek = 0;
        fseek(dataFile, seek, SEEK_SET);
        for(int iregn = 1; iregn <= 7; iregn++)
        {
            for(int i = 1; i <= ICNTX[iregn] + 1; i++)
            {
                for(int j = 1; j <= ICNTY[iregn] + 1; j++)
                {
                    g.ReadRecordFromFile(dataFile);
                    int index = IUNGRD(iregn,i,j,1);
                    int index1 = index + 1;
                    int index2 = index + 2;
                    b[index] = g.VN;
                    c[index] = g.SN;
                    b[index1] = g.VE;
                    c[index1] = g.SE;
                    b[index2] = g.VU;
                    c[index2] = g.SU;
                }
            }
        }
        fclose(dataFile); 
    }
    else if (NeededGrid[jregn] == 2)
    {
        IGRID = 2;
        FILE* dataFile = fopen("Data4.2.5B.txt", "rb");

        GridRecord g;
        long seek = 0;
        for(int i = 1; i <= 609 + 1; i++)
        {
            for(int j = 1; j <= 289 + 1; j++)
            {
                fseek(dataFile, seek, SEEK_SET);
                size_t result = fread(&g, sizeof(GridRecord), 1, dataFile);

                int index = IUNGRD(jregn,i,j,1);
                int index1 = index + 1;
                int index2 = index + 2;
                b[index] = g.VN;
                c[index] = g.SN;
                b[index1] = g.VE;
                c[index1] = g.SE;
                b[index2] = g.VU;
                c[index2] = g.SU;

                seek += sizeof(GridRecord);
            }
        }
        fclose(dataFile);
      }
}

void trans4d::GETREG(double& X0, double& YKEEP, int& JREGN)
{
    DECLARE_COMMON_BNDRY
    DECLARE_COMMON_CONST

    double Y0 = TWOPI - YKEEP;
    int IR = 0;

    if (Y0 < 0.e0) Y0 = Y0 + TWOPI;
        IR = 0;

    label_1:
    IR = IR + 1;
    if(IR > NMREGN)
    {
         JREGN = 0;
         return;
    }
    
    int IBEGIN = NPOINT[IR];
    int NUMVER = NPOINT[IR + 1] - IBEGIN;
    int NTEST = 0;

    POLYIN(X0,Y0,X[IBEGIN],Y[IBEGIN], NUMVER, NTEST);

    if(NTEST == 0) 
        goto label_1;
    JREGN = IR;
}

void trans4d::GRDVEC(int JREGN, int I, int J, double (&VEL)[2][2][3], double (&B)[800000+1])
{

// C
// C********1*********2*********3*********4*********5*********6*********7**
// C
// C NAME:        GRDVEC
// C VERSION:     9302.01   (YYMM.DD)
// C WRITTEN BY:  MR. C. RANDOLPH PHILIPP
// C PURPOSE:     THIS SUBROUTINE RETRIEVES THE APPROXIMATE VALUES OF THE
// C              GRID NODE VELOCITIES FOR GRID (I,J) 
// C              
// C  INPUT PARAMETERS FROM ARGUMENT LIST:
// C  ------------------------------------
// C JREGN        ID OF GEOGRAPHIC REGION CORRESPONDING TO GRID
// C I, J         THE COORDINATES OF LOWER LEFT CORNER OF THE GRID
// C              CONTAINING THE ABOVE POSITION
// C B            THE ARRAY CONTAINING ALL THE APPROXIMATE VALUES
// C              FOR THE ADJUSTMENT
// C
// C  OUTPUT PARAMETERS FROM ARGUMENT LIST:
// C  -------------------------------------
// C VEL          A TWO BY TWO ARRAY CONTAINING THE VELOCITY VECTORS
// C              FOR THE CORNERS OF THE GRID
// C
// C  GLOBAL VARIABLES AND CONSTANTS:
// C  -------------------------------
// C NONE
// C
// C    THIS MODULE CALLED BY:   COMVEL
// C
// C    THIS MODULE CALLS:       NONE
// C
// C    INCLUDE FILES USED:      NONE
// C
// C    COMMON BLOCKS USED:      NONE     
// C
// C    REFERENCES:  SEE RICHARD SNAY
// C
// C    COMMENTS:
// C
// C********1*********2*********3*********4*********5*********6*********7**
// C    MOFICATION HISTORY:
// C::9302.11, CRP, ORIGINAL CREATION FOR DYNAP
// C::9712.05, RAS, MODIFIED FOR HTDP (version 2.2)
// C********1*********2*********3*********4*********5*********6*********7**

    for(int II = 0; II <= 1; II++)
    {
        for(int IJ = 0; IJ <=1; IJ++)
        {
            for(int IVEC = 1; IVEC <= 3; IVEC++)
            {
                int INDEX = IUNGRD(JREGN, I + II, J + IJ, IVEC);
                VEL[II][IJ][IVEC-1] = B[INDEX];
            }
        }
    }

    return;
}

void trans4d::GRDWEI(double const& YLON, double const& YLAT, int const& JREGN, int& I, int&J, double (&WEI)[2][2] )
{
// C
// C********1*********2*********3*********4*********5*********6*********7**
// C
// C NAME:        GRDWEI
// C VERSION:     9302.01   (YYMM.DD)
// C WRITTEN BY:  MR. C. RANDOLPH PHILIPP
// C PURPOSE:     THIS SUBROUTINE RETURNS THE INDICES OF THE LOWER-LEFT
// C              HAND CORNER OF THE GRID CELL CONTAINING THE POINT
// C              AND COMPUTES NORMALIZED WEIGHTS FOR 
// C              BI-LINEAR INTERPOLATION OVER A PLANE
// C              
// C  INPUT PARAMETERS FROM ARGUMENT LIST:
// C  ------------------------------------
// C YLON         LONGITUDE OF POINT IN RADIANS, POSITIVE WEST
// C YLAT         LATITUDE OF POINT IN RADIANS, POSITIVE NORTH
// C JREGN        ID OF GEOGRAPHIC REGION CONTAINING POINT
// C
// C  OUTPUT PARAMETERS FROM ARGUMENT LIST:
// C  -------------------------------------
// C I, J         THE COORDINATES OF LOWER LEFT CORNER OF THE GRID
// C              CONTAINING THE ABOVE POSITION
// C WEI          A TWO BY TWO ARRAY CONTAINING THE NORMALIZED WEIGHTS
// C              FOR THE CORNER VECTORS
// C
// C  GLOBAL VARIABLES AND CONSTANTS:
// C  -------------------------------
// C NONE
// C
// C    THIS MODULE CALLED BY:   COMVEL
// C
// C    THIS MODULE CALLS:       NONE
// C
// C    INCLUDE FILES USED:      NONE
// C
// C    COMMON BLOCKS USED:      /CDGRID/, /CONST/
// C
// C    REFERENCES:  SEE RICHARD SNAY
// C
// C    COMMENTS:
// C
// C********1*********2*********3*********4*********5*********6*********7**
// C    MOFICATION HISTORY:
// C::9302.11, CRP, ORIGINAL CREATION FOR DYNAP
// C::9511.09, RAS, MODIFIED FOR HTDP
// C::9712.05, RAS, MODIFIED TO ACCOUNT FOR MULTIPLE GRIDS
// C********1*********2*********3*********4*********5*********6*********7**
    
// C**** COMPUTES THE WEIGHTS FOR AN ELEMENT IN A GRID

    DECLARE_COMMON_CDGRID
    DECLARE_COMMON_CONST

// C*** Convert input coordinates to degrees
    double POSX = (TWOPI - YLON) * 180.e0 / PI;
    double POSY = YLAT * 180.e0 / PI;

// C*** Obtain indices for the lower-left corner of the cell
// C*** containing the point
    double STEPX = (GRDUX[JREGN] - GRDLX[JREGN]) / ICNTX[JREGN];
    double STEPY = (GRDUY[JREGN] - GRDLY[JREGN]) / ICNTY[JREGN];
    I = IDINT((POSX - GRDLX[JREGN])/STEPX) + 1;
    J = IDINT((POSY - GRDLY[JREGN])/STEPY) + 1;

// C*** Compute the limits of the grid cell 
    double GRLX = GRDLX[JREGN] + (I - 1) * STEPX;
    double GRUX = GRLX + STEPX;       
    double GRLY = GRDLY[JREGN] + (J - 1) * STEPY;             
    double GRUY = GRLY + STEPY;

//porting note: leaving original array assignments for reference
//made array assignments 0 based since the WEI array is not global
// C*** Compute the normalized weights for the point               
//       DENOM = (GRUX - GRLX) * (GRUY - GRLY)
//       WEI(1,1) = (GRUX - POSX) * (GRUY - POSY) / DENOM
//       WEI(2,1) = (POSX - GRLX) * (GRUY - POSY) / DENOM
//       WEI(1,2) = (GRUX - POSX) * (POSY - GRLY) / DENOM
//       WEI(2,2) = (POSX - GRLX) * (POSY - GRLY) / DENOM
    double DENOM = (GRUX - GRLX) * (GRUY - GRLY);
    WEI[0][0] = (GRUX - POSX) * (GRUY - POSY) / DENOM;
    WEI[1][0] = (POSX - GRLX) * (GRUY - POSY) / DENOM;
    WEI[0][1] = (GRUX - POSX) * (POSY - GRLY) / DENOM;
    WEI[1][1] = (POSX - GRLX) * (POSY - GRLY) / DENOM;
}

void trans4d::GTOVEL(double const& YLAT, double const& YLON,  double const& EHT,
double& VN, double& VE, double& VU, double& VX, double& VY, double& VZ, int& JREGN, int const& IOPT,
double& SN, double& SE, double& SU, double& SX, double& SY, double& SZ)
{
    // *** Compute velocity in appropriate reference frame for point with
    // *** latitude YLAT (radians), longitude YLON (radians, positive west)
    // *** and ellipsoid height EHT (meters) in this reference frame.

    //*** Get reference latitude RLAT and reference longitude RLON in ITRF2014
    double RLAT, RLON, ELON;
    double X, Y, Z;
    double DATE;
    double EHT14;

    if(IOPT == 16)
    {
          RLAT = YLAT;
          RLON = YLON;
    }
    else 
    {
        ELON = -YLON;
        TOXYZ(YLAT,ELON,EHT,X,Y,Z);
        DATE = 2010.0e0;
        XTOITRF2014(X,Y,Z,RLAT,RLON,EHT14,DATE,IOPT);
    }

    // Get velocity in ITRF2014 
    GETREG(RLAT,RLON,JREGN);
    if (JREGN == 0) 
        return;

    COMVEL(RLAT,RLON,JREGN,VN,VE,VU,SN,SE,SU);

    ELON = -YLON;
    TOVXYZ(YLAT,ELON,VN,VE,VU,VX,VY,VZ);
    to_std_dev_xyz_velocity(YLAT,ELON,SN,SE,SU,SX,SY,SZ);

    // Convert velocity into another reference frame if needed
    if(IOPT != 16)
    {
         VTRANF(X,Y,Z,VX,VY,VZ, 16, IOPT);
         TOVNEU(YLAT, ELON, VX, VY, VZ, VN, VE, VU);
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
int trans4d::IUNGRD(int const& IREGN, int const& I, int const& J, int const& IVEC)
{
    DECLARE_COMMON_CDGRID

    int IUNGRD = NBASE[IREGN] +
           3 * ((J - 1) * (ICNTX[IREGN] + 1) +  (I - 1)) + IVEC;

    return IUNGRD;
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

void trans4d::PLATVL(int& IPLATE, double& X, double& Y, double& Z, double& VX, double& VY, double& VZ)
{
    // *** Compute the ITRF2014 velocity at point on plate = IPLATE
    // ***    with coordinates X, Y, Z (in meters)
    // ***    The resulting velocities--VX, VY, and VZ--will be in meters/yr
    // ***    References 
    // ***     Altamimi et al. 2017 = JGR (Paper on ITRF2014 plate motion)
    // ***     Kreemer e al. 2014 = Geochem. Geophys. & Geosyst., vol 15
    // ***     DeMets et al. 2010 = Geophysical Journal Int'l, vol 181, 
    // ***     Snay 2003 = SALIS, Vol 63, No 1 (Paper on Frames for Pacific)
    // ***     Bird 2003 = Geochem. Geophys. & Geosyst., (Paper on plate boundaries)

    // *** IPLATE = 1 --> North America (from Altamimi et al. 2017)
    // ***          2 --> Caribbean (from Kreemer et al. 2014)
    // ***          3 --> Pacific (from Altamimi et al. 2017)
    // ***          4 --> Juan de Fuca (from DeMets et al. 2010)
    // ***          5 --> Cocos (from DeMets et al. 2010)
    // ***          6 --> Mariana (from Snay, 2003)
    // ***          7 --> Philippine Sea (from Kreemer et al. 2014)
    // ***          8 --> South America (from Altamimi et al. 2017)
    // ***          9 --> Nazca (from Altamimi et al. 2017)
    // ***         10 --> Panama (from Kreemer et al. 2014)
    // ***         11 --> North Andes (from Bird 2003)

    //C++ Porting note: arrays WX, WY, & WZ are 1 based with an extra 0 in the beginning
    //making sure that anything that access their elements by index does not need -1 in it.

    // ORIGINAL FORTRAN DATA 
    //   DATA WX /0.116D-9,  -0.675D-9,-1.983D-9, 
    //  1         6.636D-9, -10.380D-9,-0.097D-9,
    //  2         9.221D-9,  -1.309D-9,-1.614D-9, 
    //  3         2.088D-9,  -1.872D-9 /
    //   DATA WY /-3.365D-9, -3.826D-9, 5.076D-9,
    //  1         11.761D-9,-14.900D-9, 0.509D-9,
    //  2         -4.963D-9, -1.459D-9,-0.679D-9, 
    //  3        -23.037D-9, -1.285D-9 /
    //   DATA WZ /-0.305D-9,  2.910D-9,-10.516D-9, 
    //  1        -10.630D-9,  9.133D-9, -1.682D-9,
    //  2        -11.554D-9, -0.679D-9,  7.868D-9, 
    //  3          6.729D-9, -0.067D-9 /  

    static const double WX[11 + 1] = { 0, 0.116e-9, -0.675e-9, -1.983e-9, 
                                        6.636e-9, -10.380e-9, -0.097e-9, 
                                        9.221e-9,  -1.309e-9, -1.614e-9,
                                        2.088e-9, -1.872e-9};

    static const double WY[11 + 1] = { 0, -3.365e-9, -3.826e-9, 5.076e-9,
                                        11.761e-9, -14.900e-9, 0.509e-9,
                                        -4.963e-9, -1.459e-9,-0.679e-9,
                                        -23.037e-9, -1.285e-9};

    static const double WZ[11 + 1] = { 0,  -0.305e-9,  2.910e-9, -10.516e-9,
                                        -10.630e-9,  9.133e-9, -1.682e-9,
                                        -11.554e-9, -0.679e-9,  7.868e-9,
                                        6.729e-9, -0.067e-9};

    if(IPLATE <= 0 || IPLATE > 11)
    {
        std::cout << " Improper plate ID in PLATVL " << IPLATE;
        exit(666);
    }

    VX = -WZ[IPLATE] * Y + WY[IPLATE] * Z;
    VY =  WZ[IPLATE] * X - WX[IPLATE] * Z;
    VZ = -WY[IPLATE] * X + WX[IPLATE] * Y;

    // *** The parameters--WX, WY, and WZ--refer to ITRF2000
    // *** for the Mariana Plate (Snay, 2003). Hence,
    // *** for this plate, VX, VY, and VZ, correspond to ITRF2000.
    // *** The following code converts these to ITRF2014 velocities for
    // *** this plate.
    if(IPLATE == 6)
    {
         VX = VX*1000.e0;
         VY = VY*1000.e0;
         VZ = VZ*1000.e0;
         VTRANF(X, Y, Z, VX, VY, VZ, 11, 15);
         VX = VX/1000.e0;
         VY = VY/1000.e0;
         VZ = VZ/1000.e0;
    }
// *** The following translations rates are added per Altamimi et al. (2012)
// *** for the other six plates
    else
    {
        VX = 0.00041e0 + VX;
        VY = 0.00022e0 + VY;
        VZ = 0.00041e0 + VZ;
    }
}

void trans4d::POLYIN(double& X0, double& Y0, double& X, double& Y, int& N, int& NPC)
{
    double* arr_X = &X;
    double* arr_Y = &Y;

    // C     SUBROUTINE TO DETERMINE IF A POINT AT (X0,Y0) IS INSIDE OR
    // C     OUTSIDE OF A CLOSED FIGURE DESCRIBED BY A SEQUENCE OF CONNECTED
    // C     STRAIGHT LINE SEGMENTS WITH VERTICES AT X, Y.
    // C
    // C     INPUT -
    // C         X0, Y0    COORDINATES OF A POINT TO BE TESTED
    // C                    Y0 corresponds to longitude and must be a number
    // C                    between 0.0 and 2*PI
    // C         X, Y      ARRAYS CONTAINING THE VERTICES, IN ORDER, OF A
    // C                   CLOSED FIGURE DESCRIBED BY STRAIGHT LINE SEGMNENTS.
    // C                   FOR EACH 'I', THE STRAIGHT LINE FROM (XI),Y(I)) TO
    // C                   TO (X(I+1),Y(I+1)), IS AN EDGE OF THE FIGURE.
    // C         N         DIMENSION OF X AND Y, NUMBER OF VERTICES, AND NUMBER
    // C                   OF STRAIGHT LINE SEGMENTS IN FIGURE.
    // C     OUTPUT -
    // C         NPC       NPC=0 WHEN X0,Y0 IS OUTSIDE OF FIGURE DESCRIBED
    // C                   BY X,Y
    // C                   NPC=1 WHEN X0,Y0 IS INSIDE FIGURE
    // C                   NPC=2 WHEN X0,Y0 IS ON BORDER OF FIGURE
    // C     METHOD -
    // C     A COUNT IS MADE OF THE NUMBER OF TIMES THE LINE FROM (X0,Y0) TO
    // C     (X0,+ INFINITY) CROSSES THE BORDER OF THE FIGURE. IF THE COUNT
    // C     IS ODD, THE POINT IS INSIDE; IF THE COUNT IS EVEN THE POINT
    // C     IS OUTSIDE.
    // C     LIMITATIONS -
    // C     NONE. THE PROGRAM LOGIC IS VALID FOR ALL CLOSED FIGURES,
    // C     NO MATTER HOW COMPLEX.
    // C     ACCURACY -
    // C     MAINTAINS FULL ACCURACY OF INPUT COORDINATES.
    // C
    //       IMPLICIT INTEGER(I-N)
    //       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    //       DIMENSION X(N),Y(N)
    //       DATA I6/6/
    int IS=0;
    int IP=0;
    int IL=0;
    double XL=0;
    double YL=0;
    int IP1=0;
    int IPN=0;
    int II=0;
    int I=0;
    NPC=0;
    // C
    // C     FIND STARTING POINT WHERE X(I).NE.X0

    label_10:
    
    //C++ Port, since IP is being used to find the starting point in
    //an array, it should be 0 initially.
    //IP = IP + 1;
    

    switch(IF_ARITHMETIC(arr_X[IP]-X0))
    {
        case -1: goto label_15;
        case 0: goto label_12;
        default: goto label_16;
    } 
    label_12:
    if(IP <= N) goto label_10;
    std::cout << "0 POLYGON INPUT ERROR - ALL POINTS ON LINE X = X0" << std::endl;
    exit(666);
    //       STOP
    label_15:
    IL = -1;
    goto label_20;

    label_16:
    IL=1;

    label_20:
    XL = arr_X[IP];
    YL = arr_Y[IP];
    // C
    // C     SET UP SEARCH LOOP
    // C
    IP1=IP+1;
    IPN=IP+N;
    //       DO 100 II=IP1,IPN
    for(II=IP1; II<=IPN; II++)
    {
        I=II;
        if(I >= N) 
        I=I-N;
        // IF(IL) 30,50,40
        switch(IF_ARITHMETIC(IL))
        {
            case -1: goto label_30;
            case 0: goto label_50;
            default: goto label_40;
        }

    
        //    30 IF(X(I)-X0) 90,32,34
        label_30:
        switch(IF_ARITHMETIC(arr_X[I]-X0))
        {
            case -1: goto label_90;
            case 0: goto label_32;
            default: goto label_34;
        }

        label_32:
        IS=-1;
        goto label_60;

        label_34:
        IL=1;
        goto label_80;
        
        label_40:
        //    40 IF(X(I)-X0) 42,44,90
        switch(IF_ARITHMETIC(arr_X[I]-X0))
        {
            case -1: goto label_42;
            case 0: goto label_44;
            default: goto label_90;
        }

        label_42:
        IL=-1;
        goto label_80;

        label_44:
        IS=1;
        goto label_60;

        label_50:
        //    50 IF(X(I)-X0) 52,55,54
        switch(IF_ARITHMETIC(arr_X[I]-X0))
        {
            case -1: goto label_52;
            case 0: goto label_55;
            default: goto label_54;
        }

        label_52:
        IL=-1;
        //       IF(IS) 90,140,80
        switch(IF_ARITHMETIC(IS))
        {
            case -1: goto label_90;
            case 0: goto label_140;
            default: goto label_80;
        }
        label_54:
        IL=1;
        //       IF(IS) 80,140,90
            switch(IF_ARITHMETIC(IS))
        {
            case -1: goto label_80;
            case 0: goto label_140;
            default: goto label_90;
        }

        label_55:
        //    55 IF(Y(I)-Y0) 57,120,58
        switch(IF_ARITHMETIC(arr_Y[I]-Y0))
        {
            case -1: goto label_57;
            case 0: goto label_120;
            default: goto label_58;
        }
        
        label_57:
        //    57 IF(YL-Y0) 90,120,120
        switch(IF_ARITHMETIC(YL-Y0))
        {
            case -1: goto label_90;
            case 0: goto label_120;
            default: goto label_120;
        }

        label_58:
        //    58 IF(YL-Y0) 120,120,90
        switch(IF_ARITHMETIC(YL-Y0))
        {
            case -1: goto label_120;
            case 0: goto label_120;
            default: goto label_90;
        }

        // C
        label_60:
        //    60 IL=0
        IL=0;
        //       IF(Y(I)-Y0) 90,120,90
            switch(IF_ARITHMETIC(arr_Y[I]-Y0))
        {
            case -1: goto label_90;
            case 0: goto label_120;
            default: goto label_90;
        }
        label_80:
        //    80 IF(YL-Y0+(Y(I)-YL)*(X0-XL)/(X(I)-XL)) 90,120,85
        switch(IF_ARITHMETIC(YL-Y0+(arr_Y[I]-YL)*(X0-XL)/(arr_X[I]-XL)))
        {
            case -1: goto label_90;
            case 0: goto label_120;
            default: goto label_85;
        }

        label_85:
        NPC=NPC+1;

        label_90:
        XL=arr_X[I];
        YL=arr_Y[I];

        label_100:
        continue;
    }

    NPC = NPC % 2;
    return;

    label_120:
    NPC=2;
    return;

    label_140:
    std::cout << "0  POLYGON LOGIC ERROR - PROGRAM SHOULD NOT REACH THIS POINT" << std::endl;
}

void trans4d::RADII(double const& YLAT, double& RADMER, double& RADPAR)
{
    // C
    // C  Computes the radius of curvature in the meridian
    // C  and the radius of curvature in a parallel of latitude
    // C

    DECLARE_COMMON_CONST
    double COSLAT = DCOS(YLAT);
    double DENOM = DSQRT(1.e0 + EPS*COSLAT*COSLAT);
    RADMER = AF/(pow(DENOM, 3));
    RADPAR = AF*COSLAT/DENOM;
    return;
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

void trans4d::TODMSS(double& val, int& id, int& im, double& s, int& isign)
{
    DECLARE_COMMON_CONST
    while(val > TWOPI)
    {
        val = val-TWOPI;
    }
 
    while(val < -TWOPI)
    {
        val = val + TWOPI;
    }

    if(val < 0)
        isign=-1;
    else
        isign=+1;
 
    s=DABS(val*RHOSEC/3600.0);
    id=IDINT(s);
    s=(s-id)*60.0;
    im=IDINT(s);
    s=(s-im)*60.0;
 
    // account for rounding error
 
      int is=IDINT(s*1.e5);
      if(is >= 6000000)
      {
        s=0.e0;
        im=im+1;
      }
      if(im >= 60){
        im=0;
        id=id+1;
      }
}

void trans4d::to_itrf2014(double const& x1, double const& y1, double const& z1,
double& x2, double& y2, double& z2, double& date, int const& jopt){

    //*** Converts cartesian coordinates in a specified reference
    //*** to ITRF2014 cartesian coordinates for the given date

    //*** (x1, y1, z1) --> input coordiates (meters)
    //*** (x2, y2, z2) --> output  ITRF2014 coordinates (meters)
    //*** date --> time (decimal years) to which the input & output
    //***          coordinates correspond
    //*** jopt --> input specifier of input reference frame

    DECLARE_COMMON_TRANPA

    int iopt;
    double dtime;
    double tranx, trany, tranz, rotnx, rotny, rotnz, ds;

    if (jopt == 0)
        iopt = 1;
    else
        iopt = jopt;

    dtime = date - refepc[iopt];
    tranx = -(tx[iopt] + dtx[iopt]*dtime);
    trany = -(ty[iopt] + dty[iopt]*dtime);
    tranz = -(tz[iopt] + dtz[iopt]*dtime);
    rotnx  = -(rx[iopt] + drx[iopt]*dtime);
    rotny  = -(ry[iopt] + dry[iopt]*dtime);
    rotnz  = -(rz[iopt] + drz[iopt]*dtime);
    ds     = 1.e0 - (scale[iopt] + dscale[iopt]*dtime);

    x2 = tranx + ds*x1 + rotnz*y1 - rotny*z1;
    y2 = trany - rotnz*x1 + ds*y1 + rotnx*z1;
    z2 = tranz + rotny*x1 - rotnx*y1 + ds*z1;
}

void trans4d::to_std_dev_xyz_velocity(double const& glat,double const& glon, double& sn, double& se,double& su,
double& sx, double& sy, double& sz)
{
    // *** Given standard deviations of velocity in north-south-up
    // *** Compute standard deviations of velocity in x-y-z
    // *** Assuming the covariances among north-south-up velocity
    // ***     components are zero

    double slat = DSIN(glat);
    double clat = DCOS(glat);
    double slon = DSIN(glon);
    double clon = DCOS(glon);

    sx = DSQRT( pow(slat*clon*sn,2)
              + pow(slon*se,2)
              + pow(clat*clon*su,2));

    sy = DSQRT( pow(slat*slon*sn,2)
              + pow(clon*se, 2)
              + pow(clat*slon*su, 2));

    sz = DSQRT( pow(clat*sn, 2) + pow(slat*su, 2) );
}

void trans4d::TOVNEU(double const& GLAT, double const& GLON, double& VX, double& VY, double& VZ, double& VN, double& VE, double& VU)
{
// *** Convert velocities from vx,vy,vz to vn,ve,vu

    double SLAT = DSIN(GLAT);
    double CLAT = DCOS(GLAT);
    double SLON = DSIN(GLON);
    double CLON = DCOS(GLON);

    VN = -SLAT*CLON*VX - SLAT*SLON*VY + CLAT*VZ;
    VE = -SLON*VX + CLON*VY;
    VU = CLAT*CLON*VX + CLAT*SLON*VY + SLAT*VZ;
}

void trans4d::TOVXYZ(double const& GLAT, double const& GLON, double& VN, double& VE, double& VU, double& VX, double& VY, double& VZ)
{
// *** Convert velocities from vn,ve,vu to vx,vy,vz
//       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
//       IMPLICIT INTEGER(I-N)

    double SLAT = DSIN(GLAT);
    double CLAT = DCOS(GLAT);
    double SLON = DSIN(GLON);
    double CLON = DCOS(GLON);

    VX = -SLAT*CLON*VN - SLON*VE + CLAT*CLON*VU;
    VY = -SLAT*SLON*VN + CLON*VE + CLAT*SLON*VU;
    VZ =  CLAT*VN + SLAT*VU;

    return;
}

void trans4d::TOXYZ(double glat, double glon, double eht, double& x, double& y, double& z)
{
    // *** compute x,y,z
    // *** ref p.17 geometric geodesy notes vol 1, osu, rapp
    
    //       implicit double precision(a-h,o-z)
    //       common/CONST/ a,f,e2,ep2,af,pi,twopi,rhosec
    DECLARE_COMMON_CONST
 
    double slat=DSIN(glat);
    double clat=DCOS(glat);
    double w=DSQRT(1.e0-E2*slat*slat);
    double en=A/w;
 
    x=(en+eht)*clat*DCOS(glon);
    y=(en+eht)*clat*DSIN(glon);
    z=(en*(1.e0-E2)+eht)*slat;

    return;
}

void trans4d::VTRANF(double& X, double& Y, double& Z, double& VX, double& VY, double& VZ, int IOPT1, int IOPT2)
{
    // *** Convert velocity from reference frame of IOPT1 to 
    // *** reference frame of IOPT2.

    DECLARE_COMMON_TRANPA

    double WX = 0;
    double WY = 0;
    double WZ = 0;
    double DS = 0;

    if(IOPT1 <= numref && IOPT2 <= numref && IOPT1 > 0 && IOPT2 > 0 )
    {
        // *** Convert from mm/yr to m/yr
        VX = VX /1000.e0;
        VY = VY / 1000.e0;
        VZ = VZ / 1000.e0;

        // *** From IOPT1 to ITRF2014
        // *** (following equations use approximations assuming
        // *** that rotations and scale change are small)
        WX = -drx[IOPT1];           
        WY = -dry[IOPT1];      
        WZ = -drz[IOPT1];      
        DS = -dscale[IOPT1];
        VX = VX - dtx[IOPT1] + DS*X + WZ*Y - WY*Z;
        VY = VY - dty[IOPT1] - WZ*X  +DS*Y + WX*Z;
        VZ = VZ - dtz[IOPT1] + WY*X - WX*Y + DS*Z;

        // *** From ITRF2014 to IOPT2 reference frame
        // *** (following equations use approximations assuming
        // ***  that rotations and scale change are small)
        WX = drx[IOPT2];
        WY = dry[IOPT2];
        WZ = drz[IOPT2];
        DS = dscale[IOPT2];
        VX = VX + dtx[IOPT2] + DS*X + WZ*Y - WY*Z;
        VY = VY + dty[IOPT2] - WZ*X + DS*Y + WX*Z;
        VZ = VZ + dtz[IOPT2] + WY*X - WX*Y + DS*Z;

        // *** FROM m/yr to mm/yr
        VX = VX * 1000.e0;
        VY = VY * 1000.e0;
        VZ = VZ * 1000.e0;

    }
    else
    {
        std::cout << "Improper reference frame in routine vtranf" << std::endl;
    }
}

void trans4d::XTOITRF2014(double& X, double& Y, double&Z ,double& RLAT, double& WLON, double& EHT14, double& DATE, int const& IOPT)
{
    // Converts X,Y,Z in specified datum to latitude and
    // longitude (in radians) and height (meters) in ITRF2014
    // datum with longitude positive west.

    DECLARE_COMMON_CONST

    double X1, Y1, Z1;
    double ELON;

    // Convert to cartesian coordinates in ITRF2014   
    if (IOPT == 16)
    {
        X1 = X;
        Y1 = Y;
        Z1 = Z;
    }
    else   
    {
        to_itrf2014(X,Y,Z,X1,Y1,Z1,DATE,IOPT);
    }
      
    // Convert to geodetic coordinates
    if(!FRMXYZ(X1,Y1,Z1,RLAT,ELON,EHT14))
        //C++ port, don't want to kill the app in the library code STOP(666);
        std::cout << "FRMXYZ failed! X:" << X << " Y:" << Y << " Z:" << Z << std::endl;

    WLON = -ELON;
    while(WLON < 0){
        WLON = WLON + TWOPI;
    }
}

