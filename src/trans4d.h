#include "utility_helpers.h"
#include "Trans4dCommon.h"

//You must change Trans4D version here if necessary
#define TRANS4D_VERSION "0.2.6"

namespace trans4d
{
   
    void InitBd();
    void InitEq();
    void InitPs();
    void InitVl();

    void InitBlockData();



    //ORIGINAL HTDP FUNCTIONS with changes for TRANS4D
    //C style Function Prototypes for Fortran Subroutines to be ported
    //ordered A-Z
    void COMPSN(double& YLATT, double& YLONT, double& HTT, double const& YLAT, double const& YLON, double const& HT, int const& MIN, double const& VN, double const& VE, double const& VU);
    void COMVEL(double& YLAT, double& YLON, int& JREGN, double& VN, double& VE, double& VU, double& SN, double& SE, double&SU);
    void DISLOC(double const& YLAT, double const& YLON, double const& STRIKE, double const& HL, double const& EQLAT,
        double const& EQLON, double const& SS, double const& DS, double const& DIP, double const& DEPTH, double const& WIDTH,
        double& DNORTH, double& DWEST, double& DUP);
    void FRIT94(double x1, double y1, double z1, double& x2, double& y2, double& z2, double& date, int jopt);
    void FRIT94_IERS(double x1, double y1, double z1, double& x2, double& y2, double& z2, double& date, int jopt);
    bool FRMXYZ(double& x, double& y, double& z, double& glat, double& glon, double& eht);
    void GETBDY();
    void GETREG(double& X0, double& YKEEP, int& JREGN);
    void GRDAMP(int const& K, int const& I, int const& J, double (&AMP)[2][2][3], double const (&PS)[18000+1]);
    void GRDCHK(double const& POSX, double const& POSY, double const& GRDLX, double const& GRDUX, double const& GRDLY, double const& GRDUY, bool& INSIDE);
    void GRDVEC(int JREGN, int I, int J, double (&VEL)[2][2][3], double (&B)[210000+1]);
    void GRDWEI(double const& YLON, double const& YLAT, int const& JREGN, int& I, int&J, double (&WEI)[2][2] );
    void GTOVEL(double const& YLAT, double const& YLON,  double const& EHT,
        double& VN, double& VE, double& VU, double& VX, double& VY, double& VZ, int& JREGN, int const& IOPT,
        double& SN, double& SE, double& SU, double& SX, double& SY, double& SZ);
    int IPSGRD(int const& IGRID, int const& I, int const& J, int const& IVEC);
    int IUNGRD(int IREGN, int I, int J, int IVEC);
    //C++ Port: Moved IYMDMJ to utility_helpers.h
    //void IYMDMJ(int& IYR, int& IMON, int& IDAY, int& MJD);
    void MODEL();
    void NEWCOR(double& YLAT, double& YLON, double& HTOLD, int& MIN1, int& MIN2, double& YLAT3, double& YLON3, double& HTNEW, 
        double& DN, double& DE, double& DU, double& VN, double& VE, double& VU);
    void OKADA(double const& X1, double const& X2, double const& XL, double const& DU, double const& W, double const& DIP,
        double& U1SS, double& U2SS, double& U3SS, double& U1DS, double& U2DS, double& U3DS);
    void OKADAW(double const& PSI, double const& ETA, double const& Q, double const& SDIP, double const& CDIP, double const& RATIO,
        double const& TWOPI, bool const& VERT, double& U1SS, double& U2SS, double& U3SS, double& U1DS, double& U2DS, double& U3DS);
    void PLATVL(int& IPLATE, double& X, double& Y, double& Z, double& VX, double& VY, double& VZ);
    void POLYIN(double& X0, double& Y0, double& X, double& Y, int& N, int& NPC);
    void PREDV(double& ylat, double& ylon, double& eht, double& date, int& iopt, int& jregn, double& vn, double& ve, double& vu);
    void PSDISP(double const& YLAT, double const& YLON, int const& MIN, double& DNORTH, double& DEAST, double& DUP);
    void PSGWEI(double const& POSX, double const& POSY, int const& K, int& I, int& J, double (&WEI)[2][2]);
    void RADII(double const& YLAT, double& RADMER, double& RADPAR);
    void RADR8T(double const& YLAT, double const& VN, double const& VE, double& VNR, double& VER);
    void SETRF();
    void SETTP();
    void TODMSS(double& val, int& id, int& im, double& s, int& isign);
    void TOIT94(double x1, double y1, double z1, double& x2, double& y2, double& z2, double date, int& jopt);
    void TOIT94_IERS(double& x1, double& y1, double& z1, double x2, double y2, double z2, double date, int& jopt);
    void to_itrf2014(double const& x1, double const& y1, double const& z1,
        double& x2, double& y2, double& z2, double& date, int const& jopt);
    void TOMNT(int& IYR, int& IMON, int& IDAY, int& IHR, int& IMN, int& MINS);
    void TOVNEU(double& GLAT, double& GLON, double& VX, double& VY, double& VZ, double& VN, double& VE, double& VU);
    void TOVXYZ(double& GLAT, double& GLON, double& VN, double& VE, double& VU, double& VX, double& VY, double& VZ);
    void TOXYZ(double glat, double glon, double eht, double& x, double& y, double& z);
    void VTRANF(double& X, double& Y, double& Z, double& VX, double& VY, double& VZ, int IOPT1, int IOPT2);
    void VTRANF_IERS(double& X, double& Y, double& Z, double& VX, double& VY, double& VZ, int IOPT1, int IOPT2);
    void XTO08(double& X, double& Y, double& Z, double& RLAT, double& RLON, double& EHTNAD, double& DATE, int& IOPT);
    void XTOITRF2014(double& X, double& Y, double&Z ,double& RLAT, double& WLON, double& EHT14, double& DATE, int const& IOPT);

};