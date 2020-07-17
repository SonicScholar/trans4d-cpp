#include "Trans4dCommon.h"

class trans4d{
    private:
        bool _blockDataInitialized;

    private:    
        void InitBd();
        void InitEq();
        void InitPs();
        void InitVl();

    public:
        trans4d():
        _blockDataInitialized(false)
        { }

        void InitBlockData();


        //todo: sort these out for trans4d
        //ORIGINAL HTDP FUNCTIONS
        //C style Function Prototypes for Fortran Subroutines to be ported
        //ordered A-Z
        void COMPSN(double& YLATT, double& YLONT, double& HTT, double const& YLAT, double const& YLON, double const& HT, int const& MIN, double const& VN, double const& VE, double const& VU);
        void COMVEL(double& YLAT, double& YLON, int& JREGN, double& VN, double& VE, double& VU);
        void DISLOC(double const& YLAT, double const& YLON, double const& STRIKE, double const& HL, double const& EQLAT,
            double const& EQLON, double const& SS, double const& DS, double const& DIP, double const& DEPTH, double const& WIDTH,
            double& DNORTH, double& DWEST, double& DUP);
        void FRIT94(double x1, double y1, double z1, double& x2, double& y2, double& z2, double& date, int jopt);
        void FRIT94_IERS(double x1, double y1, double z1, double& x2, double& y2, double& z2, double& date, int jopt);
        bool FRMXYZ(double& x, double& y, double& z, double& glat, double& glon, double& eht);
        void GETBDY();
        void GRDAMP(int const& K, int const& I, int const& J, double (&AMP)[2][2][3], double const (&PS)[18000+1]);
        void GRDCHK(double const& POSX, double const& POSY, double const& GRDLX, double const& GRDUX, double const& GRDLY, double const& GRDUY, bool& INSIDE);
        void GRDVEC(int JREGN, int I, int J, double (&VEL)[2][2][3], double (&B)[210000+1]);
        void GRDWEI(double const& YLON, double const& YLAT, int const& JREGN, int& I, int&J, double (&WEI)[2][2] );
        int IPSGRD(int const& IGRID, int const& I, int const& J, int const& IVEC);
        int IUNGRD(int IREGN, int I, int J, int IVEC);
        void IYMDMJ(int& IYR, int& IMON, int& IDAY, int& MJD);
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
        void TOIT94(double x1, double y1, double z1, double& x2, double& y2, double& z2, double date, int& jopt);
        void TOIT94_IERS(double& x1, double& y1, double& z1, double x2, double y2, double z2, double date, int& jopt);
        void TOMNT(int& IYR, int& IMON, int& IDAY, int& IHR, int& IMN, int& MINS);
        void TOVNEU(double& GLAT, double& GLON, double& VX, double& VY, double& VZ, double& VN, double& VE, double& VU);
        void TOVXYZ(double& GLAT, double& GLON, double& VN, double& VE, double& VU, double& VX, double& VY, double& VZ);
        void TOXYZ(double glat, double glon, double eht, double& x, double& y, double& z);
        void VTRANF(double& X, double& Y, double& Z, double& VX, double& VY, double& VZ, int IOPT1, int IOPT2);
        void VTRANF_IERS(double& X, double& Y, double& Z, double& VX, double& VY, double& VZ, int IOPT1, int IOPT2);
        void XTO08(double& X, double& Y, double& Z, double& RLAT, double& RLON, double& EHTNAD, double& DATE, int& IOPT);

};