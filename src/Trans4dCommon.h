/*PORTED COMMON BLOCKS */

#pragma region BNDRY
#define DECLARE_COMMON_BNDRY \
/* BNDRY */\
double (&X)[5000 + 1] = common_bndry.X;\
double (&Y)[5000 + 1] = common_bndry.Y;\
int (&NPOINT)[30 + 1] = common_bndry.NPOINT;
extern struct BNDRY
{
    double X[5000 + 1] = {0};
    double Y[5000 + 1] = {0};
    int NPOINT[30 + 1] = {0};
} common_bndry;
#pragma endregion


#pragma region CDGRID
#define NUMGRD 8
#define DECLARE_COMMON_CDGRID\
/* CDGRID */\
double (&grdlx)[NUMGRD+1] = common_cdgrid.grdlx;\
double (&grdux)[NUMGRD+1] = common_cdgrid.grdux;\
double (&grdly)[NUMGRD+1] = common_cdgrid.grdly;\
double (&grduy)[NUMGRD+1] = common_cdgrid.grduy;\
int (&icntx)[NUMGRD+1] = common_cdgrid.icntx;\
int (&icnty)[NUMGRD+1] = common_cdgrid.icnty;\
int (&nbase)[NUMGRD+1] = common_cdgrid.nbase;
extern struct CDGRID
{
  double grdlx[NUMGRD+1] = {0};
  double grdux[NUMGRD+1] = {0};
  double grdly[NUMGRD+1] = {0};
  double grduy[NUMGRD+1] = {0};
  int icntx[NUMGRD+1] = {0};
  int icnty[NUMGRD+1] = {0};
  int nbase[NUMGRD+1] = {0};
} common_cdgrid;
#pragma endregion


#pragma region CONST
#define DECLARE_COMMON_CONST \
/* CONST */\
double &A = common_const.A;\
double &F = common_const.F;\
double &E2 = common_const.E2;\
double &AF = common_const.AF;\
double &EPS = common_const.EPS;\
double &PI = common_const.PI;\
double &RHOSEC = common_const.RHOSEC;\
double &TWOPI = common_const.TWOPI;
extern struct CONST
{
    double A;
    double F;
    double E2;
    double AF;
    double EPS;
    double PI;
    double RHOSEC;
    double TWOPI;
} common_const;
#pragma endregion


#pragma region PSGRID
#define DECLARE_COMMON_PSGRID \
/* PSGRID */\
double (&PSGLX)[PSGRID::numpsg +1] = common_psgrid.PSGLX;\
double (&PSGUX)[PSGRID::numpsg +1] = common_psgrid.PSGUX;\
double (&PSGLY)[PSGRID::numpsg +1] = common_psgrid.PSGLY;\
double (&PSGUY)[PSGRID::numpsg +1] = common_psgrid.PSGUY;\
int (&ICNTPX)[PSGRID::numpsg  +1] = common_psgrid.ICNTPX;\
int (&ICNTPY)[PSGRID::numpsg +1] = common_psgrid.ICNTPY;\
int (&NBASEP)[PSGRID::numpsg  +1] = common_psgrid.NBASEP;
extern struct PSGRID
{
    static const int numpsg = 1;
    double PSGLX[numpsg +1] = {0};
    double PSGUX[numpsg +1] = {0};
    double PSGLY[numpsg+1] = {0};
    double PSGUY[numpsg +1] = {0};
    
    int ICNTPX[numpsg +1] = {0};
    int ICNTPY[numpsg+1] = {0};
    int NBASEP[numpsg +1] = {0};
} common_psgrid;
#pragma endregion


#pragma region PGRID
#define DECLARE_COMMON_PGRID \
/* PGRID */\
double (&PS)[18000 + 1] = common_pgrid.PS;
extern struct PGRID
{
    double PS[18000 + 1] = {0};
} common_pgrid;
#pragma endregion


#pragma region QPARM
#define DECLARE_COMMON_QPARM \
/* QPARM */\
double (&STRIKE)[QPARM::ndloc +1] = common_qparm.STRIKE;\
double (&HL)[QPARM::ndloc +1] = common_qparm.HL;\
double (&EQLAT)[QPARM::ndloc +1] = common_qparm.EQLAT;\
double (&EQLON)[QPARM::ndloc +1] = common_qparm.EQLON;\
double (&SSLIP)[QPARM::ndloc +1] = common_qparm.SSLIP;\
double (&DSLIP)[QPARM::ndloc +1] = common_qparm.DSLIP;\
double (&DIP)[QPARM::ndloc +1] = common_qparm.DIP;\
double (&DEPTH)[QPARM::ndloc +1] = common_qparm.DEPTH;\
double (&WIDTH)[QPARM::ndloc +1] = common_qparm.WIDTH;\
double (&EQLATR)[50+1] = common_qparm.EQLATR;\
double (&EQLONR)[50+1] = common_qparm.EQLONR;\
double (&EQRAD)[50+1] = common_qparm.EQRAD;\
int (&ITEQK)[50 +1] = common_qparm.ITEQK;\
int (&NLOC)[50+1] = common_qparm.NLOC;\
int (&NFP)[50 +1] = common_qparm.NFP;\
int (&NUMEQ) = common_qparm.NUMEQ;
extern struct QPARM
{
    static const int ndloc = 2195;
    double STRIKE[ndloc +1] = {0};
    double HL[ndloc +1] = {0};
    double EQLAT[ndloc+1] = {0};
    double EQLON[ndloc +1] = {0};
    double SSLIP[ndloc +1] = {0};
    double DSLIP[ndloc +1] = {0};
    double DIP[ndloc +1] = {0};
    double DEPTH[ndloc +1] = {0};
    double WIDTH[ndloc +1] = {0};
    double EQLATR[50 +1] = {0};
    double EQLONR[50 +1] = {0};
    double EQRAD[50 +1] = {0};
    int ITEQK[50 +1] = {0};
    int NLOC[50+1] = {0};
    int NFP[50 +1] = {0};
    int NUMEQ = 0;
} common_qparm;
#pragma endregion


#pragma region REFCON
#define DECLARE_COMMON_REFCON\
/* REFCON */\
int (&IRFCON)[29 + 1] = common_refcon.IRFCON;\
int (&JRFCON)[REFCON::numref + 1] = common_refcon.JRFCON;
extern struct REFCON
{
    static const int numref = 17;
    int IRFCON[29 + 1] =  {0};
    int JRFCON[numref + 1] = {0};
} common_refcon;
#pragma endregion


#pragma region TIMREF
#define DECLARE_COMMON_TIMREF \
/* TIMREF */\
int ITREF = common_timref.ITREF;
extern struct TIMREF
{
    int ITREF;
} common_timref;
#pragma endregion


#pragma region TRANPA
#define DECLARE_COMMON_TRANPA \
/*TRANPA*/ \
double (&tx)[TRANPA::numref + 1] = common_tranpa.tx;\
double (&ty)[TRANPA::numref + 1] = common_tranpa.ty;\
double (&tz)[TRANPA::numref + 1] = common_tranpa.tz;\
double (&dtx)[TRANPA::numref + 1] = common_tranpa.dtx;\
double (&dty)[TRANPA::numref + 1] = common_tranpa.dty;\
double (&dtz)[TRANPA::numref + 1] = common_tranpa.dtz;\
double (&rx)[TRANPA::numref + 1] = common_tranpa.rx;\
double (&ry)[TRANPA::numref + 1] = common_tranpa.ry;\
double (&rz)[TRANPA::numref + 1] = common_tranpa.rz;\
double (&drx)[TRANPA::numref + 1] = common_tranpa.drx;\
double (&dry)[TRANPA::numref + 1] = common_tranpa.dry;\
double (&drz)[TRANPA::numref + 1] = common_tranpa.drz;\
double (&scale)[TRANPA::numref + 1] = common_tranpa.scale;\
double (&dscale)[TRANPA::numref + 1] = common_tranpa.dscale;\
double (&refepc)[TRANPA::numref + 1] = common_tranpa.refepc;
extern struct TRANPA
{
    static const int numref = 17;
    double tx[numref + 1] = {0};
    double ty[numref + 1] = {0};
    double tz[numref + 1] = {0};
    double dtx[numref + 1] = {0};
    double dty[numref + 1] = {0};
    double dtz[numref + 1] = {0};
    double rx[numref + 1] = {0};
    double ry[numref + 1] = {0};
    double rz[numref + 1] = {0};
    double drx[numref + 1] = {0};
    double dry[numref + 1] = {0};
    double drz[numref + 1] = {0};
    double scale[numref + 1] = {0};
    double dscale[numref + 1] = {0};
    double refepc[numref + 1] = {0};
} common_tranpa;
#pragma endregion


#pragma region VGRID
#define DECLARE_COMMON_VGRID \
double (&b)[800000 + 1] = common_vgrid.b;
/* VGRID */
extern struct VGRID
{
  double b[800000 +1] = {0};
} common_vgrid;
#pragma endregion
