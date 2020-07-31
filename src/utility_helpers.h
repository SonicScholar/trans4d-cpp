#ifndef UTILITY_HELPERS_H
#define UTILITY_HELPERS_H
#include <math.h>

/* BUILT IN FORTRAN FUNCTIONS */
//http://userweb.eng.gla.ac.uk/peter.smart/com/com/f77-conv.htm
//http://math.hawaii.edu/~dale/190/intrinsic_fortran_functions.html
/* abs */
inline double DABS(double x){ return abs(x); }
/* arctan */
inline double DATAN(double x){ return atan(x); }
/* arctan 2*/
inline double DATAN2(double y, double x){ return atan2(y,x); }
/* convert to double */
inline double DBLE(int x) { return (double)x; }
/* cos */
inline double DCOS(double x){ return cos(x); }
/* Exponential */
inline double DEXP(double x){ return exp(x); }
/* Truncation of Real */
inline double DINT(double x) { return (double)((int)x); }
/* Natural Log */
inline double DLOG(double x) { return log(x); }
/* sin */
inline double DSIN(double x){ return sin(x); }
/* square root */
inline double DSQRT(double x){ return sqrt(x); }
/* convert to int */
inline int IDINT(double x) { return (int)x; };

/* IF_ARITHMETIC */
inline int IF_ARITHMETIC(double x)
{
    if(x < 0) return -1;
    if(x == 0) return 0;
    return 1;
}

/* Utility functions ported from HTDP and TRANS4D */
inline void IYMDMJ(int IYR, int IMON, int IDAY, int& MJD)
{
// C NAME:       IYMDMJ
// C VERSION:    Sep. 17, 2010
// C WRITTEN BY: R. SNAY (after M. SCHENEWERK)
// C PURPOSE:    CONVERT DATE TO MODIFIED JULIAN DATE 
// C
// C INPUT PARAMETERS FROM THE ARGUEMENT LIST:
// C -----------------------------------------
// C IDAY              DAY
// C IMON              MONTH
// C IYR               YEAR
// C
// C OUTPUT PARAMETERS FROM ARGUEMENT LIST:
// C --------------------------------------
// C MJD               MODIFIED JULIAN DATE 
// C
// C
// C LOCAL VARIABLES AND CONSTANTS:
// C ------------------------------
// C A                 TEMPORARY STORAGE
// C B                 TEMPORARY STORAGE
// C C                 TEMPORARY STORAGE
// C D                 TEMPORARY STORAGE
// C IMOP              TEMPORARY STORAGE
// C IYRP              TEMPORARY STORAGE
// C
// C GLOBAL VARIABLES AND CONSTANTS:
// C ------------------------------
// C
// C
// C       THIS MODULE CALLED BY: GENERAL USE
// C
// C       THIS MODULE CALLS:     DINT
// C
// C       INCLUDE FILES USED:
// C
// C       COMMON BLOCKS USED:       
// C
// C       REFERENCES:            DUFFETT-SMITH, PETER  1982, 'PRACTICAL
// C                              ASTRONOMY WITH YOUR CALCULATOR', 2ND
// C                              EDITION, CAMBRIDGE UNIVERSITY PRESS,
// C                              NEW YORK, P.9
// C
// C       COMMENTS:              THIS SUBROUTINE REQUIRES THE FULL YEAR,
// C                              I.E. 1992 RATHER THAN 92.  
// C
// C********1*********2*********3*********4*********5*********6*********7**
// C::LAST MODIFICATION
// C::8909.06, MSS, DOC STANDARD IMPLIMENTED
// C::9004.17, MSS, CHANGE ORDER YY MM DD
// C********1*********2*********3*********4*********5*********6*********7**
  int A, B, C, D;
  int IYRP = IYR;

  int IMOP;//
  if (IMON < 3) {
    IYRP = IYRP - 1;
    IMOP = IMON + 12;
  }  else {
    IMOP = IMON;
  }
// C
// C........  1.0  CALCULATION
// C
      A=  IYRP*0.01;
      B=  2 - A + DINT( A*0.25 );
      C=  365.25*IYRP;
      D=  30.6001*(IMOP + 1);
      MJD =  (B + C + D + IDAY - 679006);
// C      
    return;
}

inline void DecimalYearToMJD(double date, int& MJD)
{
    int M[12 +1];
    M[1] = 31;
    M[2] = 28;
    M[3] = 31;
    M[4] = 30;
    M[5] = 31;
    M[6] = 30;
    M[7] = 31;
    M[8] = 31;
    M[9] = 30;
    M[10] = 31;
    M[11] = 30;
    M[12] = 31;

    int IYEAR = date;
    int MJD0;
    IYMDMJ(IYEAR, 1, 1, MJD0);
    int IYEAR1 = IYEAR + 1;
    int MJD1;
    IYMDMJ(IYEAR1, 1, 1, MJD1);
    int LEAPDAY = (MJD1 - MJD0) == 366 ? 1 : 0;
    
    double REMDAY = (date - IYEAR)* (MJD1 - MJD0);
    int IBEGIN = 0;
    int ITOTAL = 31;


    if (REMDAY < ITOTAL) 
    {
        int MONTH = 1;
        int IDAY = REMDAY - IBEGIN + 1;
        IYMDMJ(IYEAR,MONTH,IDAY,MJD);
        return;
    }
    IBEGIN = ITOTAL;
    ITOTAL = ITOTAL + LEAPDAY;
    for(int i = 2; i <= 12; i++)
    {
        ITOTAL = ITOTAL + M[i];
        if (REMDAY < ITOTAL) 
        {
            int MONTH = i;
            int IDAY = REMDAY - IBEGIN + 1;
            IYMDMJ(IYEAR,MONTH,IDAY,MJD);
            return;
        }
        IBEGIN = ITOTAL;
    }
}

inline void DecimalYearToMJDMins(double date, int& MINS)
{
    int MJD;
    DecimalYearToMJD(date, MJD);
    MINS = MJD * 24 * 60;
}
#endif