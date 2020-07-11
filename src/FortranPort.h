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