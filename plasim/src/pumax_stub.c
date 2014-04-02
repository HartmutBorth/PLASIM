#include <sys/resource.h>

void pumax_dummy(void) {}

int nresources_(double *ut, double *st, long *mem, long *par, long *paf,
               long *swa, long *dr, long *dw)
{
   struct rusage ru;
   getrusage(RUSAGE_SELF,&ru);
   *ut = ru.ru_utime.tv_sec + 0.000001 * ru.ru_utime.tv_usec;
   *st = ru.ru_stime.tv_sec + 0.000001 * ru.ru_stime.tv_usec;
   *mem = ru.ru_maxrss;
   *par = ru.ru_minflt;
   *paf = ru.ru_majflt;
   *swa = ru.ru_nswap;
   *dr  = ru.ru_inblock;
   *dw  = ru.ru_oublock;
   return 1;
}

/* ------------------------------------------------ */
/* Stub routines for Absoft Compiler and others,    */
/* which require, that FORTRAN callable C-functions */
/* are written in uppercase letters only            */
/* ------------------------------------------------ */


int NRESOURCES(double *ut, double *st, long *mem, long *par, long *paf,
               long *swa, long *dr, long *dw)
{
   return nresources_(ut,st,mem,par,paf,swa,dr,dw);
}

/* ------------------------------------------------ */
/* Stub routines for IBM Compiler and others,       */
/* which require, that FORTRAN callable C-functions */
/* are written in lowercase without underscore      */
/* ------------------------------------------------ */


int nresources(double *ut, double *st, long *mem, long *par, long *paf,
              long *swa, long *dr, long *dw)
{
   return nresources_(ut,st,mem,par,paf,swa,dr,dw);
}
