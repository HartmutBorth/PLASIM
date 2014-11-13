#include <valarray>

// #define MPI_INTERFACE

using namespace std;

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// Definitions

#define NROOT 0

// Switches

int PrintLats = 1;
int PrintPoly = 1;
int PrintSigs = 1;

// Dimensions

int Lons; // longitudes
int Lats; // latitudes
int Levs; // levels

int Hors; // horizontal gridpoints (Lons * Lats)

int Trit; // triangular truncation
int Spec; // spectral coefficients (Trit+1) * (Trit+2)

int Lapp; // latitudes per process
int Hopp; // horizontal gridpoints per process
int Happ; // half latitudes per process

int ncsp = 0; // # of complex  spectral modes
int nesp = 0; // # of extented spectral modes
int nspp = 0; // # of spectral modes per process
int Procs = 1; // # of processes

// Integer

int nits  =  3; // # of initial time steps
int ntspd = 24; // # of time steps per day

// Double

#define PLARAD 6371000.0
#define WW 0.00007292
#define CV (PLARAD * WW)
#define GASCON 287.0
#define CT ((CV*CV)/GASCON)

double StepDelta           ; // length of time step (2 Pi / time steps per day)
double dtns   =    20.0    ;
double dtep   =    60.0    ;
double t0     =   250.0    ;
double tgr    =   288.0    ; // Ground Temperature in mean profile [K]
double dtrop  = 12000.0    ; // Tropopause height [m]
double dttrp  =     2.0    ; // Tropopause smoothing [K]
double GA     =     9.81   ; // Gravity
double ALR    =     0.0065 ; // Lapse rate
double plavor = 1.632993161855452; // 1 / sqrt(3/8) = planetary vorticity

// Latitudinal arrays

double *gsi ; // Gaussian abscissas (sin(lat))
double *gaw ; // Gaussian weights

// Level arrays

double *SigHalf; // Half level 
double *SigFull; // Full level
double *SigDelt; // Delta sigma
double *SigRevd; // 1.0 / (2 * delta sigma)

// Gridpoint array slices

valarray<double> dpsdmu; // d(LnPs)/d(mu)

// Some functions for a nice printout

#define COLS 72

#ifndef MIN
#define MIN(x,y) ((x)<(y)?(x):(y))
#endif

#ifndef MAX
#define MAX(x,y) ((x)<(y)?(y):(x))
#endif


void Stars(int n)  {while (n--) putchar('*');}
void ErrStars(int n)  {while (n--) putc('*',stderr);}
void Blanks(int n) {while (n--) putchar(' ');}
void Dashes(int n) {putchar('*');putchar(' ');while (--n>3) putchar('-');putchar(' ');putchar('*');}
void NewLine(void) {putchar('\n');}
void ErrNewLine(void) {putc('\n',stderr);}
void StarLine(void) {Stars(COLS); NewLine();}

/* ==================================== */
/* Abort - Print error message and exit */
/* ==================================== */

void Abort(const char *errtext)
{
   int l = strlen(errtext);
   Stars(MIN(COLS,l)); NewLine();
   puts(errtext);
   Stars(MIN(COLS,l)); NewLine();
   ErrStars(MIN(COLS,l))  ; ErrNewLine();
   fputs(errtext,stderr)  ; ErrNewLine();
   ErrStars(MIN(COLS,l))  ; ErrNewLine();
   exit(1);
}


void LeftText(const char *t, int cols)
{
   int l;
   l = strlen(t);
   if (l < 1) return;
   if (l > cols-4) puts(t);
   else 
   {
      putchar('*');
      putchar(' ');
      fputs(t,stdout);
      Blanks(cols - l - 3);
      putchar('*');
      NewLine();
   }
}


// -------------------------- MPI -------------------------------------


int mypid = NROOT;

#ifdef MPI_INTERFACE
#include <mpi.h>

void mpbrdn(double *f, int n)
{
   MPI_Bcast(f,n,MPI_DOUBLE,NROOT,MPI_COMM_WORLD);
}

void mpgadn(double *f, double *p, int n)
{
     MPI_Gather(p,n,MPI_DOUBLE,f,n,MPI_DOUBLE,NROOT,MPI_COMM_WORLD);
}

void mpscdn(double *f, double *p, int n)
{
   MPI_Scatter(f,n,MPI_DOUBLE,p,n,MPI_DOUBLE,NROOT,MPI_COMM_WORLD);
}

void mpsush(double *f, int n)
{
   double p[n];
   MPI_Reduce(f,p,n,MPI_DOUBLE,MPI_SUM,NROOT,MPI_COMM_WORLD);
   if (mypid == NROOT) memcpy(f,p,sizeof(p));
}

void mp_start(int *argc, char ***argv)
{
   MPI_Status status;
   MPI_Init(argc,argv);
   MPI_Comm_size(MPI_COMM_WORLD,&Procs );
   MPI_Comm_rank(MPI_COMM_WORLD,&mypid);

   if (Procs > 1)
   {
      Lapp = Lats / Procs;
      Happ = Lapp / 2;
      Hopp = Lons * Lapp;
      nspp = Spec / Procs;
      if (Spec % Procs) nspp++;
      nesp = nspp * Procs;
   }

   NewLine();
   if (Procs > 1)
   {
      int len;
      char Pname[Procs*256];
      char tb[80];
      MPI_Get_processor_name(Pname,&len);
      MPI_Gather(Pname,256,MPI_CHAR,Pname,256,MPI_CHARACTER,NROOT,MPI_COMM_WORLD);
      MPI_Allreduce(&len,&len,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
      if (mypid == NROOT)
      {
         Stars(len+18);
         NewLine();
         for (int i=0 ; i<Procs ; ++i)
         {
            sprintf(tb,"Process %4d: ",i);
            strncat(tb,Pname+256*i,len);
            LeftText(tb,len+18);
         }
         Stars(len+18);
         NewLine();
         NewLine();
      }
   }
}


void mp_stop(void)
{
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
}
#else
void mp_start(int *argc, char ***argv) {}
void mp_stop(void) {}
void mpbrdn(double *f, int n) {}
void mpsush(double *f, int n) {}
void mpgadn(double *f, double *p, int n)
{
   memcpy(f,p, n * sizeof(double));
}
void mpscdn(double *f, double *p, int n)
{
   memcpy(p,f, n * sizeof(double));
}



#endif


// -------------------------- PC_Var ----------------------------------

class PC_Var
{
   public:
      void Init(int code, int levs, const char *vn, const char *na);
      void Status(FILE *fp, const char *mess);
      void StatusGP(FILE *fp, const char *mess);
      void StatusFC(FILE *fp, const char *mess);
      void StatusSP(FILE *fp, const char *mess);
      void Alt2Reg(double *);
      void Reg2Alt(double *);
      void Alt2CS(double *);
      void ReadGP(void);
      void WriteGP(const char *);
      void WriteFC(const char *);
      void WriteSP(void);
      void WriteServiceGP(const char *);
      void gp2fc(void);
      void fc2gp(void);
      void sh2fc(void);
      void sh2fcdmu(double *dshdmu);
      void Section(void);
      void BroadcastSpectral();

      valarray<double>gp; // grid points
      valarray<double>fc; // fourier coefficients
      valarray<double>sh; // spherical  harmonics
      valarray<double>cs; // zonal mean cross section

   private:
      void StatusLineGP(FILE *fp, int i);
      void StatusLineFC(FILE *fp, int i);
      void StatusLineSP(FILE *fp, int i);
      void FFT_ini(void);
      void FFT_d2(double *, int);
      void FFT_i2(double *, int);
      void FFT_e2(double *);
      void FFT_s2(double *);
      void FFT_d8(double *);
      void FFT_i8(double *);
      int Code;
      int VarLevs;
      int DimG;

      long Head[ 8];
      char Varn[ 8];
      char Name[72];

      double *FFT_cos;
      double *FFT_sin;
};

void PC_Var::Init(int code, int levs, const char *vn, const char *na)
{
   Code    = code;
   VarLevs = levs;
   DimG = Hopp * VarLevs;

   strncpy(Varn,vn, 7);
   strncpy(Name,na,71);

   gp.resize(DimG,0.0);           // Allocate grid point array
   fc.resize(DimG,0.0);           // Allocate fourier coefficient array
   sh.resize(nesp*VarLevs,0.0);   // Allocate spherical harmonics array
   cs.resize(Lats*VarLevs,0.0);   // Allocate cross section array

   Head[0] = Code;
   Head[4] = Lons;
   Head[5] = Lats;

   FFT_ini();
}


void PC_Var::BroadcastSpectral(void)
{
   mpbrdn(&sh[0],nesp * VarLevs);
}


void PC_Var::Status(FILE *fp, const char *mess)
{
   int l = MAX(strlen(Varn) + strlen(Name) + strlen(mess) + 9, 16);

   fprintf(fp,"\n");
   fprintf(fp,"%2d:%s : %s [%s]\n",mypid,Varn,Name,mess);
   for (int i=0 ; i<l ; ++i) putc('=',fp); putc('\n',fp);
   fprintf(fp,"%2d:Lons =%10d\n",mypid,Lons);
   fprintf(fp,"%2d:Lapp =%10d\n",mypid,Lapp);
   fprintf(fp,"%2d:Levs =%10d\n",mypid,VarLevs);
   fprintf(fp,"%2d:Trit =%10d\n",mypid,Trit);
   fprintf(fp,"%2d:Spec =%10d\n",mypid,Spec);
   fprintf(fp,"%2d:Hopp =%10d\n",mypid,Hopp);
   fprintf(fp,"%2d:DimG =%10d\n",mypid,DimG);
   if (gp.size() > 0)
   {   
      fprintf(fp,"%2d:Min  =%10.2lf\n",mypid,gp.min());
      fprintf(fp,"%2d:Max  =%10.2lf\n",mypid,gp.max());
   }
   fprintf(fp,"\n");
}


void PC_Var::StatusLineGP(FILE *fp, int i)
{
   if (Procs > 1) fprintf(fp,"%2d:",mypid);
   fprintf(fp,"gp[%2d][0:5] = %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf\n",
           i/Lons,gp[i],gp[i+1],gp[i+2],gp[i+3],gp[i+4],gp[i+5]);
}


void PC_Var::StatusGP(FILE *fp, const char *mess)
{
   int i;
   int l = MAX(strlen(Varn) + strlen(Name) + strlen(mess) + 9, 16);

   fprintf(fp,"\n");
   fprintf(fp,"%2d:%s : %s [%s]\n",mypid,Varn,Name,mess);
   for (i=0 ; i<l ; ++i) putc('=',fp); putc('\n',fp);

   StatusLineGP(fp,          0);
   StatusLineGP(fp,     Lons  );
   StatusLineGP(fp,Hopp-Lons*2);
   StatusLineGP(fp,Hopp-Lons  );
}


void PC_Var::StatusLineFC(FILE *fp, int i)
{
   if (Procs > 1) fprintf(fp,"%2d:",mypid);
   fprintf(fp,"fc[%2d][0:5] = %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf\n",
           i/Lons,fc[i],fc[i+1],fc[i+2],fc[i+3],fc[i+4],fc[i+5]);
}


void PC_Var::StatusFC(FILE *fp, const char *mess)
{
   int i;
   int l = MAX(strlen(Varn) + strlen(Name) + strlen(mess) + 9, 16);

   fprintf(fp,"\n");
   fprintf(fp,"%2d:%s : %s [%s]\n",mypid,Varn,Name,mess);
   for (i=0 ; i<l ; ++i) putc('=',fp); putc('\n',fp);

   StatusLineFC(fp,          0);
   StatusLineFC(fp,     Lons  );
   StatusLineFC(fp,Hopp-Lons*2);
   StatusLineFC(fp,Hopp-Lons  );
}


void PC_Var::StatusLineSP(FILE *fp, int i)
{
   if (Procs > 1) fprintf(fp,"%2d:",mypid);
   fprintf(fp,"sh[%7d] = %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf\n",
           i,sh[i],sh[i+1],sh[i+2],sh[i+3],sh[i+4],sh[i+5]);
}


void PC_Var::StatusSP(FILE *fp, const char *mess)
{
   int i;
   int l = MAX(strlen(Varn) + strlen(Name) + strlen(mess) + 9, 16);

   fprintf(fp,"\n");
   fprintf(fp,"%2d:%s : %s [%s]\n",mypid,Varn,Name,mess);
   for (i=0 ; i<l ; ++i) putc('=',fp); putc('\n',fp);
   StatusLineSP(fp,          0);
   StatusLineSP(fp,i  = 2 * (Trit+1));
   StatusLineSP(fp,i += 2 *  Trit   );
   StatusLineSP(fp,i += 2 * (Trit-1));
}


void PC_Var::Alt2Reg(double *alt)
{
   int l,n,s,e,o;
   double reg[Hors];

   l = Lons * sizeof(double);

   for (e=0,o=Lons,n=0,s=Hors-Lons ; e<Hors ; e+=2*Lons,o+=2*Lons,n+=Lons,s-=Lons)
   {
      memcpy(&reg[n],&alt[e],l);
      memcpy(&reg[s],&alt[o],l);
   }
   memcpy(alt,reg,Hors * sizeof(double));
}


void PC_Var::Reg2Alt(double *reg)
{
   int l,n,s,e,o;
   double alt[Hors];

   l = Lons * sizeof(double);

   for (e=0,o=Lons,n=0,s=Hors-Lons ; e<Hors ; e+=2*Lons,o+=2*Lons,n+=Lons,s-=Lons)
   {
      memcpy(alt+e,reg+n,l);
      memcpy(alt+o,reg+s,l);
   }
   memcpy(reg,alt,sizeof(alt));
}


// Transform a latitude array from alternating to regular order

void PC_Var::Alt2CS(double *alt)
{
   int n,s,e,o;
   double reg[Lats];

   for (e=0,o=1,n=0,s=Lats-1 ; e<Lats ; e+=2,o+=2,n++,s--)
   {
      reg[n] = alt[e];
      reg[s] = alt[o];
   }
   memcpy(alt,reg,Lats * sizeof(double));
}


void PC_Var::ReadGP(void)
{
   char filename[80];
   char tb[80];
   FILE *fp;
   int h[8];
   valarray<double>f(Hors);

   if (mypid == NROOT)
   {
      sprintf(filename,"N%3.3d_surf_%4.4d.sra",Lats,Code);
      fp = fopen(filename,"r");
      if (fp)
      {
         for (int i=0 ; i<   8 ; ++i) fscanf(fp,"%d",h+i);
         for (int i=0 ; i<Hors ; ++i) fscanf(fp,"%lf",&f[i]);
         fclose(fp);
         Reg2Alt(&f[0]);
         Stars(28);
         NewLine();
         sprintf(tb,"Read: %s",filename);
         LeftText(tb,28);
         Dashes(28);
         NewLine();
         sprintf(tb,"Min : %10.2lf",f.min());
         LeftText(tb,28);
         sprintf(tb,"Max : %10.2lf",f.max());
         LeftText(tb,28);
         Stars(28);
      }
      else
      {
         f = 0.0;
      }
   }
   if (Procs == 1) gp = f;
   else           mpscdn(&f[0],&gp[0],Hopp);
}


void PC_Var::WriteFC(const char *prefix)
{
   char filename[80];
   char tb[80];
   FILE *fp;
   valarray<double>f(Hors);

   if (Procs == 1) f = fc;
   else            mpgadn(&f[0],&fc[0],Hopp);
   if (mypid == NROOT)
   {
      Alt2Reg(&f[0]);
      sprintf(filename,"%sF%3.3d_surf_%4.4d.sra",prefix,Lats,Code);
      fp = fopen(filename,"w");
      for (int i=0 ; i<   8 ; ++i) fprintf(fp,"%10ld",Head[i]);
      fputc('\n',fp);
      for (int i=0 ; i<Hors ; ++i)
      {
         fprintf(fp,"%10.2lf",f[i]);
         if (i % 8 == 7) fputc('\n',fp);
      }
      fclose(fp);
   }
}


void PC_Var::WriteGP(const char *prefix)
{
   char filename[80];
   char tb[80];
   FILE *fp;
   valarray<double>f(Hors);

   if (Procs == 1) f = gp;
   else           mpgadn(&f[0],&gp[0],Hopp);
   if (mypid == NROOT)
   {
      Alt2Reg(&f[0]);
      sprintf(filename,"%sF%3.3d_surf_%4.4d.sra",prefix,Lats,Code);
      fp = fopen(filename,"w");
      for (int i=0 ; i<   8 ; ++i) fprintf(fp,"%10ld",Head[i]);
      fputc('\n',fp);
      for (int i=0 ; i<Hors ; ++i)
      {
         fprintf(fp,"%10.2lf",f[i]);
         if (i % 8 == 7) fputc('\n',fp);
      }
      fclose(fp);
   }
}


void PC_Var::WriteServiceGP(const char *prefix)
{
   char filename[80];
   char tb[80];
   FILE *fp;
   int hfcw;
   int afcw;

   valarray<double>f(Hors);

   if (Procs == 1) f = gp;
   else           mpgadn(&f[0],&gp[0],Hopp);
   if (mypid == NROOT)
   {
      hfcw = sizeof(Head);
      afcw = Hors * sizeof(double);
      Alt2Reg(&f[0]);
      sprintf(filename,"%sF%3.3d_surf_%4.4d.srv",prefix,Lats,Code);
      fp = fopen(filename,"w");
      fwrite(&hfcw,4,1,fp);
      fwrite(&Head[0],sizeof(long),8,fp);
      fwrite(&hfcw,4,1,fp);
      fwrite(&afcw,4,1,fp);
      fwrite(&f[0],sizeof(double),Hors,fp);
      fwrite(&afcw,4,1,fp);
      fclose(fp);
   }
}


void PC_Var::WriteSP(void)
{
   int i;
   char filename[80];
   FILE *fp;
   long h[8];
   valarray<double>f(nesp);

   if (mypid == NROOT)
   {
      memcpy(h,Head,sizeof(h));
      h[4] = Spec;
      h[5] =    1;

      sprintf(filename,"T%3.3d_surf_%4.4d.sra",Trit,Code);
      fp = fopen(filename,"w");
      for (i=0 ; i<   8 ; ++i) fprintf(fp,"%10ld",h[i]);
      fputc('\n',fp);
      for (i=0 ; i<Spec ; ++i)
      {
         fprintf(fp,"%10.2lf",sh[i]);
         if (i % 8 == 7) fputc('\n',fp);
      }
      if (i % 8) fputc('\n',fp);
      fclose(fp);
   }
}


// Extract a zonal mean cross section from an array of fourier coefficients
 
void PC_Var::Section(void)
{
   valarray<double>p(Lapp); // Zonal means for one level and process

   for (int jlev=0 ; jlev < VarLevs ; ++jlev)
   {
      p = fc[slice(jlev*Hopp,Lapp,Lons)]; // Extract 0th. fourier coefficient
      mpgadn(&cs[jlev*Lats],&p[0],Lapp);  // Gather
      if (mypid == NROOT) Alt2CS(&cs[jlev*Lats]);
   }
}


// -------------------------- FFT -------------------------------------


// =======
// FFT_ini
// =======

void PC_Var::FFT_ini(void)
{
   FFT_cos = new double[Lons/2];
   FFT_sin = new double[Lons/2];

   for (int k=0 ; k < Lons/2 ; ++k)
   {
      FFT_cos[k] = cos((k * 2.0 * M_PI) / Lons);
      FFT_sin[k] = sin((k * 2.0 * M_PI) / Lons);
   }
}


// ======
// FFT_d8
// ======

void PC_Var::FFT_d8(double *a)
{
   int i0,i1,i2,i3,i4,i5,i6,i7,la;
   double a0p4,a1p5,a2p6,a3p7,a5m1,a7m3,a0m4,a6m2;
   double a0p4p2p6,a1p5p3p7,a7m3p5m1,a7m3m5m1;
   double z,zsin45;

   la = Lons / 8;
   z  = 1.0 / Lons; 
   zsin45 = z * M_SQRT1_2;

   for (i0 = 0 ; i0 < la ; ++i0)
   {
      i1 = i0 + la;
      i2 = i1 + la;
      i3 = i2 + la;
      i4 = i3 + la;
      i5 = i4 + la;
      i6 = i5 + la;
      i7 = i6 + la;
   
      a0p4 =  a[i0] + a[i4];
      a1p5 =  a[i1] + a[i5];
      a2p6 =  a[i2] + a[i6];
      a3p7 =  a[i3] + a[i7];
      a5m1 =  a[i5] - a[i1];
      a7m3 =  a[i7] - a[i3];
      a0m4 = (a[i0] - a[i4]) * z;
      a6m2 = (a[i6] - a[i2]) * z;
   
      a0p4p2p6 = a0p4 + a2p6;
      a1p5p3p7 = a1p5 + a3p7;
      a7m3p5m1 = (a7m3 + a5m1) * zsin45;
      a7m3m5m1 = (a7m3 - a5m1) * zsin45;
   
      a[i0] = z * (a0p4p2p6 + a1p5p3p7);
      a[i7] = z * (a0p4p2p6 - a1p5p3p7);
      a[i3] = z * (a0p4 - a2p6);
      a[i4] = z * (a3p7 - a1p5);
      a[i1] = a0m4 + a7m3m5m1;
      a[i5] = a0m4 - a7m3m5m1;
      a[i2] = a7m3p5m1 + a6m2;
      a[i6] = a7m3p5m1 - a6m2;
   }
}


// ======
// FFT_e2
// ======

void PC_Var::FFT_e2(double *a)
{
   int i,ja,jb;
   double co,si,a1p3,a3m1;
   double c[Lons];

   c[0] = a[0] + a[1];
   c[1] = 0.0 ;

   for (i=2,ja=2,jb=Lons-2 ; i <= Lons-6 ; i+=4,ja+=2,jb-=2)
   {
      co = FFT_cos[ja/2];
      si = FFT_sin[ja/2];
      a1p3 = co * a[i+1] + si * a[i+3];
      a3m1 = co * a[i+3] - si * a[i+1];
      c[ja  ] = a[i] + a1p3;
      c[jb  ] = a[i] - a1p3;
      c[ja+1] = a3m1 + a[i+2];
      c[jb+1] = a3m1 - a[i+2];
   }

   c[ja  ] =  a[Lons-2];
   c[ja+1] = -a[Lons-1];

   memcpy(a,c,sizeof(c));
}


// ======
// FFT_d2
// ======

void PC_Var::FFT_d2(double *a, int la)
{
   int k,l,ja,jb,ia,i0,i1,i2,i3;
   double co,si,a1p3,a3m1;
   double c[Lons];

   for (ja=0,jb=Lons-la ; ja<la ; ++ja,++jb)
   {
      c[ja] = a[ja] + a[ja+la];
      c[jb] = a[ja] - a[ja+la];
   }

   ja = la;
   jb = Lons - 3 * la;
   for (k=la,ia=2*la ; k <= (Lons-2)/4 ; k+=la,ia+=4*la)
   {
      co = FFT_cos[k];
      si = FFT_sin[k];

      for (i0=ia ; i0<ia+la ; ++i0,++ja,++jb)
      {
         i1 = i0 + la;
         i2 = i1 + la;
         i3 = i2 + la;

         a1p3 = co * a[i1] + si * a[i3];
         a3m1 = co * a[i3] - si * a[i1];

         c[ja   ] = a[i0] + a1p3;
         c[jb   ] = a[i0] - a1p3;
         c[ja+la] = a3m1 + a[i2];
         c[jb+la] = a3m1 - a[i2];
      }
      ja += la;
      jb -= 3 * la;
   }

   for (i0=ia ; i0<ia+la ; ++i0,++ja)
   {
      c[ja   ] =  a[i0];
      c[ja+la] = -a[i0+la];
   }
   memcpy(a,c,sizeof(c));
}


// =====
// gp2fc
// =====

void PC_Var::gp2fc(void)
{
   fc = gp;
   for (int j=0 ; j < DimG ; j += Lons)
   {
      FFT_d8(&fc[j]);
      for (int k=Lons>>4 ; k>1 ; k>>=1) FFT_d2(&fc[j],k);
      FFT_e2(&fc[j]);
   }
}


// ======
// FFT_i8
// ======

void PC_Var::FFT_i8(double *a)
{
   int i0,i1,i2,i3,i4,i5,i6,i7,la;
   double a0p7,a0m7,a1p5,a1m5,a2p6,a2m6;
   double a0p7p3,a0p7m3,a0m7p4,a0m7m4;
   double a1m5p2p6,a1m5m2p6;

   la = Lons / 8;

   for (i0=0 ; i0<la ; ++i0)
   {
      i1 = i0 + la;
      i2 = i1 + la;
      i3 = i2 + la;
      i4 = i3 + la;
      i5 = i4 + la;
      i6 = i5 + la;
      i7 = i6 + la;
   
      a0p7 = a[i0] + a[i7];
      a0m7 = a[i0] - a[i7];
      a1p5 = a[i1] + a[i5];
      a1m5 = a[i1] - a[i5];
      a2p6 = a[i2] + a[i6];
      a2m6 = a[i2] - a[i6];

      a0p7p3   = a0p7 + a[i3];
      a0p7m3   = a0p7 - a[i3];
      a0m7p4   = 2.0 * (a0m7 + a[i4]);
      a0m7m4   = 2.0 * (a0m7 - a[i4]);
      a1m5p2p6 = M_SQRT2 * (a1m5 + a2p6);
      a1m5m2p6 = M_SQRT2 * (a1m5 - a2p6);

      a[i0]  = 2.0 * (a0p7p3 + a1p5);
      a[i2]  = 2.0 * (a0p7m3 - a2m6);
      a[i4]  = 2.0 * (a0p7p3 - a1p5);
      a[i6]  = 2.0 * (a0p7m3 + a2m6);

      a[i1]  = a0m7m4 + a1m5m2p6;
      a[i3]  = a0m7p4 - a1m5p2p6;
      a[i5]  = a0m7m4 - a1m5m2p6;
      a[i7]  = a0m7p4 + a1m5p2p6;
   }
}


// ======
// FFT_s2
// ======

void PC_Var::FFT_s2(double *a)
{
   int j,ia,ib;
   double co,si,amb,apb;
   double c[Lons];

   c[1] = c[0] = 0.5 * a[0];

   for (j=2,ia=2,ib=Lons-2 ; j <= Lons-6 ; j+=4,ia+=2,ib-=2)
   {
      co = FFT_cos[ia/2];
      si = FFT_sin[ia/2];
      amb = a[ia  ] - a[ib  ];
      apb = a[ia+1] + a[ib+1];
      c[j  ] = a[ia  ] + a[ib  ];
      c[j+2] = a[ia+1] - a[ib+1];
      c[j+1] = co * amb - si * apb;
      c[j+3] = si * amb + co * apb;
   }
   c[Lons-2] =  a[ia  ];
   c[Lons-1] = -a[ia+1];

   memcpy(a,c,sizeof(c));
}


// ======
// FFT_i2
// ======

void PC_Var::FFT_i2(double *a, int la)
{
   int j,k,a0,a1,a2,a3,b0,b1,b2,b3;
   double co,si,a0m2,a1p3;
   double b[Lons];

   for (j=0 ; j < la ; ++j)
   {
      b[j   ] = a[j] + a[Lons-la+j];
      b[j+la] = a[j] - a[Lons-la+j];
   }

   for (k=la ; k < Lons/4 ; k+=la)
   {   
      a1 = k  *  2;
      a0 = a1 - la;
      a3 = Lons  - a1;
      a2 = a3 - la;

      b0 = a0 + a0;
      b1 = b0 + la;
      b2 = b1 + la;
      b3 = b2 + la;

      co = FFT_cos[k];
      si = FFT_sin[k];

      for (j=0 ; j < la ; ++j)
      {   
         a0m2    = a[a0+j] - a[a2+j];
         a1p3    = a[a1+j] + a[a3+j];
         b[b0+j] = a[a0+j] + a[a2+j];
         b[b2+j] = a[a1+j] - a[a3+j];
         b[b1+j] = co * a0m2 - si * a1p3;
         b[b3+j] = co * a1p3 + si * a0m2;
      }   
   }   

   a1 = Lons  /  2;
   a0 = a1 - la;
   b1 = Lons  - la;
   b0 = b1 - la;

   for (j=0 ; j < la ; ++j)
   {   
      b[b0+j] =  a[a0+j];
      b[b1+j] = -a[a1+j];
   }   

   memcpy(a,b,sizeof(b));
}


// =====
// fc2gp
// =====


void PC_Var::fc2gp(void)
{
   gp = fc;
   for (int j=0 ; j < DimG ; j += Lons)
   {
      FFT_s2(&gp[j]);
      for (int k=2 ; k < Lons/8 ; k<<=1) FFT_i2(&gp[j],k);
      FFT_i8(&gp[j]);
   }
}


// -------------------------- GAUSS -----------------------------------

#define GAUSS_ITER 50
#define GAUSS_EPS  1.0e-16

double ql(int k, double p)
{
   double z0,z1,z2,z3,z4;
   int j;
   z0 = acos(p);
   z1 = 1.0;
   z2 = 0.0;
   for (j=k ; j >= 0 ; j-=2)
   {
      z3  = z1 * cos(z0 * j);
      z2 += z3;
      z4  = (k-j+1) * (k+j) * 0.5;
      z1 *= z4 / (z4 + (j-1));
   }
   if (k % 2 == 0) z2 -= 0.5 * z3;

   z0 = M_SQRT2;
   for (j=1 ; j <= k ; ++j)
   {
      z0 *= sqrt(1.0 - 0.25 / (j*j));
   }
   return z0 * z2;
}


double qld(int k, double p)
{
   double z;

   z = p * ql(k,p) - sqrt((k + k + 1.0) / (k + k - 1.0)) * ql(k-1,p);
   return (p * p - 1.0) / (k * z);
}


// ======
// inigau
// ======

void inigau(int klat,double *pz0, double *pzw)
{
   //  klat          Number of Gaussian latitudes
   //  pz0[klat]     Gaussian abscissas
   //  pzw[klat]     Gaussian weights
   int jlat ;     // Latitudinal loop index
   int jiter;     // Iteration loop index
   double z0,z1,z2,z3,z4,z5;
   
   // Compute Gaussian abscissas & weights for alternating lats
   
   z0 = M_PI / (2*klat+1);
   z1 = 1.0  / (klat*klat*8);
   z4 = 2.0  / (klat*klat);

   for (jlat=0 ; jlat < klat ; jlat+=2)
   {
      z2 = z0 * (jlat + 1.5);
      z2 = cos(z2 + z1 / tan(z2));
      jiter = 0;
      for (jiter=0 ; jiter < GAUSS_ITER ; ++jiter)
      {
         z3  = ql(klat,z2) * qld(klat,z2);
         z2 -= z3;
         if (fabs(z3) < GAUSS_EPS) break; // converged
      }
      z5 = ql(klat-1,z2) / sqrt(klat - 0.5);
      pz0[jlat] = z2;
      pz0[jlat+1] = -z2;
      pzw[jlat] = pzw[jlat+1] = z4 * (1.0 - z2 * z2) / (z5 * z5);
   }
}


// -------------------------- LEGENDRE --------------------------------


// ************************
// * Legendre Polynomials *
// ************************

double *qi; // P(m,n) = Associated Legendre Polynomials
double *qj; // Q(m,n) = Used for d/d(mu)
double *qc; // P(m,n) * gwd              used in fc2sh
double *qe; // Q(mn,) * gwd / cos2       used in mktend
double *qm; // P(m,n) * gwd / cos2 * m   used in mktend
double *qq; // P(m,n) * gwd / cos2 * n * (n+1) / 2  "
double *qu; // P(m,n) / (n*(n+1)) * m    used in dz2uv
double *qv; // Q(m,n) / (n*(n+1))        used in dz2uv

// ======
// legini
// ======

void legini(void)
{

   int jlat;  // Latitude;
   int lm;    // mode index
   int m;     // zonal wavenumber
   int n;     // total wavenumber
   int j;     // polynom index

   double amsq;
   double z1;
   double z2;
   double z3;
   double f1m;
   double f2m;
   double znn1;
   double zsin;    // sin
   double zcsq;    // cos2
   double zcos;    // cos
   double zgwd;    // gw
   double zgwdcsq; // gw / cos2
 
   double zpli[ncsp];
   double zpld[ncsp];

   FILE *fp;

   if (PrintPoly && mypid == NROOT) fp = fopen("poly21.p+","w");
   qi = new double[ncsp * Happ];
   qj = new double[ncsp * Happ];
   qc = new double[ncsp * Happ];
   qe = new double[ncsp * Happ];
   qm = new double[ncsp * Happ];
   qq = new double[ncsp * Happ];
   qu = new double[ncsp * Happ];
   qv = new double[ncsp * Happ];

   for (jlat=0,j=0 ; jlat < Lapp ; jlat+=2)
   {
      // set p(0,0) and p(0,1);

      zgwd    = gaw[jlat];            // gaussian weight - from inigau
      zsin    = gsi[jlat];            // sin(phi) - from inigau
      zcsq    = 1.0 - zsin * zsin;    // cos(phi) squared
      zgwdcsq = zgwd / zcsq;          // weight / cos squared
      zcos    = sqrt(zcsq);           // cos(phi)
      f1m     = sqrt(1.5);
      zpli[0] = M_SQRT1_2;            // sqrt(0.5)
      zpli[1] = f1m * zsin;
      zpld[0] = 0.0;
      lm      = 1;

      // loop over wavenumbers;

      for (m=0 ; m <= Trit ; ++m)
      {
         if (m > 0)
         {
            lm++;
            f2m = -f1m * sqrt(zcsq / (m+m));
            f1m =  f2m * sqrt(m+m + 3.0);
            zpli[lm] = f2m;
            if (lm < ncsp-1)
            {
               zpld[  lm] =  -m * f2m * zsin;
               zpli[++lm] =       f1m * zsin;
            } // (lm < ncsp-1)
         } // (m > 0)

         amsq = m * m;

         for (n=m+2 ; n <= Trit ; ++n)
         {
            z1 = sqrt(((n-1)*(n-1) - amsq) / (4*(n-1)*(n-1)-1));
            z2 = zsin * zpli[lm] - z1 * zpli[lm-1];
            zpld[  lm] = (1-n) * z2 + n * z1 * zpli[lm-1];
            zpli[++lm] = z2 * sqrt((4*n*n-1) / (n*n-amsq));
         } // n

         if (lm < ncsp) // mode (m,Trit);
         {
            z3 = sqrt((Trit*Trit-amsq) / (4*Trit*Trit-1));
            zpld[lm]=-Trit*zsin*zpli[lm] + (Trit+Trit+1)*zpli[lm-1]*z3;
         }
         else           // mode (Trit,Trit);
            zpld[lm]=-Trit*zsin*zpli[lm];
      } // m

      for (m=0,lm=0 ; m <= Trit ; ++m)
      {
         for (n=m ; n <= Trit ; ++n,++lm,++j)
         {
              znn1 = 0.0;
              if (n > 0) znn1 = 1.0 / (n*(n+1));
              qi[j] = zpli[lm];
              qj[j] = zpld[lm];
              qc[j] = zpli[lm] * zgwd;
              qu[j] = zpli[lm] * znn1 * m;
              qv[j] = zpld[lm] * znn1;
              qe[j] = zpld[lm] * zgwdcsq;
              qq[j] = zpli[lm] * zgwdcsq * n * (n+1) * 0.5;
              qm[j] = zpli[lm] * zgwdcsq * m;
              if (PrintPoly && mypid == NROOT)
                 fprintf(fp,"%3d %3d %3d %3d %12.5lf %12.5lf %12.5lf %12.5lf\n",
                         jlat/2,lm,m,n,qi[j],qj[j],qq[j],qc[j]);
         } // n
      } // m
   } // jlat
   if (PrintPoly && mypid == NROOT) fclose(fp);
}


// =====
// fc2sh
// =====

void fc2sh(double *fc, double *sh)
{
   int l;      // Index for latitude
   int m;      // Index for zonal wavenumber
   int n;      // Index for total wavenumber
   double *p;  // pointer to sh
   double *q;  // pointer to polynomial
   double frs,fis,fra,fia;
   
   memset(sh,0,nesp * sizeof(double));
   
   for (l=0,q=qc ; l < Happ ; ++l)
   {
      for (m=0,p=sh ; m <= Trit ; ++m)
      {
         frs = fc[2*(l*Lons+m)  ] + fc[Lons+2*(l*Lons+m)  ];
         fis = fc[2*(l*Lons+m)+1] + fc[Lons+2*(l*Lons+m)+1];
         fra = fc[2*(l*Lons+m)  ] - fc[Lons+2*(l*Lons+m)  ];
         fia = fc[2*(l*Lons+m)+1] - fc[Lons+2*(l*Lons+m)+1];

         for (n=m ; n < Trit ; n+=2)
         {
            *p++ += *q * frs; // real symm
            *p++ += *q * fis; // imag symm
            ++q;
            *p++ += *q * fra; // real anti
            *p++ += *q * fia; // imag anti
            ++q;
         }
         if (n == Trit)
         {
            *p++ += *q * frs; // real symm
            *p++ += *q * fis; // imag symm
            ++q;
         }
      } // m
   } // l
   mpsush(sh,nesp);
   mpbrdn(sh,nesp);
}


// =====
// sh2fc
// =====

void PC_Var::sh2fc(void)
{
   int v;      // index for level
   int l;      // index for latitude
   int m;      // index for zonal wavenumber
   int n;      // index for total wavenumber
   double *f;  // pointer to fc
   double *p;  // pointer to sh
   double *q;  // pointer to polynomial
   double fnr,fni,fsr,fsi;
   
   memset(&fc[0],0,DimG * sizeof(double));
   
   for (v = 0 ; v < VarLevs ; ++v)
   {
      for (l=0,q=qi ; l < Happ ; ++l)
      {
         p = &sh[0] + v * Spec;
         f = &fc[0] + v * Hopp + 2 * Lons * l;
         for (m=0 ; m <= Trit ; ++m)
         {
            fnr = 0.0;
            fni = 0.0;
            fsr = 0.0;
            fsi = 0.0;
   
            for (n=m ; n < Trit ; n+=2)
            {
               fnr += *q * *p;
               fsr += *q * *p;
               ++p;
               fni += *q * *p;
               fsi += *q * *p;
               ++p;
               ++q;
               fnr += *q * *p;
               fsr -= *q * *p;
               ++p;
               fni += *q * *p;
               fsi -= *q * *p;
               ++p;
               ++q;
            }
            if (n == Trit)
            {
               fnr += *q * *p;
               fsr += *q * *p;
               ++p;
               fni += *q * *p;
               fsi += *q * *p;
               ++p;
               ++q;
            }
            f[     0] = fnr;
            f[     1] = fni;
            f[Lons  ] = fsr;
            f[Lons+1] = fsi;
            f += 2;
         } // m
      } // l
   } // v
}


// ========
// sh2fcdmu
// ========

void PC_Var::sh2fcdmu(double *fcdmu)
{
   int l;      // index for latitude
   int m;      // index for zonal wavenumber
   int n;      // index for total wavenumber
   double *f;  // pointer to fc
   double *p;  // pointer to sh
   double *q;  // pointer to polynomial
   double fnr,fni,fsr,fsi;
   
   memset(&fcdmu[0],0,Hopp * sizeof(double));
   
   for (l=0,q=qj ; l < Happ ; ++l)
   {
      p = &sh[0];
      f = &fcdmu[0] + 2 * Lons * l;
      for (m=0 ; m <= Trit ; ++m)
      {
         fnr = 0.0;
         fni = 0.0;
         fsr = 0.0;
         fsi = 0.0;

         for (n=m ; n < Trit ; n+=2)
         {
            fnr += *q * *p;
            fsr -= *q * *p;
            ++p;
            fni += *q * *p;
            fsi -= *q * *p;
            ++p;
            ++q;
            fnr += *q * *p;
            fsr += *q * *p;
            ++p;
            fni += *q * *p;
            fsi += *q * *p;
            ++p;
            ++q;
         }
         if (n == Trit)
         {
            fnr += *q * *p;
            fsr -= *q * *p;
            ++p;
            fni += *q * *p;
            fsi -= *q * *p;
            ++p;
            ++q;
         }
         f[     0] = fnr;
         f[     1] = fni;
         f[Lons  ] = fsr;
         f[Lons+1] = fsi;
         f += 2;
      } // m
   } // l
}


// ======
// dz2uv
// ======

void dz2uv(double *d, double *z, double *u, double *v)
{
   int e;      // index for level
   int l;      // index for latitude
   int m;      // index for zonal wavenumber
   int n;      // index for total wavenumber
   double *uf; // pointer to fc of u
   double *vf; // pointer to fc of v
   double *dp; // pointer to sh of divergence
   double *zp; // pointer to sh of vorticity
   double *uq; // pointer to u-polynomial
   double *vq; // pointer to v-polynomial
   double unr,uni,usr,usi;
   double vnr,vni,vsr,vsi;
   double udr,udi,vdr,vdi;
   double uzr,uzi,vzr,vzi;
   double zsave;
   
   memset(u,0,Hopp * Levs * sizeof(double));
   memset(v,0,Hopp * Levs * sizeof(double));
   
   for (e = 0 ; e < Levs ; ++e)
   {
      zsave = z[2+e*Spec];   // Save mode(0,1) of vorticity
      z[2+e*Spec] -= plavor; // Convert to relative vorticity

      for (l=0 , uq=qu , vq=qv ; l < Happ ; ++l)
      {
         dp = d + e * Spec;
         zp = z + e * Spec;
         uf = u + e * Hopp + 2 * Lons * l;
         vf = v + e * Hopp + 2 * Lons * l;
         for (m=0 ; m <= Trit ; ++m)
         {
            unr = 0.0;
            uni = 0.0;
            usr = 0.0;
            usi = 0.0;
            vnr = 0.0;
            vni = 0.0;
            vsr = 0.0;
            vsi = 0.0;
   
            for (n=m ; n < Trit ; n+=2)
            {
               // symmetric mode (m+n) even

               udr = *uq * *dp;
               vzr = *vq * *zp;
               uzr = *uq * *zp;
               vdr = *vq * *dp;
               ++dp;
               ++zp;
               udi = *uq * *dp;
               vzi = *vq * *zp;
               uzi = *uq * *zp;
               vdi = *vq * *dp;
               ++dp;
               ++zp;
               ++uq;
               ++vq;
               unr += vzr + udi;
               uni += vzi - udr;
               usr -= vzr + udi;
               usi -= vzi - udr;
               vnr += uzi - vdr;
               vni -= uzr - vdi;
               vsr += uzi + vdr;
               vsi -= uzr + vdi;

               // antisymmetric mode (m+n) odd
               
               udr = *uq * *dp;
               vzr = *vq * *zp;
               uzr = *uq * *zp;
               vdr = *vq * *dp;
               ++dp;
               ++zp;
               udi = *uq * *dp;
               vzi = *vq * *zp;
               uzi = *uq * *zp;
               vdi = *vq * *dp;
               ++dp;
               ++zp;
               ++uq;
               ++vq;
               unr += vzr + udi;
               uni += vzi - udr;
               usr += vzr - udi;
               usi += vzi + udr;
               vnr += uzi - vdr;
               vni -= uzr - vdi;
               vsr -= uzi - vdr;
               vsi += uzr - vdi;
            }
            if (n == Trit)
            {
               udr = *uq * *dp;
               vzr = *vq * *zp;
               uzr = *uq * *zp;
               vdr = *vq * *dp;
               ++dp;
               ++zp;
               udi = *uq * *dp;
               vzi = *vq * *zp;
               uzi = *uq * *zp;
               vdi = *vq * *dp;
               ++dp;
               ++zp;
               ++uq;
               ++vq;
               unr += vzr + udi;
               uni += vzi - udr;
               usr -= vzr + udi;
               usi -= vzi - udr;
               vnr += uzi - vdr;
               vni -= uzr - vdi;
               vsr += uzi + vdr;
               vsi -= uzr + vdi;
            }
            uf[     0] = unr;
            uf[     1] = uni;
            uf[Lons  ] = usr;
            uf[Lons+1] = usi;
            vf[     0] = vnr;
            vf[     1] = vni;
            vf[Lons  ] = vsr;
            vf[Lons+1] = vsi;
            uf += 2;
            vf += 2;
         } // m
      } // l
      z[2+e*Spec] = zsave; // Restore absolute vorticity
   } // e
}


// ------------------------- main -------------------------------------

PC_Var C; // Constant part of Tr
PC_Var D; // Divergence
PC_Var M; // Modulated part of Tr
PC_Var O; // Orography
PC_Var P; // Ln Ps
PC_Var T; // Temperature
PC_Var Z; // Vorticity
PC_Var U; // Zonal wind
PC_Var V; // Meridional wind

void Usage(void)
{
   printf("\n\nUsage: puma.x <Lats> <Levs>\n\n");
   exit(1);
}


void process_arguments(int argc, char *argv[])
{
   if (argc < 3) Usage();

   Lats = atoi(argv[1]);
   Levs = atoi(argv[2]);
   Lons = Lats * 2;
   Hors = Lons * Lats;
   Hopp = Lons * Lats;
   Trit = (Lons - 1) / 3;
   Spec = (Trit + 1) * (Trit + 2);
   ncsp = Spec / 2;
   nesp = Spec;
   nspp = nesp;
   Lapp = Lats;
   Happ = Lapp / 2;

   if (Lats != 32) Abort("Illegal value for <Lats>");
   if (Levs <   5) Abort("Illegal value for <Levs>");
}

   
void resource_stats(void)
{
   char tb[80];
   double ut,st;
   struct rusage ru;
   getrusage(RUSAGE_SELF,&ru);
   ut = ru.ru_utime.tv_sec + 0.000001 * ru.ru_utime.tv_usec;
   st = ru.ru_stime.tv_sec + 0.000001 * ru.ru_stime.tv_usec;

   Stars(37);
   NewLine();
   sprintf(tb,"User     time: %10.3lf seconds",ut);
   LeftText(tb,37);
   sprintf(tb,"System   time: %10.3lf seconds",st);
   LeftText(tb,37);
   sprintf(tb,"Total    time: %10.3lf seconds",ut+st);
   LeftText(tb,37);
   if (ru.ru_maxrss > 0) 
   {    
      sprintf(tb,"Memory  usage: %10.3lf MBytes",0.000001 * ru.ru_maxrss);
      LeftText(tb,37);
   }    
   if (ru.ru_minflt > 0) 
   {    
      sprintf(tb,"Page reclaims: %10ld page",ru.ru_minflt);
      if (ru.ru_minflt != 1) strcat(tb,"s");
      LeftText(tb,37);
   }    
   if (ru.ru_majflt > 0) 
   {    
      sprintf(tb,"Page   faults: %10ld page",ru.ru_majflt);
      if (ru.ru_majflt != 1) strcat(tb,"s");
      LeftText(tb,37);
   }    
   if (ru.ru_nswap > 0) 
   {    
      sprintf(tb,"Swaps        : %10ld",ru.ru_nswap);
      LeftText(tb,37);
   }    
   if (ru.ru_inblock > 0) 
   {    
      sprintf(tb,"Disk     read: %10ld block",ru.ru_inblock);
      if (ru.ru_inblock != 1) strcat(tb,"s");
      LeftText(tb,37);
   }    
   if (ru.ru_oublock > 0) 
   {    
      sprintf(tb,"Disk    write: %10ld block",ru.ru_oublock);
      if (ru.ru_oublock != 1) strcat(tb,"s");
      LeftText(tb,37);
   }    
   Stars(37);
   NewLine();
}


void PrintGauss(void)
{
   printf("********************************************\n");
   printf("*  Lat |    Deg |    Sin |    Cos | Weight *\n");
   printf("********************************************\n");
   for (int lat=0 ; lat < Lats ; ++lat)
   {
      printf("* %4d | %6.2lf | %6.3lf | %6.3lf | %6.3lf *\n",
             lat,180.0*asin(gsi[lat])/M_PI,gsi[lat],
             sqrt(1.0-gsi[lat]*gsi[lat]),gaw[lat]);
   }
   printf("********************************************\n\n");
}


void PrintSigmaArrays(void)
{
   printf("********************************************\n");
   printf("*  Lev |   Half |   Full |  Delta |   1/2d *\n");
   printf("********************************************\n");
   for (int lev=0 ; lev < Levs ; ++lev)
   {
      printf("* %4d | %6.3lf | %6.3lf | %6.3lf | %6.3lf *\n",
             lev,SigHalf[lev],SigFull[lev],SigDelt[lev],SigRevd[lev]);
   }
   printf("********************************************\n\n");
}


void Init_Vertical(void)
{
   int lev;

   SigHalf = new double[Levs];
   SigFull = new double[Levs];
   SigDelt = new double[Levs];
   SigRevd = new double[Levs];

   // equidistant half sigma levels

   SigHalf[Levs-1] = 1.0;
   for (lev = 1 ; lev < Levs ; ++lev)
      SigHalf[lev-1] = double(lev) / double(Levs);

   // delta sigma and (1.0 / two delta sigma)

   SigDelt[0] = SigHalf[0];
   SigRevd[0] = 0.5 / SigDelt[0];
   for (lev = 1 ; lev < Levs ; ++lev)
   {
      SigDelt[lev] = SigHalf[lev] - SigHalf[lev-1];
      SigRevd[lev] = 0.5 / SigDelt[lev];
   }

   // full sigma levels

   SigFull[0] = 0.5 * SigHalf[0];
   for (lev = 1 ; lev < Levs ; ++lev)
      SigFull[lev] = (SigHalf[lev] + SigHalf[lev-1]) * 0.5;

   if (PrintSigs) PrintSigmaArrays();
}
   

// *************************************************************
// * Set up the restoration temperature fields C and M         *
// * for aqua planet conditions.                               *
// * The temperature at sigma = 1 is <tgr>, entered in kelvin. *
// * The lapse rate is ALR (standard = 6.5 [K/m] below the     *
// * tropopause and zero above.                                *
// * The tropopause height is defined by <dtrop>.              *
// * The smoothing ot the tropopause depends on <dttrp>.       *
// ************************************************************* 

void Init_Spectral(void)
{
   int    lev;           // level index
   double r_sqrt6;       // factor for mode (0,1)
   double r_sqrtx;       // factor for mode (0,2)
   double HeiLev;        // Height      of current  level
   double HeiPrv;        // Height      of previous level
   double SigPrv;        // Sigma       of previous level
   double TemLay;        // Temperature of current  layer
   double TemPrv;        // Temperature of previous level
   double TemTrp;        // Temperature of tropopause
   double LogSig;        // Logarithm of (sigma previous / sigma current)
   double RdivG;         // Gas constant divided by gravity
   double SigTrp;        // Sigma value of tropopause
   double TemGra;        // Horizontal    temperature gradient
   double TemPro[Levs];  // Mean vertical temperature  profile

   r_sqrt6 =  1.0 / sqrt(6.0);
   r_sqrtx = -2.0 / 3.0 * sqrt(0.4);

   // Initialize with bottom values (lowest half level)

   SigPrv = 1.0; //  sigma value
   TemPrv = tgr; //  Temperature [K]
   HeiPrv = 0.0; //  Height      [m]
   RdivG  = GASCON / GA; // Gas constant divided by gravity

   for (lev = Levs-1 ; lev >= 0 ; --lev) // from bottom to top of atmosphere
   {
      LogSig = log(SigPrv / SigFull[lev]);

      // Perform two iterative steps for smoothing

      HeiLev = HeiPrv + (TemPrv * RdivG) * LogSig;
      TemLay = 0.5 * ALR * (HeiLev - dtrop);
      TemTrp = tgr - dtrop * ALR + sqrt(TemLay*TemLay + dttrp*dttrp) - TemLay;

      HeiLev = HeiPrv + (0.5 * (TemPrv + TemTrp) * RdivG) * LogSig;
      TemLay = 0.5 * ALR * (HeiLev - dtrop);
      TemTrp = tgr - dtrop * ALR + sqrt(TemLay*TemLay + dttrp*dttrp) - TemLay;

      HeiPrv += (0.5 * (TemTrp + TemPrv) * RdivG) * LogSig;
      TemPro[lev] = TemPrv = TemTrp;
      SigPrv = SigFull[lev];
   }

// *********************************************************************
// * loop to set array TemGra - this controls temperature gradients as *
// * a function of sigma in tres. it is a sine wave from one at        *
// * sigma = 1 to zero at SigTrp (sigma at the tropopause).            *
// *********************************************************************

   SigTrp = pow((tgr-dtrop*ALR) / tgr , GA / (ALR*GASCON));

// now the latitudinal variation in tres is set up ( this being in terms
// of a deviation from t0 which is usually constant with height)

   for (lev = 0 ; lev < Levs ; ++lev)
   {
      TemGra = sin(0.5 * M_PI * (SigFull[lev] - SigTrp) / (1.0 - SigTrp));
      if (TemGra < 0.0) TemGra = 0.0;
      C.sh[    lev * nesp] = M_SQRT2 * (TemPro[lev]/CT - t0);  // mode (0,0) mean
      M.sh[2 + lev * nesp] = r_sqrt6 * dtns * TemGra;          // mode (0,1) N-S
      C.sh[4 + lev * nesp] = r_sqrtx * dtep * TemGra;          // mode (0,2) E-P
printf("Temp[%d] = %7.3lf\n",lev,TemPro[lev]);
   }

   O.sh /= CV * CV;    // descale orography
   C.sh /= CT;         // descale constant  part of Tr
   M.sh /= CT;         // descale modulated part of Tr
}


void prolog(void)
{
   StepDelta = 2.0 * M_PI / ntspd; // length of time step

   gsi = new double[Lats];
   gaw = new double[Lats];
   if (mypid == NROOT)
   {
      inigau(Lats,gsi,gaw);
      if (PrintLats) PrintGauss();
   }
   mpscdn(gsi,gsi,Lapp);
   mpscdn(gaw,gaw,Lapp);
   legini();

   dpsdmu.resize(Hopp,0.0);

   O.Init(129,   1,"O","Orography"   );
   P.Init(152,   1,"P","Pressure"    );
   T.Init(130,Levs,"T","Temperature" );
   D.Init(155,Levs,"D","Divergence"  );
   Z.Init(138,Levs,"Z","Vorticity"   );
   U.Init(131,Levs,"U","Zonal wind"  );
   V.Init(132,Levs,"V","Merid. wind" );
   C.Init(153,Levs,"C","Tr Constant" );
   M.Init(154,Levs,"M","Tr Modulated");

   O.ReadGP(); // Read orograpy gridpoint array

   P.Status(stdout,"after Init");
   T.Status(stdout,"after Init");
   D.Status(stdout,"after Init");
   Z.Status(stdout,"after Init");

   Init_Vertical();
   if (mypid == NROOT)
   {
      Init_Spectral();
   }
   C.BroadcastSpectral();
   M.BroadcastSpectral();
   C.StatusSP(stdout,"Tr Const");
   M.StatusSP(stdout,"Tr Modul");
}


void epilog(void)
{
   if (mypid == NROOT)
   {
      resource_stats();
   }
}


// =========
// gridpoint
// =========

void gridpoint(void)
{
   // Transfrom from spherical harmonics to fourier coefficients

   D.sh2fc();
   T.sh2fc();
   P.sh2fc();
   Z.sh2fc();

   // Compute d(LnPs) / d(mu) 

   P.sh2fcdmu(&dpsdmu[0]);

   // Compute fourier coefficients of u * cos(phi) and v * cos(phi)
   // from spherical harmonics of divergence and vorticity

   dz2uv(&D.sh[0],&Z.sh[0],&U.fc[0],&V.fc[0]);

   // Compute zonal mean cross sections for diagnostics and GUI

   U.Section();
   V.Section();
   T.Section();
   
}


// ======
// master
// ======

void master(void)
{
   while (nits > 0) // number of initial time steps (nits)
   {
      StepDelta = (2.0 * M_PI / ntspd) / (1 << nits--) ;

      gridpoint(); // compute nonlinear tendencies

      //   call makebm
      //   call spectral
   }

   StepDelta = 2.0 * M_PI / ntspd;

   O.gp2fc();
   O.StatusFC(stdout,"after gp2fc");
   O.WriteFC("A");

   fc2sh(&O.fc[0],&O.sh[0]);
   O.StatusSP(stdout,"after fc2sh");
   O.WriteSP();

   O.sh2fc();
   O.StatusFC(stdout,"after sh2fc");
   O.WriteFC("B");

   O.fc2gp();
   O.StatusGP(stdout,"after fc2gp");
   O.WriteGP("C");
}


int main(int argc, char *argv[])
{
   process_arguments(argc,argv);

   mp_start(&argc,&argv);
   prolog();
   master();
   epilog();
   mp_stop();

   return 0;
}

