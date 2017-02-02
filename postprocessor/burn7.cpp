/*  c++ -O2 -o burn7.x burn7.cpp -lm -lnetcdf_c++ -lnetcdf */

#define NETCDF_OUTPUT

#define V0 "burn 7.7 (02-Feb-2017)"
#define V1 "KlimaCampus"
#define V2 "Usage: burn7 [-help|-c|-d|-m|-n|-r|-s] <modelfile> <resultfile>"
#define V3 "New: option <-g> writes Grads ctl for service plotting"
#define V4 "New: comments (#) are allowed in namelist file"

/* ============= */
/* include files */
/* ============= */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/resource.h>
#include <vector>
#include <valarray>

#ifdef OPEN_MP
#include <omp.h>
#endif

#ifdef NETCDF_OUTPUT
#include <netcdfcpp.h>
#endif

using namespace std;

#define LONG long long

#ifndef M_PI
#define M_PI 3.1415926535
#endif

#ifndef M_SQRT2
#define M_SQRT2 1.4142136
#endif

#ifndef MAX
#define MAX(v1,v2) ((v1) > (v2) ? (v1) : (v2))
#endif

#ifndef min
#define min(v1,v2) ((v1) < (v2) ? (v1) : (v2))
#endif

#ifndef abs
#define abs(x) ((x) >= 0 ? (x) : -(x))
#endif

#ifndef TRUE
#define TRUE  1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#define LEV_SURFACE      1
#define LEV_99          99
#define LEV_ISOBARIC   100
#define LEV_MEANSEA    102
#define LEV_ALTITUDE   103
#define LEV_HEIGHT     105
#define LEV_SIGMA      107
#define LEV_HYBRID     109
#define LEV_GROUND     112

#define REP_REGULAR      0
#define REP_GAUSS        4
#define REP_SPECTRAL    50

#define MAX_HOURS        4
#define MAX_LEVELS      99

#define MAX_NVCT      (MAX_LEVELS * 2 + 2)

#define L_TIMES_RHOH2O (-333700000.0)
#define EARTH_RADIUS 6371220.0
#define MARS_RADIUS  3400000.0
#define SQUARE_RADIUS (-PlanetRadius * PlanetRadius)
#define EARTH_GRAV 9.80665
#define MARS_GRAV  3.728
#define RG   (1.0 / Grav)
#define MARS_RD     (189.0 )

/* ************************************** */
/* Thermodynamical constants adopted from */
/* ECMWF IFS-Code                         */
/* ************************************** */

#define RKBOL  (1.380658e-23)
#define RNAVO  (6.0221367e+23)
#define R      (RKBOL * RNAVO)
#define RMD    (28.9644)
#define RMV    (18.0153)
#define EARTH_RD (1000. * R / RMD)
#define RV     (1000. * R / RMV)
#define RCPD   (3.5 * RD)
#define RCPV   (4.0 * RV)
#define RETV   (RV / RD - 1.)
#define RCW    (4218.)
#define RCS    (2106.)
#define RTT    (273.16)
#define RLVTT  (2.5008e+6)
#define RLSTT  (2.8345e+6)
#define RESTT  (611.14)
#define RLAPSE (0.0065)

#ifdef NETCDF_OUTPUT

NcFile *NetFile;

NcVar *LonVar;
NcVar *LatVar;
NcVar *LevVar;
NcVar *TimVar;

NcDim *LonDim;
NcDim *LatDim;
NcDim *LevDim;
NcDim *TimDim;

#endif

int     SaveMemory =  0; /* Switch on for dynamic memory usage */
int     PolyCreate =  0; /* Create polynomials files for hires T1365 and more */
int     PolyDisk   =  0; /* Read polynomials from disk */
int     GaussGrid  = -1; /* 1: use Gaussian grid, 0: use regular grid */
int     DPM        =  0; /* Days Per Month */
int     DPY        =  0; /* Days Per Year */
int     DayDivisor =  0; /* Use for day adjustment if more than 99 days per month */

char    VerType; /* s=Sigma   p=Pressure */
char    HorType; /* s=Spherical   f=Fourier   z=Zonal Mean   g=Gauss Grid */
char   *filename;
char    ifile[256]; // name of input file (PUMA-II format)
char    ofile[256]; // Name of output file (Service or NetCDF format)
char    gfile[256]; // Name of Grads control file
char    rfile[256]; // Name of Grads data file (for zonal means only)
#define MAX_NL 40960
char    namelist[MAX_NL];
char    AllocName[256];

int    PumaCode;
int    PumaLevel;
int    RepGrib;
int    Debug;
int    Dim3FC     ;
int    Dim3SP     ;
int    Dim3GG     ;
int    Dim3GP     ;
int    DimFC      ; // Dimension of fourier array
int    DimGP      ; // Dimension of output grid
int    DimGG      ; // Dimension of Gauss  grid
int    DimAB      ; // Dimension of array buffer
int    DimSP      ;
int    DimSP_half ;
int    CoreBigEndian ; // Do we run on a big endian machine ?
int    FileBigEndian ;
int    Endian  = 0   ; /* Marker for reverse endian format */
int    LongFCW = 0   ; /* Flag for 64bit (1) or 32bit (0) FCW's */
int    RealSize = 4  ; /* Size of real data (4 = float) (8 = double) */
int    HeadSize = 32 ; // 32:Service single   64:Service double
int    EndOfMonth    ;
int    EndOfTerm     ;
int    Fouriers      ;
int    HumInfo       ; // Flag for humidity info issued
double   Grav          = EARTH_GRAV;
double   SigmaTop      = 0.0;
int    NetCDF        ;
int    GaussianOutput = 1;
int    Grads         ;
int    HeaderWords   ; /* Length of file header in 32-bit words */
int    Gats     ;
int    Lats     ;
int    AllLevs     = 0   ; // # of sigma levels in data file
int    SigLevs     = 0   ; // # of sigma levels used
int    SingleLevel = 0   ; // Set to true for SAM/SOM models
double LevelFactor = 1.0 ; // Multiplier for head(2)
int    LevelType     ;
int    Gons    ;
int    Lons    ;
int    Cyclical = 0 ;  //   1 = Cyclical completion (Lons from 0 to 360)
double   L_times_rhoH2O = L_TIMES_RHOH2O;
int    mars          ;
int    Mean          ;
int    MeanCount     ; // Count terms during month
int    Multi         ;
int    FirstMonth = 1;
int    LastMonth = 12;
double   PlanetRadius= EARTH_RADIUS;
double   RD          = EARTH_RD;
int    Spectral = FALSE;
int    TermCount     ;
int    MonthCount    ;
int    OutputCount   ;
int    Truncation = 0;
int    Waves         ;
int    SpecialUV     ;

struct tm NewDate;
struct tm OldDate;
int    NewMonth;
int    OldMonth;

// Some functions for a nice printout

#define COLS 72

void Stars(int n)  {while (n--) putchar('*');}
void ErrStars(int n)  {while (n--) putc('*',stderr);}
void Blanks(int n) {while (n--) putchar(' ');}
void Dashes(int n) {while (n--) putchar('-');}
void NewLine(void) {putchar('\n');}
void ErrNewLine(void) {putc('\n',stderr);}
void StarLine(void) {Stars(COLS); NewLine();}

/* ==================================== */
/* Abort - Print error message and exit */
/* ==================================== */

void Abort(const char *errtext)
{
   Stars(min(80,strlen(errtext))); NewLine();
   puts(errtext);                  NewLine();
   Stars(min(80,strlen(errtext))); NewLine();
   ErrStars(min(80,strlen(errtext))); ErrNewLine();
   fputs(errtext,stderr);             ErrNewLine();
   ErrStars(min(80,strlen(errtext))); ErrNewLine();
   exit(1);
}

void BlankLine(void)
{
   putchar('*');
   Blanks(COLS-2);
   putchar('*');
   NewLine();
}

void DashLine(void)
{
   putchar('*');
   putchar(' ');
   Dashes(COLS-4);
   putchar(' ');
   putchar('*');
   NewLine();
}

void CenterText(const char *t)
{
   int i,j,l;
   l = strlen(t);
   if (l < 1) return;
   if (l > COLS-4) puts(t);
   else
   {
      i = (COLS - 4 - l) / 2;
      j = (COLS - 4 - l - i);
      putchar('*');
      Blanks(i+1);
      fputs(t,stdout);
      Blanks(j+1);
      putchar('*');
      NewLine();
   }
}

void LeftText(const char *t)
{
   int l;
   l = strlen(t);
   if (l < 1) return;
   if (l > COLS-4) puts(t);
   else
   {
      putchar('*');
      putchar(' ');
      fputs(t,stdout);
      Blanks(COLS - l - 3);
      putchar('*');
      NewLine();
   }
}


#define MAX_ID_LEN  8
#define MAX_NA_LEN 32
#define MAX_UN_LEN 16

class Control
{
   public:
   int    readit      ;
   int    needed      ;
   int    selected    ;
   int    detected    ;
   int    hlev        ;
   int    plev        ;
   int    loff        ;
   int    twod        ;
   int    code        ;
   valarray<double>hsp;
   valarray<double>hfc;
   valarray<double>hgp;
   valarray<double>pgp;
   valarray<double>mgp;
   valarray<double>pfc;
   valarray<double>psp;
   char Id[MAX_ID_LEN];
   char Na[MAX_NA_LEN];
   char Un[MAX_UN_LEN];
#ifdef NETCDF_OUTPUT
   NcVar *NetVar     ;
#endif

   void Status(void);
   void Init(const char* Idf, const char *Name, const char *Units, int TwoD);
   void SetHSpec(int Hlev, int Plev, int Twod);
   void SetHFour(int Hlev, int Plev, int Twod);
   void SetHGrid(int Hlev, int Plev, int Twod);
   void SetPGrid(void);
   void SetPFour(void);
   void SetPSpec(void);
};

void Control::Status(void)
{
   printf("\nStatus for code %3d: %s\n",code,Id);
   printf("-------------------------\n");
   printf("needed:   %5d\n",needed);
   printf("selected: %5d\n",selected);
   printf("detected: %5d\n",detected);
   printf("hlev:     %5d\n",hlev);
   printf("plev:     %5d\n",plev);
   printf("twod:     %5d\n",twod);
   printf("hsp:      %5ld\n",hsp.size());
   printf("hfc:      %5ld\n",hfc.size());
   printf("hgp:      %5ld\n",hgp.size());
   printf("pgp:      %5ld\n",pgp.size());
   printf("mean:     %5ld\n",mgp.size());
   printf("pfc:      %5ld\n",pfc.size());
   printf("psp:      %5ld\n",psp.size());
   if (hgp.size()) printf("mean of hgp: %16.10lf\n",hgp.sum() / hgp.size());
   if (pgp.size()) printf("mean of pgp: %16.10lf\n",pgp.sum() / pgp.size());
};

void Control::Init(const char* Idf, const char *Name, const char *Units, int TwoD)
{  
   strncpy(Id,Idf  ,MAX_ID_LEN-1);
   strncpy(Na,Name ,MAX_NA_LEN-1);
   strncpy(Un,Units,MAX_UN_LEN-1);
   twod = TwoD;
}

void Control::SetHSpec(int Hlev, int Plev, int Twod)
{
   int OldSize;
   int NewSize;

   hlev  = Hlev;
   plev  = Plev;
   twod  = Twod;

   OldSize = hsp.size();
   NewSize = Hlev * DimSP;

   if (NewSize == OldSize) return;

   hsp.resize(NewSize);

   if (Debug)
   {
      char tb[COLS];
      if (OldSize == 0)
      sprintf(tb,"Alloc:   %p %6s[%3d].hsp   %6ld double",&hsp[0],Id,code,hsp.size());
      else
      sprintf(tb,"Realloc: %p %6s[%3d].hsp   %6ld double",&hsp[0],Id,code,hsp.size());
      LeftText(tb);
   }
}

void Control::SetHFour(int Hlev, int Plev, int Twod)
{
   if (hfc.size() == Hlev * DimFC) return;

   hlev  = Hlev;
   plev  = Plev;
   twod  = Twod;

   hfc.resize(hlev * DimFC);

   if (Debug)
   {
      char tb[COLS];
      sprintf(tb,"Alloc:   %p %6s[%3d].hfc   %6ld double",&hfc[0],Id,code,hfc.size());
      LeftText(tb);
   }
}

void Control::SetPFour(void)
{
   if (pfc.size() == plev * DimFC) return;

   pfc.resize(plev * DimFC);

   if (Debug)
   {
      char tb[COLS];
      sprintf(tb,"Alloc:   %p %6s[%3d].pfc   %6ld double",&pfc[0],Id,code,pfc.size());
      LeftText(tb);
   }
}

void Control::SetHGrid(int Hlev, int Plev, int Twod)
{
   if (hgp.size() == Hlev * DimGP) return;

   hlev  = Hlev;
   plev  = Plev;
   twod  = Twod;
   hgp.resize(hlev * DimGP);

   if (Debug)
   {
      char tb[COLS];
      sprintf(tb,"Alloc:   %p %6s[%3d].hgp   %6d double",&hgp[0],Id,code,hlev*DimGP);
      LeftText(tb);
   }
}
   
void Control::SetPGrid(void)
{
   if (twod && hgp.size())
   {
      pgp.resize(DimGP);
      pgp = hgp;
      hgp.resize(0);
      if (Debug)
      {
         char tb[COLS];
         sprintf(tb,"Moved:   %p %6s[%3d].hgp to pgp",&pgp[0],Id,code);
         LeftText(tb);
      }
   }
   else if (pgp.size() != plev * DimGP)
   {
      pgp.resize(plev * DimGP);
      if (Debug)
      {
         char tb[COLS];
         sprintf(tb,"Alloc:   %p %6s[%3d].pgp   %6d double",&pgp[0],Id,code,plev*DimGP);
         LeftText(tb);
      }
   }
}
   
void Control::SetPSpec(void)
{
   int OldSize;
   int NewSize;

   OldSize = psp.size();
   NewSize = plev * DimSP;

   if (NewSize == OldSize) return;

   psp.resize(NewSize);

   if (Debug)
   {
      char tb[COLS];
      if (OldSize == 0)
      sprintf(tb,"Alloc:   %p %6s[%3d].psp   %6ld double",&psp[0],Id,code,psp.size());
      else
      sprintf(tb,"Realloc: %p %6s[%3d].psp   %6ld double",&psp[0],Id,code,psp.size());
      LeftText(tb);
   }
}

#define CODES 512
#define MAXCODES (CODES+5)

#define GEOSCODE 129
#define    TCODE 130
#define    UCODE 131
#define    VCODE 132
#define   SHCODE 133
#define   PSCODE 134
#define    WCODE 135
#define   WZCODE 137
#define  VORCODE 138
#define  STRCODE 148
#define  VELCODE 149
#define  SLPCODE 151
#define LNPSCODE 152
#define  DIVCODE 155
#define    ZCODE 156
#define   RHCODE 157

int SpecialCodes[4] = {DIVCODE,VORCODE,STRCODE,VELCODE};

Control All[MAXCODES];

Control *Geopotential = &All[129];
Control  *Temperature = &All[130];
Control       *u_wind = &All[131];
Control       *v_wind = &All[132];
Control     *Humidity = &All[133];
Control           *Ps = &All[134];
Control        *Omega = &All[135];
Control       *w_wind = &All[137];
Control    *Vorticity = &All[138];
Control           *Ts = &All[139];
Control      *StreamF = &All[148];
Control      *VeloPot = &All[149];
Control          *SLP = &All[151];
Control         *LnPs = &All[152];
Control   *Divergence = &All[155];
Control *GeopotHeight = &All[156];
Control    *Rhumidity = &All[157];
Control        *speed = &All[259];
Control       *precip = &All[260];
Control      *net_top = &All[261];
Control      *net_bot = &All[262];
Control     *net_heat = &All[263];
Control    *net_water = &All[264];
Control       *sw_atm = &All[268];
Control       *lw_atm = &All[269];
Control      *net_atm = &All[270];
Control  *surf_runoff = &All[271];
Control        *dpsdx = &All[273];
Control        *dpsdy = &All[274];
Control  *fresh_water = &All[275];

Control    *HalfPress = &All[277];
Control    *FullPress = &All[278];
Control       *ThetaH = &All[279];
Control       *ThetaF = &All[280];

int  *vert_index;
valarray<double>Orography;
double *poli;
double *pol2;
double *pliu;
double *pliv;

char polin[80]; /* filenames for polynomial files */
char pol2n[80];
char pliun[80];
char plivn[80];

FILE *polif;
FILE *pol2f;
FILE *pliuf;
FILE *plivf;

int OutLevs;  /* number of requested output   levels  */
int nrpl;     /* number of requested pressure levels  */
int nrml;     /* number of requested model    levels  */
int nrqh;

#define ATTR_MAX 20
int nattr;    /* number of global attributes (NetCDF) */
char AttrNam[ATTR_MAX][80];
char AttrVal[ATTR_MAX][80];

int nvct = 0;
double vct[MAX_NVCT];

int DaysPerYear = 360;
double DataStep = 0.0;
int DeltaDy;
int DeltaHr;
int DeltaMn;

vector<int>HeadSt(8,0); // First  header vector
vector<int>HeadIn(8,0); // Input  header vector
vector<int>HeadOu(8,0); // Output header vector

#define HYB_SPEC  0
#define HYB_FOUR  1
#define HYB_ZONM  2
#define HYB_GRID  3
#define PRE_GRID  4
#define PRE_FOUR  5
#define PRE_ZONM  6
#define PRE_SPEC  7

int OutRep = HYB_SPEC;   // Output Representation 

int hours[MAX_HOURS+1];
double level[MAX_LEVELS+1];
int  LevelFound[MAX_LEVELS+1];

double hPa[MAX_LEVELS];
int  mol[MAX_LEVELS];
int  mom[MAX_LEVELS];        /* Mask for selected levels */
int LevelUnits[MAX_LEVELS];
int SigmaF[MAX_LEVELS];

double *Record_double; 
float  *Record_float;  // Buffer for FORTRAN records
int    *Record_int;    // All share the same space
unsigned short  *Record_short;
char   *Record_char;

double *CosPhi;        /* Cosine of Phi (Latitude)  */
double *RevCosPhi;     /* 1.0 / CosPhi              */
double *DerivationFactor;

double *wfc;     /* FFT work array */
double *wgp;     /* FFT work array */
int   ifax[10]; /* FFT factorization */

FILE  *fpi;
FILE  *gp;
 
class RegLon // class for equidistant longituinal array
{
   char   Name[16];        // Array name
   double DeltaLam;        // Distance of Longitudes

public:
    int     Lons;           // Number of latitudes
    double *Lam;            // Coordinate in degrees

    RegLon(int, const char *);   // Constructor
   ~RegLon(void);          // Destructor
    void PrintArray(void); // Print table
};

RegLon::RegLon(int n, const char *name)
{
   int jlon;

   Lons = n;
   DeltaLam = 360.0 / Lons;
   Lam = new double[Lons+1];
   for (jlon=0 ; jlon < Lons ; ++jlon)
   {
       Lam[jlon] = jlon * DeltaLam;
   }
   Lam[Lons] = 360.0;
   strncpy(Name,name,15);
}

RegLon::~RegLon(void)
{
   delete Lam;
}

void RegLon::PrintArray(void)
{
   int jlon;

   printf("*******************************\n");
   printf("* %16.16s  Longitude *\n",Name);
   printf("*******************************\n");
   for (jlon=0 ; jlon < Lons ; ++jlon)
      printf("* %16d %10.4f *\n",jlon,Lam[jlon]);
   printf("*******************************\n");
}

class RegLat // class for equidistant latitudinal array
{
   char   Name[16];        // Array name
   double DeltaPhi;        // Distance of latitudes
   double FirstPhi;        // First latitude

public:
    int     Lats;          // Number of latitudes
    double *Phi;           // Coordinate in degrees
    double *gmu;           // Sine of phi
    double *gwt;           // Gaussian weight

    RegLat(int, const char *);   // Constructor
   ~RegLat(void);          // Destructor
    void PrintArray(void); // Print table
};

RegLat::RegLat(int n, const char *name)
{
   int jlat;

   Lats = n;
   DeltaPhi = 180.0 / (Lats - 1);
   Phi = new double[Lats];
   gmu = new double[Lats];
   gwt = new double[Lats];

   if (Lats & 1) /* odd -> start with DeltaPhi */
   {
      DeltaPhi = 180.0 / (Lats + 1);
      FirstPhi =  90.0 - DeltaPhi;
   }
   else /* even -> start with 0.5 DeltaPhi */
   {
      DeltaPhi = 180.0 / Lats;
      FirstPhi =  90.0 - 0.5 * DeltaPhi;
   }

   for (jlat=0 ; jlat < Lats ; ++jlat)
   {
       Phi[jlat] = FirstPhi - jlat * DeltaPhi;
       gmu[jlat] = sin(Phi[jlat] * M_PI / 180.0);
       gwt[jlat] = 0.0;
   }
   strncpy(Name,name,15);
}

RegLat::~RegLat(void)
{
   delete Phi;
   delete gmu;
   delete gwt;
}

void RegLat::PrintArray(void)
{
   int jlat;

   printf("*******************************\n");
   printf("* %16.16s   Latitude *\n",Name);
   printf("*******************************\n");
   for (jlat=0 ; jlat < Lats ; ++jlat)
      printf("* %16d %10.4f *\n",jlat,Phi[jlat]);
   printf("*******************************\n");
}

class GauLat : public RegLat
{
public:
    GauLat(int, const char *); // Constructor
   ~GauLat(void);        // Destructor
private:
   void IniGau(void);
   double nlp(int, double);
   double dlp(int, double);
};

GauLat::GauLat(int n, const char *name) : RegLat(n, name)
{
   int jlat;
   IniGau();
   for (jlat = 0 ; jlat < Lats ; ++jlat)
      Phi[jlat] = 180.0 * asin(gmu[jlat]) / M_PI;
}

GauLat::~GauLat(void)
{
}

/* ============================== */
/* Calculate Legendre Polynomials */
/* ============================== */

double GauLat::nlp(int k, double p) // After Nodorp (1988)
{
   int j;
   double z0,z1,z2,z3,z4;

   z0 = acos(p);
   z1 = 1.0;
   z2 = 0.0;
   z3 = 0.0;
   
   for (j=k ; j >= 0 ; j-=2)
   {
      z3  = z1 * cos(j * z0);
      z2 += z3;
      z4  = (k-j+1) * (k+j)/2;
      z1 *= z4 / (z4+j-1);
   }
   if (k % 2 == 0) z2 -= 0.5 * z3;

   z0 = sqrt(2.0);
   for (j=1; j <= k ; ++j)
      z0 *= sqrt(1.0 - 0.25/(j*j));
   return (z0*z2);
}

/* ============================================= */
/* Calculate Derivatives of Legendre Polynomials */
/* ============================================= */

double GauLat::dlp(int k, double p) // After Nodorp (1988)
{
   double z0,z3,z4;

   if (!k) return 0.0;

   z0 = 1.0 / (p*p - 1.0);
   z3 = sqrt((k+k+1.0)/(k+k-1.0));
   z4 = p * nlp(k,p) - z3 * nlp(k-1,p);

   return(k*z0*z4);
}

void GauLat::IniGau(void)
{
   int jlat,Iter;
   double z0, z1, z2, z3;
   double eps=1.e-15;

   z0 = M_PI / (2 * Lats + 1);
   z1 = 1.0  / (8 * Lats * Lats);
   for (jlat=0 ; jlat < Lats/2 ; ++jlat) // North to Equator
   {
      z2 = (2 * jlat + 1.5) * z0;
      z2 = cos(z2 + z1 / tan(z2));
      Iter = 0;
      do
      {
         z3  = nlp(Lats,z2) / dlp(Lats,z2);
         z2 -= z3;
      } while ((z3 < -eps || z3 > eps) && ++Iter < 1000);
      gmu[jlat]     =  z2;
      gmu[Lats-1-jlat] = -z2;
   }
   z1 = 2.0 / (Lats * Lats);
   for (jlat=0 ; jlat < Lats/2 ; ++jlat) // North to Equator
   {
      z0 = nlp(Lats-1,gmu[jlat]) / sqrt(Lats - 0.5);
      z0 = z0 * z0;
      gwt[jlat] = z1 * (1.0 - gmu[jlat] * gmu[jlat]) / z0;
      gwt[Lats-1-jlat] = gwt[jlat];
   }
}

RegLon *Outlon;
RegLat *Outlat;
GauLat *Gaulat;

#ifdef NETCDF_OUTPUT
#define TITLE        "PUMA/PLASIM DATA"
#define HISTORY      "Created by PumaBurner 7.4"


void NetOpen(char *NetFileName)
{
   int yy,mm,dd;
   int jlev;
   int jvar;
   int londim;

   double *Outlev;

   const char *title=TITLE;
   const char *conv="CF-1.0";
   const char *hist=HISTORY;

   char cale[80];
   char t_unit[80];
   int BaseDate[4];

   Outlev = new double[OutLevs];

   if (VerType == 's') // sigma
   {
      for (jlev = 0 ; jlev < OutLevs ; ++jlev)
         Outlev[jlev] = 0.5 *
         (vct[SigLevs+mol[jlev]]+vct[SigLevs+mol[jlev]+1]);
   }
   else // pressure levels [hPa]
   {
      for (jlev = 0 ; jlev < OutLevs ; ++jlev) Outlev[jlev] = hPa[jlev];
   }

   BaseDate[0] = yy = HeadSt[2] / 10000;
   BaseDate[1] = mm = HeadSt[2] / 100 % 100;
   BaseDate[2] = dd = HeadSt[2] % 100;
   BaseDate[3] = 0;
   if (Mean) 
   {
      sprintf(t_unit,"months since %04d-%02d-%02d 00:00:00",yy,mm,dd);
      sprintf(cale,"360_day");
   }
   else
   {
      sprintf(t_unit,"days since %04d-%02d-%02d 00:00:00",yy,mm,dd);
      if (DaysPerYear != 365) sprintf(cale,"%d_day",DaysPerYear);
      else                    sprintf(cale,"proleptic_gregorian");
   }

   /* Create NetCDF file */

   NetFile = new NcFile(NetFileName,NcFile::Replace,NULL,0,NcFile::Offset64Bits);
   if (!NetFile->is_valid()) Abort("Could not open NetCDF file");

   /* Define dimensions */

   if (OutRep == HYB_ZONM || OutRep == PRE_ZONM) londim = 1;
   else                            londim = Lons + Cyclical;

   LonDim = NetFile->add_dim("lon" , londim);
   LatDim = NetFile->add_dim("lat" , Lats  );
   LevDim = NetFile->add_dim("lev" , OutLevs  );
   TimDim = NetFile->add_dim("time"        );

   LonVar = NetFile->add_var("lon" , ncDouble, LonDim);
   LatVar = NetFile->add_var("lat" , ncDouble, LatDim);
   LevVar = NetFile->add_var("lev" , ncDouble, LevDim);
   TimVar = NetFile->add_var("time", ncDouble, TimDim);

   LonVar->add_att("axis"         ,"X"            );
   LonVar->add_att("long_name"    ,"longitude"    );
   LonVar->add_att("standard_name","longitude"    );
   LonVar->add_att("units"        ,"degrees_east" );

   LatVar->add_att("axis"         ,"Y"            );
   LatVar->add_att("long_name"    ,"latitude"     );
   LatVar->add_att("standard_name","latitude"     );
   LatVar->add_att("units"        ,"degrees_north");

   if (VerType == 's') // sigma level
   {
      LevVar->add_att("axis"         ,"Z"                          );
      LevVar->add_att("long_name"    ,"sigma at layer midpoints"   );
      LevVar->add_att("standard_name","atmosphere_sigma_coordinate");
      LevVar->add_att("positive"     ,"down"                       );
      LevVar->add_att("units"        ,"level"                      );
   }
   else // pressure level
   {
      LevVar->add_att("axis"         ,"Z"            );
      LevVar->add_att("long_name"    ,"pressure"     );
      LevVar->add_att("standard_name","pressure"     );
      LevVar->add_att("units"        ,"hPa"          );
   }

   TimVar->add_att("calendar"     ,cale           );
   TimVar->add_att("units"        ,t_unit         );
   TimVar->add_att("base_date", 4 ,BaseDate       );

   
   NetFile->add_att("title"      , title);
   NetFile->add_att("history"    , hist );
   NetFile->add_att("Conventions", conv );

   for (jvar = 0 ; jvar < nattr ; ++jvar)
   {
      NetFile->add_att(AttrNam[jvar],AttrVal[jvar]);
   }
   
   LonVar->put(Outlon->Lam,londim );
   LatVar->put(Outlat->Phi,Lats   );
   LevVar->put(Outlev     ,OutLevs);

}


void NetVarDefine(void)
{
   int jvar;

   for (jvar = 0 ; jvar < CODES ; ++jvar)
   if (All[jvar].selected)
   {
      if (All[jvar].twod)
      {
         if (RealSize == 8)
            All[jvar].NetVar = NetFile->add_var(All[jvar].Id,ncDouble,TimDim,LatDim,LonDim);
         else
            All[jvar].NetVar = NetFile->add_var(All[jvar].Id,ncFloat ,TimDim,LatDim,LonDim);
      }
      else
      {
         if (RealSize == 8)
            All[jvar].NetVar = NetFile->add_var(All[jvar].Id,ncDouble,TimDim,LevDim,LatDim,LonDim);
         else
            All[jvar].NetVar = NetFile->add_var(All[jvar].Id,ncFloat ,TimDim,LevDim,LatDim,LonDim);
      }
      All[jvar].NetVar->add_att("long_name"    ,All[jvar].Na);
      All[jvar].NetVar->add_att("standard_name",All[jvar].Na);
      All[jvar].NetVar->add_att("units"        ,All[jvar].Un);
      All[jvar].NetVar->add_att("code"         ,    jvar    );
      if (GaussianOutput)
      All[jvar].NetVar->add_att("grid_type"    ,"gaussian"  );
   }
}


void NetBuffer(double *d, float *f)
{
   for (int jdim = 0 ; jdim < DimGP ; ++jdim) *f++ = *d++;
}
   

void NetScale(float *f, int dim, double s)
{
   int j;

   for (j = 0 ; j < dim ; ++j) *f++ *= s;
}


void NetScale(double *f, int dim, double s)
{
   int j;

   for (j = 0 ; j < dim ; ++j) *f++ *= s;
}


void NetWrite32(int code, double *Var)
{
   int jlev;

   if (All[code].twod)
   {
      NetBuffer(Var,Record_float);    
      if (code==SLPCODE || code==PSCODE)
         NetScale(Record_float,Lats*(Lons+Cyclical),0.01);
      All[code].NetVar->set_cur(OutputCount);
      All[code].NetVar->put(Record_float,1,Lats,Lons+Cyclical);
   }
   else
   {
      for (jlev = 0 ; jlev < OutLevs ; ++jlev, Var += DimGP)
      {
         NetBuffer(Var,Record_float);
         if (code==WCODE)
            NetScale(Record_float,Lats*(Lons+Cyclical),0.01);
         All[code].NetVar->set_cur(OutputCount,jlev);
         All[code].NetVar->put(Record_float,1,1,Lats,Lons+Cyclical);
      }
   }
}


void NetWrite64(int code, double *Var)
{
   int jlev;

   if (All[code].twod)
   {
      memcpy(Record_double,Var,DimGP * sizeof(double));
      if (code==SLPCODE || code==PSCODE)
         NetScale(Record_double,Lats*(Lons+Cyclical),0.01);
      All[code].NetVar->set_cur(OutputCount);
      All[code].NetVar->put(Record_double,1,Lats,Lons+Cyclical);
   }
   else
   {
      for (jlev = 0 ; jlev < OutLevs ; ++jlev, Var += DimGP)
      {
         memcpy(Record_double,Var,DimGP * sizeof(double));
         if (code==WCODE)
            NetScale(Record_double,Lats*(Lons+Cyclical),0.01);
         All[code].NetVar->set_cur(OutputCount,jlev);
         All[code].NetVar->put(Record_double,1,1,Lats,Lons+Cyclical);
      }
   }
}


void NetWriteSection(int code, double *Var)
{
   int jvar,jlev;
   double *vp;

   for (jlev = 0 ; jlev < OutLevs ; ++jlev)
   {
      vp = Var + jlev * DimFC;
      if (code==SLPCODE || code==PSCODE || code==WCODE) 
      {
         for (jvar = 0 ; jvar < Lats ; ++jvar)
         {
            Record_float[jvar] = vp[jvar] * 0.01;
         }
      }
      else
      {
         for (jvar = 0 ; jvar < Lats ; ++jvar)
         {
            Record_float[jvar] = vp[jvar];
         }
      }
      All[code].NetVar->set_cur(OutputCount,jlev);
      All[code].NetVar->put(Record_float,1,1,Lats,1);
   }
}


void NetClose(void)
{
   int jt;
   double * Outtim;
   Outtim = new double[OutputCount];
   for (jt = 0 ; jt < OutputCount ; jt++) Outtim[jt] = jt;
   if (Mean == 0 && (DataStep < 0.99999 || DataStep > 1.00001))
   for (jt = 0 ; jt < OutputCount ; jt++) Outtim[jt] *= DataStep;

   TimVar->put(Outtim,OutputCount);
   delete Outtim;
   delete NetFile;
}

#endif


class ServiceGrid
{
   int HeadControl;
   FILE  *File;
   float *FloatBuffer;
public:
   int h4;  // head[4] = 1st dimension
   int h5;  // head[5] = 2nd dimension
   int Dim; // total dimension
   int FieldControl;
   ServiceGrid(FILE *, int, int);
  ~ServiceGrid(void);
   void Write(int *, double *);
   void WriteCode(int code, double *field, int twod);
   void Write_hspec(void);
   void Write_pspec(void);
   void Write_hfour(void);
   void Write_pfour(void);
   void Write_hgrid(void);
   void Write_pgrid(void);
   void Clear_hspec(void);
   void Clear_pspec(void);
   void Clear_hfour(void);
   void Clear_pfour(void);
   void Clear_hgrid(void);
   void Clear_pgrid(void);
};
   
ServiceGrid::ServiceGrid(FILE *fd, int Dim1, int Dim2)
{
   h4           = Dim1;
   h5           = Dim2;
   Dim          = Dim1 * Dim2;
   HeadControl  =   8 * RealSize;
   FieldControl = Dim * RealSize;
   File         = fd;
   FloatBuffer  = new float[Dim];
}

ServiceGrid::~ServiceGrid(void)
{
   delete FloatBuffer;
}

void ServiceGrid::Write(int *Head, double *Array)
{
   int i,j;
   LONG lhc,lfc;

   if (Debug) // Check for NaN (Not A Number)
   {
      for (j=0 ; j<Dim ; ++j)
      {
         if (isnan(Array[j]))
         {
            printf("\n Head: ");
            for (i=0 ; i<6 ; i++) printf(" %d",Head[i]);
            printf("\n Array[%d] = NaN\n",j);
            exit(1);
         }
      }
   }
   if (RealSize == 4)
   {
      for (j=0 ; j<Dim ; ++j) FloatBuffer[j] = Array[j];
   
      if (LongFCW)
      {
         lhc = HeadControl;
         lfc = FieldControl;
         fwrite(&lhc        ,sizeof(lhc)  ,1  ,File);
         fwrite( Head       ,sizeof(int)  ,8  ,File);
         fwrite(&lhc        ,sizeof(lhc)  ,1  ,File);
         fwrite(&lfc        ,sizeof(lfc)  ,1  ,File);
         fwrite( FloatBuffer,sizeof(float),Dim,File);
         fwrite(&lfc        ,sizeof(lfc)  ,1  ,File);
      }
      else
      {
         fwrite(&HeadControl ,sizeof(int)  ,1  ,File);
         fwrite( Head        ,sizeof(int)  ,8  ,File);
         fwrite(&HeadControl ,sizeof(int)  ,1  ,File);
         fwrite(&FieldControl,sizeof(int)  ,1  ,File);
         fwrite( FloatBuffer ,sizeof(float),Dim,File);
         fwrite(&FieldControl,sizeof(int)  ,1  ,File);
      }
   }
   else // RealSize == 8
   {
      if (LongFCW)
      {
         lhc = HeadControl;
         lfc = FieldControl;
         fwrite(&lhc        ,sizeof(lhc)   ,1  ,File);
         fwrite( Head       ,sizeof(int)   ,8  ,File);
         fwrite(&lhc        ,sizeof(lhc)   ,1  ,File);
         fwrite(&lfc        ,sizeof(lfc)   ,1  ,File);
         fwrite( Array      ,sizeof(double),Dim,File);
         fwrite(&lfc        ,sizeof(lfc)  ,1  ,File);
      }
      else
      {
         long LongHead[8];
         for (int i=0 ; i<8 ; ++i) LongHead[i] = Head[i];
         fwrite(&HeadControl ,sizeof(int)   ,1  ,File);
         fwrite( LongHead    ,RealSize      ,8  ,File);
         fwrite(&HeadControl ,sizeof(int)   ,1  ,File);
         fwrite(&FieldControl,sizeof(int)   ,1  ,File);
         fwrite( Array       ,RealSize      ,Dim,File);
         fwrite(&FieldControl,sizeof(int)   ,1  ,File);
      }
   }
}
   
void ServiceGrid::WriteCode(int code, double *field, int twod)
{
   int jlev;

   if (field == NULL)
   {
      fprintf(stderr, "*** Error in ServiceGrid::WriteCode\n");
      fprintf(stderr, "    Code %d not defined for this OutRep\n",code);
      fprintf(stderr, "    Skipping this record\n\n");
      return;
   }

   HeadOu[0] = code;
   HeadOu[4] = h4;
   HeadOu[5] = h5;
   if (twod)
   {
      HeadOu[1] = 0;
      HeadOu[7] = 0;
      Write(&HeadOu[0],field);
   }
   else for (jlev = 0; jlev < OutLevs; jlev++)
   {
      HeadOu[1] = LevelUnits[jlev];
      HeadOu[7] = SigmaF[jlev];
      Write(&HeadOu[0], field + jlev * Dim);
   }
}


void ServiceGrid::Write_hspec(void)
{
   for (int code = 0; code < CODES; code++)
      if (All[code].selected)
         WriteCode(code,&All[code].hsp[0],All[code].twod);
}


void ServiceGrid::Write_pspec(void)
{
   for (int code = 0; code < CODES; code++)
      if (All[code].selected)
         WriteCode(code,&All[code].psp[0],All[code].twod);
}


void ServiceGrid::Clear_hspec(void)
{
   for (int code = CODES-1 ; code >= 0 ; --code)
   if (All[code].hsp.size() && !All[code].twod)
   {
      if (Debug)
      {
         char tb[COLS];
         sprintf(tb,"Clear:   %p %6s[%3d].hsp",
                 &All[code].hsp[0],All[code].Id,code);
         LeftText(tb);
      }
      All[code].hsp.resize(0);  
   }
}

void ServiceGrid::Clear_hfour(void)
{
   for (int code = CODES-1 ; code >= 0 ; --code)
   if (All[code].hfc.size() && !All[code].twod)
   {
      if (Debug)
      {
         char tb[COLS];
         sprintf(tb,"CLear:   %p %6s[%3d].hfc",&All[code].hfc[0],All[code].Id,code);
         LeftText(tb);
      }
      All[code].hfc.resize(0);
   }
}

void ServiceGrid::Clear_hgrid(void)
{
   for (int code = CODES-1 ; code >= 0 ; --code)
   if (All[code].hgp.size())
   {
      if (Debug)
      {
         char tb[COLS];
         sprintf(tb,"CLear:   %p %6s[%3d].hgp",&All[code].hgp[0],All[code].Id,code);
         LeftText(tb);
      }
      All[code].hgp.resize(0);
   }
}

void ServiceGrid::Clear_pgrid(void)
{
   for (int code = CODES-1 ; code >= 0 ; --code)
   if (All[code].pgp.size())
   {
      if (Debug)
      {
         char tb[COLS];
         sprintf(tb,"CLear:   %p %6s[%3d].pgp",&All[code].pgp[0],All[code].Id,code);
         LeftText(tb);
      }
      All[code].pgp.resize(0);
   }
}

void ServiceGrid::Clear_pfour(void)
{
   for (int code = CODES-1 ; code >= 0 ; --code)
   if (All[code].pfc.size())
   {
      if (Debug)
      {
         char tb[COLS];
         sprintf(tb,"CLear:   %p %6s[%3d].pfc",&All[code].pfc[0],All[code].Id,code);
         LeftText(tb);
      }
      All[code].pfc.resize(0);
   }
}

void ServiceGrid::Clear_pspec(void)
{
   for (int code = CODES-1 ; code >= 0 ; --code)
   if (All[code].psp.size())
   {
      if (Debug)
      {
         char tb[COLS];
         sprintf(tb,"CLear:   %p %6s[%3d].psp",&All[code].psp[0],All[code].Id,code);
         LeftText(tb);
      }
      All[code].psp.resize(0);
   }
}

void ServiceGrid::Write_hfour(void)
{
   int code;

   for (code = 0; code < CODES; code++)
   if (All[code].selected)
      WriteCode(code,&All[code].hfc[0],All[code].twod);
}


void ServiceGrid::Write_pfour(void)
{
   int code;

   for (code = 0; code < CODES; code++)
   if (All[code].selected)
      WriteCode(code,&All[code].pfc[0],All[code].twod);
}


void ServiceGrid::Write_hgrid(void)
{
   int code;

   for (code = 0; code < CODES; code++)
   if (All[code].selected)
   {
      if (All[code].hgp.size() == 0)
      {
         fprintf(stderr, "*** Error in ServiceGrid::Write_hgrid\n");
         fprintf(stderr, "    Code %d is not available on sigma level\n",code);
         fprintf(stderr, "    Skipping this record\n\n");
         return;
      }
#ifdef NETCDF_OUTPUT
      if (NetCDF)
      {
         if (RealSize == 8) NetWrite64(code,&All[code].hgp[0]);
         else               NetWrite32(code,&All[code].hgp[0]);
      }
      else 
#endif
         WriteCode(code,&All[code].hgp[0],All[code].twod);
   }
   OutputCount++;
}


void ServiceGrid::Write_pgrid(void)
{
   int code;

   for (code = 0; code < CODES; code++)
   if (All[code].selected)
   {
#ifdef NETCDF_OUTPUT
      if (NetCDF)
      {
         if (RealSize == 8) NetWrite64(code,&All[code].pgp[0]);
         else               NetWrite32(code,&All[code].pgp[0]);
      }
      else 
#endif
         WriteCode(code,&All[code].pgp[0],All[code].twod);
   }
   OutputCount++;
}


class ServiceSect : public ServiceGrid
{
   double *Buffer;
public:
   ServiceSect(FILE *, int, int);
  ~ServiceSect(void);
   void WriteCode(int code, double *field, int twod);
   void Write_hfour(void);
   void Write_pfour(void);
};

ServiceSect::ServiceSect(FILE *fd, int Dim1, int Dim2) :
             ServiceGrid(fd,Dim1,Dim2)
{
   Buffer = new double[Dim1*Dim2];
}

ServiceSect::~ServiceSect(void)
{
   delete Buffer;
}

void ServiceSect::WriteCode(int code, double *field, int twod)
{
   int jlev;

   HeadOu[0] = code;
   HeadOu[1] =   -1;
   HeadOu[7] =    0;
   if (twod)
   {
      h5 = 1;
      memcpy(Buffer,field,Lats * sizeof(double));
   }
   else
   {
      h5 = OutLevs;
      for (jlev=0 ; jlev < OutLevs ; ++jlev)
         memcpy(Buffer+jlev*Lats,field+jlev*DimFC,Lats*sizeof(double));
   }
   HeadOu[4] = h4;
   HeadOu[5] = h5;
   Dim = h4 * h5;
   FieldControl = Dim * sizeof(float);
   Write(&HeadOu[0],Buffer);
}

void ServiceSect::Write_hfour(void)
{
   int code;

   for (code = 0; code < CODES; code++)
   if (All[code].selected)
   {
#ifdef NETCDF_OUTPUT
      if (NetCDF)
         NetWriteSection(code,&All[code].hfc[0]);
      else 
#endif
         WriteCode(code,&All[code].hfc[0],All[code].twod);
   }
   OutputCount++;
}


void ServiceSect::Write_pfour(void)
{
   int code;

   for (code = 0; code < CODES; code++)
   if (All[code].selected)
   {
#ifdef NETCDF_OUTPUT
      if (NetCDF)
         NetWriteSection(code,&All[code].pfc[0]);
      else 
#endif
         WriteCode(code,&All[code].pfc[0],All[code].twod);
   }
   OutputCount++;
}


ServiceGrid *HybSpec;
ServiceGrid *HybFour;
ServiceGrid *HybGrid;
ServiceSect *HybSect;

void HeadToDate(vector<int>Head, struct tm *Date)
{      
   Date->tm_year = Head[2] / 10000;
   Date->tm_mon  = Head[2] / 100 % 100;
   Date->tm_mday = Head[2] % 100;
   Date->tm_hour = Head[3] / 100;
   Date->tm_min  = Head[3] % 100;
}


/* ============================================== */
/* DoubleZero -  Set array of type double to zero */
/* ============================================== */

inline void DoubleZero(double *field, int words)
{
   memset((char *)field,0,words * sizeof(double));
}

/* ====================== */
/* Fast Fourier Transform */
/* ====================== */

/* constants in math.h :
#define M_E         2.71828182845904523536
#define M_LOG2E     1.44269504088896340736
#define M_LOG10E    0.434294481903251827651
#define M_LN2       0.693147180559945309417
#define M_LN10      2.30258509299404568402
#define M_PI        3.14159265358979323846
#define M_PI_2      1.57079632679489661923
#define M_PI_4      0.785398163397448309616
#define M_1_PI      0.318309886183790671538
#define M_2_PI      0.636619772367581343076
#define M_1_SQRTPI  0.564189583547756286948
#define M_2_SQRTPI  1.12837916709551257390
#define M_SQRT2     1.41421356237309504880
#define M_SQRT_2    0.707106781186547524401
*/

#define QUA 0.25
#define QT5 0.559016994374947

#define S36 0.587785252292473
#define S60 0.866025403784437
#define S72 0.951056516295154

#define SQ2 0.707106781186547524401

#define D60 (S60+S60)

#define FORK for(k=la;k<=kstop;k+=la){
#define LOOP for(l=0;l<la;++l){i=ibase;j=jbase;for(ijk=0;ijk<lot;++ijk){
#define ENDL i+=inc3;j+=inc4;}ibase++;jbase++;}

void fft_set(int n)
{
   int j,k,nfax;
   nfax = 0;
   for (k = 0; k < 9; ++k) ifax[k] = 0;
   ifax[9] = n;
   if    (n % 8 == 0) {ifax[++nfax] = 8; n /= 8;}
   while (n % 6 == 0) {ifax[++nfax] = 6; n /= 6;}
   while (n % 5 == 0) {ifax[++nfax] = 5; n /= 5;}
   while (n % 4 == 0) {ifax[++nfax] = 4; n /= 4;}
   while (n % 3 == 0) {ifax[++nfax] = 3; n /= 3;}
   if    (n % 2 == 0) {ifax[++nfax] = 2; n /= 2;}
   ifax[0] = nfax;
   for (k = 0; k < nfax/2; ++k)
   {
      j = ifax[k+1];
      ifax[k+1] = ifax[nfax-k];
      ifax[nfax-k] = j;
   }
}


int rpassc(double *a, double *b, double *c, double *d,
      int inc3, int inc4,
      int lot , int n   , int ifac, int la  )
{
/*
   rpassc' - performs one pass through data as part;
   of multiple real fft (fourier synthesis) routine;

   a is first real input vector
   b is equivalent to a + la
   c is first real output vector
   d is equivalent to c + ifac * la
   inc3 is the increment between input vectors a;
   inc4 is the increment between output vectors c;
   lot is the number of vectors;
   n is the length of the vectors;
   ifac is the current factor of n;
   la is the product of previous factors;
   ierr is an error indicator:;
   0 - pass completed without error;
   2 - ifac not catered for;
   3 - ifac only catered for if la=n/ifac;
*/

  int i0,i1,i2,i3,i4,i5,i6,i7;
  int j0,j1,j2,j3,j4,j5,j6,j7;
  int ia,ib,ic,id,ie,iF;
  int ja,jb,jc,jd,je,jf;
  int i,j,k,ijk,l,m;
  int ibase,jbase;
  int iink,jink;
  int jump;
  int kstop;

  double c1,c2,c3,c4,c5;
  double s1,s2,s3,s4,s5;
  double kpidn;
  double angle;
  double qqrt5;
  double ssin36;
  double ssin72;
  double pin;

  double a10,a11,a20,a21;
  double b10,b11,b20,b21;

  m     = n  / ifac;
  iink  = la;
  jink  = la;
  jump  = (ifac-1) * jink;
  kstop = (n-ifac) / (2*ifac);
  pin   = 2.0 * M_PI / n;
  ibase = 0;
  jbase = 0;

  switch (ifac)
  {
    case 2:
    {
      double a0m1,b0p1;

      i0 = j0 = 0;
      i1 = i0 + (m+m-la);
      j1 = j0 + jink;
      if (la != m)
      {
        LOOP
          c[j0+j] = a[i0+i] + a[i1+i];
          c[j1+j] = a[i0+i] - a[i1+i];
        ENDL
        i0    += iink;
        iink  += iink;
        i1    -= iink;
        ibase  = 0;
        jbase += jump;
        jump  += jump + jink;

        if (i0 != i1)
        {
          FORK
            angle = k * pin;
            c1 = cos(angle);
            s1 = sin(angle);
            ibase = 0;
            LOOP
              a0m1 = a[i0+i] - a[i1+i];
              b0p1 = b[i0+i] + b[i1+i];

              c[j0+j] = a[i0+i] + a[i1+i];
              d[j0+j] = b[i0+i] - b[i1+i];
              c[j1+j] = c1 * a0m1 - s1 * b0p1;
              d[j1+j] = s1 * a0m1 + c1 * b0p1;
            ENDL
            i0    += iink;
            i1    -= iink;
            jbase += jump;
          } /* End FORK */
          if (i0 > i1) return 0;
        } /* End (i0 != i1) */
        ibase = 0;
        LOOP
          c[j0+j] =  a[i0+i];
          c[j1+j] = -b[i0+i];
        ENDL
           }
           else /* (la != m) */ {
        LOOP
          c[j0+j] = 2.0 * (a[i0+i] + a[i1+i]);
          c[j1+j] = 2.0 * (a[i0+i] - a[i1+i]);
        ENDL
      }
      return 0;
    }

    case 3: {
      double afa1,a1p2,a1m2,a0mm,a0mp;
      double bfa1,b1p2,b1m2,b0mm,b0mp;

      i0 = j0 = 0 ;
      i1 = i0 + (m+m-la);
      i2 = i1;
      j1 = j0 + jink;
      j2 = j1 + jink;

      if (la != m) {
        LOOP
          afa1 = a[i0+i] - 0.5 * a[i1+i];
          bfa1 =           S60 * b[i1+i];

          c[j0+j] = a[i0+i] + a[i1+i];
          c[j1+j] = afa1 - bfa1;
          c[j2+j] = afa1 + bfa1;
        ENDL
        i0    += iink;
        iink  += iink;
        i1    += iink;
        i2    -= iink;
        jbase += jump;
        jump  += jump + jink;

        if (i0 != i2) {
          FORK
            angle = k * pin;
            c1 = cos(angle);
            s1 = sin(angle);
            angle += angle;
            c2 = cos(angle);
            s2 = sin(angle);
            ibase = 0;
            LOOP
              a1p2 = a[i0+i] - 0.5 * (a[i1+i] + a[i2+i]);
              b1m2 = b[i0+i] - 0.5 * (b[i1+i] - b[i2+i]);
              a1m2 =           S60 * (a[i1+i] - a[i2+i]);
              b1p2 =           S60 * (b[i1+i] + b[i2+i]);

              a0mm = a1p2 - b1p2;
              a0mp = a1p2 + b1p2;
              b0mm = b1m2 - a1m2;
              b0mp = b1m2 + a1m2;

              c[j0+j] = a[i0+i] + a[i1+i] + a[i2+i];
              d[j0+j] = b[i0+i] + b[i1+i] - b[i2+i];
              c[j1+j] = c1 * a0mm - s1 * b0mp;
              d[j1+j] = s1 * a0mm + c1 * b0mp;
              c[j2+j] = c2 * a0mp - s2 * b0mm;
              d[j2+j] = s2 * a0mp + c2 * b0mm;
            ENDL
            i0    += iink;
            i1    += iink;
            i2    -= iink;
            jbase += jump;
          } /* End FORK */
          if (i0 > i2) return 0;
        } /* End (i0 != i2) */
        ibase=0;
        LOOP
          a0mp = 0.5 * a[i0+i];
          b0mp = S60 * b[i0+i];

          c[j0+j] = a[i0+i] + a[i1+i];
          c[j1+j] = a0mp - a[i1+i] - b0mp;
          c[j2+j] = a[i1+i] - a0mp - b0mp;
        ENDL
      }
      else /* (la != m) */ {
        LOOP
          a0mp = 2.0 * a[i0+i] - a[i1+i];
          b0mp = D60 * b[i1+i];

          c[j0+j] = 2.0 * (a[i0+i] + a[i1+i]);
          c[j1+j] = a0mp - b0mp;
          c[j2+j] = a0mp + b0mp;
        ENDL
      }
      return 0;
    }

    case 4: {
      double a0m1,a0p2,a1p3,a0m2,a1m3,a0p2ma1p3,a0m2pb1p3,a0m2mb1p3;
      double b0p1,b0p2,b1p3,b0m2,b1m3,b0p2pa1m3,b0p2ma1m3,b0m2mb1m3;

      i0 = j0 = 0;
      i1 = i3 = i0 + (m+m-la);
      i2 = i1 + (m+m);
      j1 = j0 + jink;
      j2 = j1 + jink;
      j3 = j2 + jink;

      if (la != m) {
        LOOP
          a0p2 = a[i0+i] + a[i2+i];
          a0m2 = a[i0+i] - a[i2+i];

          c[j0+j] = a0p2 + a[i1+i];
          c[j1+j] = a0m2 - b[i1+i];
          c[j2+j] = a0p2 - a[i1+i];
          c[j3+j] = a0m2 + b[i1+i];
        ENDL
        i0    += iink;
        iink  += iink;
        i1    += iink;
        i2    -= iink;
        i3    -= iink;
        jbase += jump;
        jump  += jump + jink;

        if (i1 != i2) {
          FORK
            angle = kpidn = k * pin;
            c1 = cos(angle);
            s1 = sin(angle);
            angle += kpidn;
            c2 = cos(angle);
            s2 = sin(angle);
            angle += kpidn;
            c3 = cos(angle);
            s3 = sin(angle);
            ibase=0;
            LOOP
              a0p2 = a[i0+i] + a[i2+i];
              a0m2 = a[i0+i] - a[i2+i];
              a1p3 = a[i1+i] + a[i3+i];
              a1m3 = a[i1+i] - a[i3+i];
              b0p2 = b[i0+i] + b[i2+i];
              b0m2 = b[i0+i] - b[i2+i];
              b1p3 = b[i1+i] + b[i3+i];
              b1m3 = b[i1+i] - b[i3+i];

              a0p2ma1p3 = a0p2 - a1p3;
              a0m2pb1p3 = a0m2 + b1p3;
              a0m2mb1p3 = a0m2 - b1p3;
              b0p2pa1m3 = b0p2 + a1m3;
              b0p2ma1m3 = b0p2 - a1m3;
              b0m2mb1m3 = b0m2 - b1m3;

              c[j0+j] = a0p2 + a1p3;
              d[j0+j] = b0m2 + b1m3;
              c[j2+j] = c2 * a0p2ma1p3 - s2 * b0m2mb1m3;
              d[j2+j] = s2 * a0p2ma1p3 + c2 * b0m2mb1m3;
              c[j1+j] = c1 * a0m2mb1p3 - s1 * b0p2pa1m3;
              d[j1+j] = s1 * a0m2mb1p3 + c1 * b0p2pa1m3;
              c[j3+j] = c3 * a0m2pb1p3 - s3 * b0p2ma1m3;
              d[j3+j] = s3 * a0m2pb1p3 + c3 * b0p2ma1m3;
            ENDL
            i0    += iink;
            i1    += iink;
            i2    -= iink;
            i3    -= iink;
            jbase += jump;
          } /* End FORK */
          if (i1 > i2) return 0;
        } /* End (i1 != i2) */
        ibase=0;
        LOOP
          a0m1 = a[i0+i] - a[i1+i];
          b0p1 = b[i0+i] + b[i1+i];

          c[j0+j] =  a[i0+i] + a[i1+i];
          c[j2+j] =  b[i1+i] - b[i0+i];

          c[j1+j] =  SQ2 * (a0m1 - b0p1);
          c[j3+j] = -SQ2 * (a0m1 + b0p1);
        ENDL
      }
      else /* (la != m) */ {
        LOOP
          a0p2 = a[i0+i] + a[i2+i];
          a0m2 = a[i0+i] - a[i2+i];

          c[j0+j] = 2.0 * (a0p2 + a[i1+i]);
          c[j1+j] = 2.0 * (a0m2 - b[i1+i]);
          c[j2+j] = 2.0 * (a0p2 - a[i1+i]);
          c[j3+j] = 2.0 * (a0m2 + b[i1+i]);
        ENDL
      }
      return 0;
    }

    case 5: {
      double a1p2,a1m2,a0mm,a0mp,b136,b172,b236,b272;

      i0 = j0 = 0;
      i1 = i4 = i0 + (m+m-la);
      i2 = i3 = i1 + (m+m);
      j1 = j0 + jink;
      j2 = j1 + jink;
      j3 = j2 + jink;
      j4 = j3 + jink;

      if (la != m) {
        LOOP
          a1p2 = QUA * (a[i1+i] + a[i2+i]);
          a1m2 = QT5 * (a[i1+i] - a[i2+i]);

          a0mp = a[i0+i] - a1p2 + a1m2;
          a0mm = a[i0+i] - a1p2 - a1m2;

          b136 = b[i1+i] * S36;
          b172 = b[i1+i] * S72;
          b236 = b[i2+i] * S36;
          b272 = b[i2+i] * S72;

          c[j0+j] = a[i0+i] + a[i1+i] + a[i2+i];
          c[j1+j] = a0mp - b172 - b236;
          c[j2+j] = a0mm - b136 + b272;
          c[j3+j] = a0mm + b136 - b272;
          c[j4+j] = a0mp + b172 + b236;
        ENDL
        i0    += iink;
        iink  += iink;
        i1    += iink;
        i2    += iink;
        i3    -= iink;
        i4    -= iink;
        jbase += jump;
        jump  += jump + jink;

        if (i1 != i3) {
          FORK
            angle = kpidn = k * pin;
            c1 = cos(angle);
            s1 = sin(angle);
            angle += kpidn;
            c2 = cos(angle);
            s2 = sin(angle);
            angle += kpidn;
            c3 = cos(angle);
            s3 = sin(angle);
            angle += kpidn;
            c4 = cos(angle);
            s4 = sin(angle);
            ibase=0;
            LOOP
              a10=(a[i0+i]-0.25*((a[i1+i]+a[i4+i])+(a[i2+i]+a[i3+i])))
              +QT5*((a[i1+i]+a[i4+i])-(a[i2+i]+a[i3+i]));
              a20=(a[i0+i]-0.25*((a[i1+i]+a[i4+i])+(a[i2+i]+a[i3+i])))
              -QT5*((a[i1+i]+a[i4+i])-(a[i2+i]+a[i3+i]));
              b10=(b[i0+i]-0.25*((b[i1+i]-b[i4+i])+(b[i2+i]-b[i3+i])))
              +QT5*((b[i1+i]-b[i4+i])-(b[i2+i]-b[i3+i]));
              b20=(b[i0+i]-0.25*((b[i1+i]-b[i4+i])+(b[i2+i]-b[i3+i])))
              -QT5*((b[i1+i]-b[i4+i])-(b[i2+i]-b[i3+i]));

              a11=S72*(b[i1+i]+b[i4+i])+S36*(b[i2+i]+b[i3+i]);
              a21=S36*(b[i1+i]+b[i4+i])-S72*(b[i2+i]+b[i3+i]);
              b11=S72*(a[i1+i]-a[i4+i])+S36*(a[i2+i]-a[i3+i]);
              b21=S36*(a[i1+i]-a[i4+i])-S72*(a[i2+i]-a[i3+i]);

              c[j0+j]=a[i0+i]+((a[i1+i]+a[i4+i])+(a[i2+i]+a[i3+i]));
              d[j0+j]=b[i0+i]+((b[i1+i]-b[i4+i])+(b[i2+i]-b[i3+i]));
              c[j1+j]=c1*(a10-a11)-s1*(b10+b11);
              d[j1+j]=s1*(a10-a11)+c1*(b10+b11);
              c[j4+j]=c4*(a10+a11)-s4*(b10-b11);
              d[j4+j]=s4*(a10+a11)+c4*(b10-b11);
              c[j2+j]=c2*(a20-a21)-s2*(b20+b21);
              d[j2+j]=s2*(a20-a21)+c2*(b20+b21);
              c[j3+j]=c3*(a20+a21)-s3*(b20-b21);
              d[j3+j]=s3*(a20+a21)+c3*(b20-b21);
            ENDL
            i0    += iink;
            i1    += iink;
            i2    += iink;
            i3    -= iink;
            i4    -= iink;
            jbase += jump;
          } /* End FORK */
          if (i1 > i3) return 0;
        } /* End (i1 != i3) */
        ibase=0;
        LOOP
          c[j0+j] = a[i0+i] + a[i1+i] + a[i2+i];
          c[j1+j] = (QT5  * (a[i0+i]-a[i1+i])
             + (0.25  * (a[i0+i]+a[i1+i])-a[i2+i]))
             - (S36 *  b[i0+i]+S72*b[i1+i]);
          c[j4+j] =-(QT5  * (a[i0+i]-a[i1+i])
             + (0.25  * (a[i0+i]+a[i1+i])-a[i2+i]))
             - (S36 *  b[i0+i]+S72*b[i1+i]);
          c[j2+j] = (QT5  * (a[i0+i]-a[i1+i])
             - (0.25  * (a[i0+i]+a[i1+i])-a[i2+i]))
             - (S72 *  b[i0+i]-S36*b[i1+i]);
          c[j3+j] =-(QT5  * (a[i0+i]-a[i1+i])
             - (0.25  * (a[i0+i]+a[i1+i])-a[i2+i]))
             - (S72 *  b[i0+i]-S36*b[i1+i]);
        ENDL
      }  else {
        qqrt5  = 2.0 * QT5 ;
        ssin36 = 2.0 * S36;
        ssin72 = 2.0 * S72;
        LOOP
          c[j0+j]= 2.0 *(a[i0+i]+a[i1+i]+a[i2+i]);
          c[j1+j]=(2.0 *(a[i0+i]-0.25*(a[i1+i]+a[i2+i]))
            +qqrt5*(a[i1+i]-a[i2+i]))-(ssin72*b[i1+i]+ssin36*b[i2+i]);
          c[j2+j]=(2.0 *(a[i0+i]-0.25*(a[i1+i]+a[i2+i]))
            -qqrt5*(a[i1+i]-a[i2+i]))-(ssin36*b[i1+i]-ssin72*b[i2+i]);
          c[j3+j]=(2.0 *(a[i0+i]-0.25*(a[i1+i]+a[i2+i]))
            -qqrt5*(a[i1+i]-a[i2+i]))+(ssin36*b[i1+i]-ssin72*b[i2+i]);
          c[j4+j]=(2.0 *(a[i0+i]-0.25*(a[i1+i]+a[i2+i]))
            +qqrt5*(a[i1+i]-a[i2+i]))+(ssin72*b[i1+i]+ssin36*b[i2+i]);
        ENDL
      }
      return 0;
    }

    case 6: {
      ia = 0;
      ib = ia+(2*m-la);
      ic = ib+2*m;
      id = ic+2*m;
      ie = ic;
      iF = ib;
      ja = 0;
      jb = ja+jink;
      jc = jb+jink;
      jd = jc+jink;
      je = jd+jink;
      jf = je+jink;

      if (la != m) /* go to 690 */ {
        LOOP
          c[ja+j] =  (a[ia+i]+a[id+i]) +    (a[ib+i]+a[ic+i]) ;
          c[jd+j] =  (a[ia+i]-a[id+i]) -    (a[ib+i]-a[ic+i]) ;
          c[jb+j] = ((a[ia+i]-a[id+i]) +0.5*(a[ib+i]-a[ic+i]))
                     - S60*(b[ib+i]+b[ic+i]);
          c[jf+j] = ((a[ia+i]-a[id+i]) +0.5*(a[ib+i]-a[ic+i]))
                     + S60*(b[ib+i]+b[ic+i]);
          c[jc+j] = ((a[ia+i]+a[id+i]) -0.5*(a[ib+i]+a[ic+i]))
                     - S60*(b[ib+i]-b[ic+i]);
          c[je+j] = ((a[ia+i]+a[id+i]) -0.5*(a[ib+i]+a[ic+i]))
                     + S60*(b[ib+i]-b[ic+i]);
        ENDL
        ia    += iink;
        iink  += iink;
        ib    += iink;
        ic    += iink;
        id    -= iink;
        ie    -= iink;
        iF    -= iink;
        jbase += jump;
        jump  += jump+jink;

        if (ic != id) /* go to 660 */ {
          FORK
            angle = kpidn = k * pin;
            c1 = cos(angle);
            s1 = sin(angle);
            angle += kpidn;
            c2 = cos(angle);
            s2 = sin(angle);
            angle += kpidn;
            c3 = cos(angle);
            s3 = sin(angle);
            angle += kpidn;
            c4 = cos(angle);
            s4 = sin(angle);
            angle += kpidn;
            c5 = cos(angle);
            s5 = sin(angle);
            ibase=0;
            LOOP
              a11 = a[ie+i] + a[ib+i] + a[ic+i] + a[iF+i];
              a20 = a[ia+i] + a[id+i] - 0.5 * a11;
              a21 = S60*((a[ie+i]+a[ib+i])-(a[ic+i]+a[iF+i]));
              b11 = b[ib+i] - b[ie+i] + b[ic+i] - b[iF+i];
              b20 = b[ia+i] - b[id+i] - 0.5 * b11;
              b21 = S60*((b[ib+i]-b[ie+i])-(b[ic+i]-b[iF+i]));

              c[ja+j] = a[ia+i] + a[id+i] + a11;
              d[ja+j] = b[ia+i] - b[id+i] + b11;
              c[jc+j] = c2*(a20-b21)-s2*(b20+a21);
              d[jc+j] = s2*(a20-b21)+c2*(b20+a21);
              c[je+j] = c4*(a20+b21)-s4*(b20-a21);
              d[je+j] = s4*(a20+b21)+c4*(b20-a21);

              a11 = (a[ie+i]-a[ib+i]) + (a[ic+i]-a[iF+i]);
              b11 = (b[ie+i]+b[ib+i]) - (b[ic+i]+b[iF+i]);
              a20 = (a[ia+i]-a[id+i]) - 0.5 * a11;
              a21 = S60 * ((a[ie+i]-a[ib+i]) - (a[ic+i]-a[iF+i]));
              b20 = (b[ia+i]+b[id+i]) + 0.5 * b11;
              b21 = S60 * ((b[ie+i]+b[ib+i]) + (b[ic+i]+b[iF+i]));

              c[jd+j] = c3*(a[ia+i] - a[id+i] + a11)
                 - s3*(b[ia+i] + b[id+i] - b11);
              d[jd+j] = s3*(a[ia+i] - a[id+i] + a11)
                 + c3*(b[ia+i] + b[id+i] - b11);
              c[jb+j] = c1*(a20-b21)-s1*(b20-a21);
              d[jb+j] = s1*(a20-b21)+c1*(b20-a21);
              c[jf+j] = c5*(a20+b21)-s5*(b20+a21);
              d[jf+j] = s5*(a20+b21)+c5*(b20+a21);
            ENDL
            ia    += iink;
            ib    += iink;
            ic    += iink;
            id    -= iink;
            ie    -= iink;
            iF    -= iink;
            jbase += jump;
          }
          if (ic > id) return 0;
        }
        ibase=0;
        LOOP
          c[ja+j]= a[ib+i] + (a[ia+i] + a[ic+i]);
          c[jd+j]= b[ib+i] - (b[ia+i] + b[ic+i]);
          c[jb+j]= (S60*(a[ia+i]-a[ic+i]))-(0.5*(b[ia+i]+b[ic+i])+b[ib+i]);
          c[jf+j]=-(S60*(a[ia+i]-a[ic+i]))-(0.5*(b[ia+i]+b[ic+i])+b[ib+i]);
          c[jc+j]=  S60*(b[ic+i]-b[ia+i]) +(0.5*(a[ia+i]+a[ic+i])-a[ib+i]);
          c[je+j]=  S60*(b[ic+i]-b[ia+i]) -(0.5*(a[ia+i]+a[ic+i])-a[ib+i]);
        ENDL
      }  else {
        LOOP
          c[ja+j]=(2.0*(a[ia+i]+a[id+i]))+(2.0*(a[ib+i]+a[ic+i]));
          c[jd+j]=(2.0*(a[ia+i]-a[id+i]))-(2.0*(a[ib+i]-a[ic+i]));
          c[jb+j]=(2.0*(a[ia+i]-a[id+i])+(a[ib+i]-a[ic+i]))
            -(D60*(b[ib+i]+b[ic+i]));
          c[jf+j]=(2.0*(a[ia+i]-a[id+i])+(a[ib+i]-a[ic+i]))
            +(D60*(b[ib+i]+b[ic+i]));
          c[jc+j]=(2.0*(a[ia+i]+a[id+i])-(a[ib+i]+a[ic+i]))
            -(D60*(b[ib+i]-b[ic+i]));
          c[je+j]=(2.0*(a[ia+i]+a[id+i])-(a[ib+i]+a[ic+i]))
            +(D60*(b[ib+i]-b[ic+i]));
        ENDL
      }
      return 0;
    }

    case 8: {
      double a0p7,a1p5,a2p6,p073,p074,p152;
      double a0m7,a1m5,a2m6,m073,m074,m152;

      if (la != m) return 3;
      i0  = 0;
      i1  = i0 + iink;
      i2  = i1 + iink;
      i3  = i2 + iink;
      i4  = i3 + iink;
      i5  = i4 + iink;
      i6  = i5 + iink;
      i7  = i6 + iink;
      j0  = 0;
      j1  = j0 + jink;
      j2  = j1 + jink;
      j3  = j2 + jink;
      j4  = j3 + jink;
      j5  = j4 + jink;
      j6  = j5 + jink;
      j7  = j6 + jink;

      LOOP
        a0p7 = a[i0+i] + a[i7+i];
        a0m7 = a[i0+i] - a[i7+i];
        a1p5 = a[i1+i] + a[i5+i];
        a1m5 = a[i1+i] - a[i5+i];
        a2p6 = a[i2+i] + a[i6+i];
        a2m6 = a[i2+i] - a[i6+i];

        p073 = a0p7 + a[i3+i];
        m073 = a0p7 - a[i3+i];

        p074 = 2.0 * (a0m7 + a[i4+i]);
        m074 = 2.0 * (a0m7 - a[i4+i]);

        p152 = M_SQRT2 * (a1m5 + a2p6);
        m152 = M_SQRT2 * (a1m5 - a2p6);

        c[j0+j] = 2.0 * (p073 + a1p5);
        c[j4+j] = 2.0 * (p073 - a1p5);
        c[j2+j] = 2.0 * (m073 - a2m6);
        c[j6+j] = 2.0 * (m073 + a2m6);

        c[j1+j] = m074 + m152;
        c[j5+j] = m074 - m152;
        c[j3+j] = p074 - p152;
        c[j7+j] = p074 + p152;
      ENDL
    }
  }
  return 0;
}

int qpassc(double *a, double *b, double *c, double *d,
           int inc3, int inc4,
           int lot , int n   , int ifac, int la  )
{
/*
     qpassc - performs one pass through data as part;
     of multiple real fft (fourier analysis) routine;

     a is first real input vector;
     b is equivalent to a + ifac * la
     c is first real output vector;
     d is equivalent to c + la
     inc3 is the increment between input vectors a;
     inc4 is the increment between output vectors c;
     lot is the number of vectors;
     n is the length of the vectors;
     ifac is the current factor of n;
     la is the product of previous factors;
     qpassc returns an error indicator:;
       0 - pass completed without error;
       2 - ifac not catered for;
       3 - ifac only catered for if la=n/ifac;
*/

  int i0,i1,i2,i3,i4,i5,i6,i7;
  int j0,j1,j2,j3,j4,j5,j6,j7;
  int ia,ib,ic;
  int ja,jb,jc;
  int i,j,k;
  int ibase,jbase;
  int iink,jink;
  int ijk;
  int jump;
  int kstop;
  int l;
  int m;

  double a0,a1,a2,a3;
  double b0,b1,b2,b3;
  double c1,c2,c3,c4,c5;
  double s1,s2,s3,s4,s5;
  double w,x,y,z;
  double angle,kpidn,pin;

  m     = n  / ifac;
  iink  = la;
  jink  = la;
  jump  = (ifac-1) * iink;
  kstop = (n-ifac) / (2*ifac);
  pin   = 2.0 * M_PI / n;
  ibase = 0;
  jbase = 0;

  switch (ifac) {
    case 2: {
      i0 = j0 = 0;
      i1 = i0 + iink;
      j1 = j0 + (m+m-la);
      if (la != m) {
        LOOP
          c[j0+j] = a[i0+i] + a[i1+i];
          c[j1+j] = a[i0+i] - a[i1+i];
        ENDL
        j0    += jink;
        jink  += jink;
        j1    -= jink;
        ibase += jump;
        jump  += jump + iink;

        if (j0 != j1) {
          FORK
            angle = k * pin;
            c1 = cos(angle);
            s1 = sin(angle);
            jbase = 0;
            LOOP
              c[j0+j] = a[i0+i] + c1 * a[i1+i] + s1 * b[i1+i];
              c[j1+j] = a[i0+i] - c1 * a[i1+i] - s1 * b[i1+i];
              d[j0+j] = c1 * b[i1+i] - s1 * a[i1+i] + b[i0+i];
              d[j1+j] = c1 * b[i1+i] - s1 * a[i1+i] - b[i0+i];
            ENDL
            j0    += jink;
            j1    -= jink;
            ibase += jump;
          } /* End FORK */
          if (j0 > j1) return 0;
        } /* End (i0 != i1) */
        jbase = 0;
        LOOP
          c[j0+j] =  a[i0+i];
          d[j1+j] = -a[i1+i];
        ENDL
      }
      else /* (la != m) */ {
        z = 1.0 / n;
        LOOP
          c[j0+j] = z * (a[i0+i] + a[i1+i]);
          c[j1+j] = z * (a[i0+i] - a[i1+i]);
        ENDL
      }
      return 0;
    }

    case 3: {
      ia = 0;
      ib = ia + iink;
      ic = ib + iink;

      ja = 0;
      jb = ja + (m+m-la);
      jc = jb;

      if (la != m) /* else 390 */ {
        LOOP
          c[ja+j] = a[ia+i] + a[ib+i] + a[ic+i];
          c[jb+j] = a[ia+i] - 0.5 * (a[ib+i] + a[ic+i]);
          d[jb+j] = S60 * (a[ic+i] - a[ib+i]);
        ENDL
        ja    += jink;
        jink  += jink;
        jb    += jink;
        jc    -= jink;
        ibase += jump;
        jump  += jump + iink;

        if (ja != jc) /* else  360 */ {
          FORK
            angle = k * pin;
            c1 = cos(angle);
            s1 = sin(angle);
            angle += angle;
            c2 = cos(angle);
            s2 = sin(angle);
            jbase = 0;
            LOOP
              a1 = c1 * a[ib+i] + s1 * b[ib+i] + c2 * a[ic+i] + s2 * b[ic+i];
              b1 = c1 * b[ib+i] - s1 * a[ib+i] + c2 * b[ic+i] - s2 * a[ic+i];
              a2 = a[ia+i] - 0.5 * a1;
              b2 = b[ia+i] - 0.5 * b1;
              a3 = S60 * (c1 * a[ib+i] + s1 * b[ib+i] - c2 * a[ic+i] - s2 * b[ic+i]);
              b3 = S60 * (c1 * b[ib+i] - s1 * a[ib+i] - c2 * b[ic+i] + s2 * a[ic+i]);

              c[ja+j] = a[ia+i] + a1;
              d[ja+j] = b[ia+i] + b1;
              c[jb+j] = a2 + b3;
              d[jb+j] = b2 - a3;
              c[jc+j] = a2 - b3;
              d[jc+j] =-b2 - a3;
            ENDL
            ja    += jink;
            jb    += jink;
            jc    -= jink;
            ibase += jump;
          } /* End FORK */
          if (ja > jc) return 0;
        } /* End (ia != ic) */
        jbase = 0;
        LOOP
        /* soweit */
          c[ja+j] = a[ia+i] + 0.5 * (a[ib+i] - a[ic+i]);
          d[ja+j] =-S60 * (a[ib+i] + a[ic+i]);
          c[jb+j] = a[ia+i] - a[ib+i] + a[ic+i];
        ENDL
           }
           else /* (la != m) */ {
        z = 1.0   / n;
        y = S60 / n;
        LOOP
          c[ja+j] = z * (a[ia+i] + a[ib+i] + a[ic+i]);
          c[jb+j] = z * (a[ia+i] - 0.5 * (a[ib+i] + a[ic+i]));
          d[jb+j] = y * (a[ic+i] - a[ib+i]);
        ENDL
      }
      return 0;
    }

    case 4: {
      double a0p2,a1p3;

      i0 = 0;
      i1 = i0 + iink;
      i2 = i1 + iink;
      i3 = i2 + iink;
      j0 = 0;
      j1 = j0 + (m+m-la);
      j2 = j1 + (m+m   );
      j3 = j1;

      if (la != m) /*else go to 490 */ {
        LOOP
          a0p2 = a[i0+i] + a[i2+i];
          a1p3 = a[i1+i] + a[i3+i];

          c[j0+j] = a0p2 + a1p3;
          c[j2+j] = a0p2 - a1p3;

          c[j1+j] = a[i0+i] - a[i2+i];
          d[j1+j] = a[i3+i] - a[i1+i];
        ENDL
        j0    += jink;
        jink  += jink;
        j1    += jink;
        j2    -= jink;
        j3    -= jink;
        ibase += jump;
        jump  += jump + iink;

        if (j1 != j2) /* else go to 460; */ {
          FORK
            angle = kpidn = k * pin;
            c1 = cos(angle);
            s1 = sin(angle);
            angle += kpidn;
            c2 = cos(angle);
            s2 = sin(angle);
            angle += kpidn;
            c3 = cos(angle);
            s3 = sin(angle);
            jbase=0;
            LOOP
              a0 = a[i0+i] + c2 * a[i2+i] + s2 * b[i2+i];
              a2 = a[i0+i] - c2 * a[i2+i] - s2 * b[i2+i];
              b0 = b[i0+i] + c2 * b[i2+i] - s2 * a[i2+i];
              b2 = b[i0+i] - c2 * b[i2+i] + s2 * a[i2+i];

              a1 = c1*a[i1+i] + s1*b[i1+i] + c3*a[i3+i] + s3*b[i3+i];
              a3 = c1*a[i1+i] + s1*b[i1+i] - c3*a[i3+i] - s3*b[i3+i];
              b1 = c1*b[i1+i] - s1*a[i1+i] + c3*b[i3+i] - s3*a[i3+i];
              b3 = c1*b[i1+i] - s1*a[i1+i] - c3*b[i3+i] + s3*a[i3+i];

              c[j0+j] = a0 + a1;
              c[j2+j] = a0 - a1;
              d[j0+j] = b0 + b1;
              d[j2+j] = b1 - b0;
              c[j1+j] = a2 + b3;
              c[j3+j] = a2 - b3;
              d[j1+j] = b2 - a3;
              d[j3+j] =-b2 - a3;
            ENDL
            j0    += jink;
            j1    += jink;
            j2    -= jink;
            j3    -= jink;
            ibase += jump;
          } /* End FORK */
          if (j1 > j2) return 0;
        } /* End (i1 != i2) */
        jbase=0;
        LOOP
          c[j0+j] = a[i0+i] + SQ2 * (a[i1+i] - a[i3+i]);
          c[j1+j] = a[i0+i] - SQ2 * (a[i1+i] - a[i3+i]);
          d[j0+j] =-a[i2+i] - SQ2 * (a[i1+i] + a[i3+i]);
          d[j1+j] = a[i2+i] - SQ2 * (a[i1+i] + a[i3+i]);
        ENDL
      }
      else /* (la != m) */ {
        z = 1.0 / n;
        LOOP
          a0p2 = a[i0+i] + a[i2+i];
          a1p3 = a[i1+i] + a[i3+i];

          c[j0+j] = z * (a0p2 + a1p3);
          c[j2+j] = z * (a0p2 - a1p3);
          c[j1+j] = z * (a[i0+i] - a[i2+i]);
          d[j1+j] = z * (a[i3+i] - a[i1+i]);
        ENDL
      }
      return 0;
    }

    case 5: {
      double a1p4,a2p3,b1p4,b2p3,a025,b025,asps,bsps,a0pq,b0pq;
      double a1m4,a2m3,b1m4,b2m3,aqrt,bqrt,asms,bsms,a0mq,b0mq;

      i0 = 0;
      i1 = i0 + iink;
      i2 = i1 + iink;
      i3 = i2 + iink;
      i4 = i3 + iink;
      j0 = 0;
      j1 = j0 + (m+m-la);
      j2 = j1 + (m+m);
      j3 = j2;
      j4 = j1;

      if (la != m) /* else go to 590; */ {
        LOOP
          a1p4 = a[i1+i] + a[i4+i];
          a1m4 = a[i1+i] - a[i4+i];
          a2p3 = a[i2+i] + a[i3+i];
          a2m3 = a[i2+i] - a[i3+i];

          a025 = a[i0+i] - 0.25 * (a1p4 + a2p3);
          aqrt =           QT5 * (a1p4 - a2p3);

          c[j0+j] = a[i0+i] + a1p4 + a2p3;
          c[j1+j] = a025 + aqrt;
          c[j2+j] = a025 - aqrt;
          d[j1+j] = -S72 * a1m4 - S36 * a2m3;
          d[j2+j] = -S36 * a1m4 + S72 * a2m3;
        ENDL
        j0    += jink;
        jink  += jink;
        j1    += jink;
        j2    += jink;
        j3    -= jink;
        j4    -= jink;
        ibase += jump;
        jump  += jump + iink;

        if (j1 != j3) {
          FORK
            angle = kpidn = k * pin;
            c1 = cos(angle);
            s1 = sin(angle);
            angle += kpidn;
            c2 = cos(angle);
            s2 = sin(angle);
            angle += kpidn;
            c3 = cos(angle);
            s3 = sin(angle);
            angle += kpidn;
            c4 = cos(angle);
            s4 = sin(angle);
            jbase=0;
            LOOP
              a1p4 = c1*a[i1+i] + s1*b[i1+i] + c4*a[i4+i] + s4*b[i4+i];
              a1m4 = c1*a[i1+i] + s1*b[i1+i] - c4*a[i4+i] - s4*b[i4+i];
              a2p3 = c2*a[i2+i] + s2*b[i2+i] + c3*a[i3+i] + s3*b[i3+i];
              a2m3 = c2*a[i2+i] + s2*b[i2+i] - c3*a[i3+i] - s3*b[i3+i];
              b1p4 = c1*b[i1+i] - s1*a[i1+i] + c4*b[i4+i] - s4*a[i4+i];
              b1m4 = c1*b[i1+i] - s1*a[i1+i] - c4*b[i4+i] + s4*a[i4+i];
              b2p3 = c2*b[i2+i] - s2*a[i2+i] + c3*b[i3+i] - s3*a[i3+i];
              b2m3 = c2*b[i2+i] - s2*a[i2+i] - c3*b[i3+i] + s3*a[i3+i];

              a025 = a[i0+i] - 0.25 * (a1p4 + a2p3);
              aqrt =           QT5  * (a1p4 - a2p3);
              b025 = b[i0+i] - 0.25 * (b1p4 + b2p3);
              bqrt =           QT5  * (b1p4 - b2p3);

              a0pq = a025 + aqrt;
              a0mq = a025 - aqrt;
              b0pq = b025 + bqrt;
              b0mq = b025 - bqrt;

              asps = S72 * a1m4 + S36 * a2m3;
              asms = S36 * a1m4 - S72 * a2m3;
              bsps = S72 * b1m4 + S36 * b2m3;
              bsms = S36 * b1m4 - S72 * b2m3;

              c[j0+j] = a[i0+i] + a1p4 + a2p3;
              c[j1+j] = a0pq + bsps;
              c[j2+j] = a0mq + bsms;
              c[j3+j] = a0mq - bsms;
              c[j4+j] = a0pq - bsps;
              d[j0+j] = b[i0+i] + b1p4 + b2p3;
              d[j1+j] = b0pq - asps;
              d[j2+j] = b0mq - asms;
              d[j3+j] =-b0mq - asms;
              d[j4+j] =-b0pq - asps;
            ENDL
            j0    += jink;
            j1    += jink;
            j2    += jink;
            j3    -= jink;
            j4    -= jink;
            ibase += jump;
          } /* End FORK */
          if (j1 > j3) return 0;
        } /* End (jb != jd) */
        jbase=0;
        LOOP
          a1p4 = a[i1+i] + a[i4+i];
          a1m4 = a[i1+i] - a[i4+i];
          a2p3 = a[i2+i] + a[i3+i];
          a2m3 = a[i2+i] - a[i3+i];

          a025 = a[i0+i] + 0.25 * (a1m4 - a2m3);
          aqrt =           QT5 * (a1m4 + a2m3);

          c[j0+j] = a025 + aqrt;
          c[j1+j] = a025 - aqrt;
          c[j2+j] = a[i0+i] - a1m4 + a2m3;
          d[j0+j] = -S36 * a1p4 - S72 * a2p3;
          d[j1+j] = -S72 * a1p4 + S36 * a2p3;

        ENDL
      }  else {
        z = 1.0   / n;
        y = QT5  / n;
        x = S36 / n;
        w = S72 / n;

        LOOP
          a1p4 = a[i1+i] + a[i4+i];
          a1m4 = a[i1+i] - a[i4+i];
          a2p3 = a[i2+i] + a[i3+i];
          a2m3 = a[i2+i] - a[i3+i];

          a025 = z * (a[i0+i] - 0.25 * (a1p4 + a2p3));
          aqrt = y * (a1p4 - a2p3);

          c[j0+j] = z * (a[i0+i] + a1p4 + a2p3);
          c[j1+j] = a025 + aqrt;
          c[j2+j] = a025 - aqrt;
          d[j1+j] = -w * a1m4 - x * a2m3;
          d[j2+j] =  w * a2m3 - x * a1m4;
        ENDL
      }
      return 0;
    }

    case 6: {
      double ab1a,ab2a,ab3a,ab4a,ab5a;
      double ab1b,ab2b,ab3b,ab4b,ab5b;
      double a0p3,a1p4,a1p5,a2p4,a2p5;
      double a0m3,a1m4,a1m5,a2m4,a2m5;
      double b1p4,b2p5;
      double b1m4,b2m5;
      double ap05,bp05,ap60,bp60;
      double am05,bm05,am60,bm60;

      i0 = 0;
      i1 = i0 + iink;
      i2 = i1 + iink;
      i3 = i2 + iink;
      i4 = i3 + iink;
      i5 = i4 + iink;
      j0 = 0;
      j1 = j0 + (m+m-la);
      j2 = j1 + (m+m);
      j3 = j2 + (m+m);
      j4 = j2;
      j5 = j1;

      if (la != m) {
        LOOP
          a0p3 = a[i0+i] + a[i3+i];
          a0m3 = a[i0+i] - a[i3+i];
          a1p4 = a[i1+i] + a[i4+i];
          a1m4 = a[i1+i] - a[i4+i];
          a2p5 = a[i2+i] + a[i5+i];
          a2m5 = a[i2+i] - a[i5+i];

          c[j0+j] = a0p3 + a1p4 + a2p5;
          c[j3+j] = a0m3 + a2m5 - a1m4;

          c[j1+j] = a0m3 - 0.5 * (a2m5 - a1m4);
          c[j2+j] = a0p3 - 0.5 * (a1p4 + a2p5);

          d[j1+j] = S60 * (-a2m5 - a1m4);
          d[j2+j] = S60 * ( a2p5 - a1p4);
        ENDL
        j0    += jink;
        jink  += jink;
        j1    += jink;
        j2    += jink;
        j3    -= jink;
        j4    -= jink;
        j5    -= jink;
        ibase += jump;
        jump  += jump+iink;

        if (j2 != j3) {
          FORK
            angle = kpidn = k * pin;
            c1 = cos(angle);
            s1 = sin(angle);
            angle += kpidn;
            c2 = cos(angle);
            s2 = sin(angle);
            angle += kpidn;
            c3 = cos(angle);
            s3 = sin(angle);
            angle += kpidn;
            c4 = cos(angle);
            s4 = sin(angle);
            angle += kpidn;
            c5 = cos(angle);
            s5 = sin(angle);
            jbase = 0;
            LOOP
              ab1a = c1 * a[i1+i] + s1 * b[i1+i];
              ab1b = c1 * b[i1+i] - s1 * a[i1+i];
              ab2a = c2 * a[i2+i] + s2 * b[i2+i];
              ab2b = c2 * b[i2+i] - s2 * a[i2+i];
              ab3a = c3 * a[i3+i] + s3 * b[i3+i];
              ab3b = c3 * b[i3+i] - s3 * a[i3+i];
              ab4a = c4 * a[i4+i] + s4 * b[i4+i];
              ab4b = c4 * b[i4+i] - s4 * a[i4+i];
              ab5a = c5 * a[i5+i] + s5 * b[i5+i];
              ab5b = c5 * b[i5+i] - s5 * a[i5+i];

              a1p4 = ab1a + ab4a;
              a1m4 = ab1a - ab4a;
              a2p5 = ab2a + ab5a;
              a2m5 = ab2a - ab5a;

              b1p4 = ab1b + ab4b;
              b1m4 = ab1b - ab4b;
              b2p5 = ab2b + ab5b;
              b2m5 = ab2b - ab5b;

              ap05 = a[i0+i] + ab3a - 0.5 * (a1p4 + a2p5);
              bp05 = b[i0+i] + ab3b - 0.5 * (b1p4 + b2p5);
              am05 = a[i0+i] - ab3a - 0.5 * (a2m5 - a1m4);
              bm05 =-b[i0+i] + ab3b - 0.5 * (b1m4 - b2m5);

              ap60 = S60 * ( a2p5 - a1p4);
              bp60 = S60 * ( b2p5 - b1p4);
              am60 = S60 * (-a2m5 - a1m4);
              bm60 = S60 * (-b2m5 - b1m4);

              c[j0+j] = a[i0+i] + ab3a + a1p4 + a2p5;
              d[j0+j] = b[i0+i] + ab3b + b1p4 + b2p5;
              c[j1+j] = am05 - bm60;
              d[j1+j] = am60 - bm05;
              c[j2+j] = ap05 - bp60;
              d[j2+j] = ap60 + bp05;
              c[j3+j] = a[i0+i] - ab3a - a1m4 + a2m5;
              d[j3+j] =-b[i0+i] + ab3b + b1m4 - b2m5;
              c[j4+j] = ap05 + bp60;
              d[j4+j] = ap60 - bp05;
              c[j5+j] = am05 + bm60;
              d[j5+j] = am60 + bm05;
            ENDL
            j0    += jink;
            j1    += jink;
            j2    += jink;
            j3    -= jink;
            j4    -= jink;
            j5    -= jink;
            ibase += jump;
          }
          if (j2 > j3) return 0;
        }
        jbase = 0;
        LOOP
          a1p5 = a[i1+i] + a[i5+i];
          a1m5 = a[i1+i] - a[i5+i];
          a2p4 = a[i2+i] + a[i4+i];
          a2m4 = a[i2+i] - a[i4+i];

          c[j0+j] = a[i0+i] + 0.5 * a2m4 + S60 * a1m5;
          d[j0+j] =-a[i3+i] - 0.5 * a1p5 - S60 * a2p4;
          c[j1+j] = a[i0+i] - a2m4;
          d[j1+j] = a[i3+i] - a1p5;
          c[j2+j] = a[i0+i] + 0.5 * a2m4 - S60 * a1m5;
          d[j2+j] =-a[i3+i] - 0.5 * a1p5 + S60 * a2p4;
        ENDL
      }  else {
        z = 1.0 / n;
        y = S60 / n;
        LOOP
          a0p3 = a[i0+i] + a[i3+i];
          a0m3 = a[i0+i] - a[i3+i];
          a1p4 = a[i1+i] + a[i4+i];
          a1m4 = a[i1+i] - a[i4+i];
          a2p5 = a[i2+i] + a[i5+i];
          a2m5 = a[i2+i] - a[i5+i];

          c[j0+j] = z * (a0p3 + a1p4 + a2p5);
          c[j3+j] = z * (a0m3 + a2m5 - a1m4);

          c[j1+j] = z * (a0m3 - 0.5 * (a2m5 - a1m4));
          c[j2+j] = z * (a0p3 - 0.5 * (a1p4 + a2p5));

          d[j1+j] = y * (-a2m5 - a1m4);
          d[j2+j] = y * ( a2p5 - a1p4);
        ENDL
      }
      return 0;
    }

    case 8: {
      double a0p4,a1p5,a2p6,a3p7;
      double a0m4,a1m5,a2m6,a3m7;

      if (la != m) return 3;
      i0 = 0;
      i1 = i0 + iink;
      i2 = i1 + iink;
      i3 = i2 + iink;
      i4 = i3 + iink;
      i5 = i4 + iink;
      i6 = i5 + iink;
      i7 = i6 + iink;
      j0 = 0;
      j1 = j0 + jink;
      j2 = j1 + jink;
      j3 = j2 + jink;
      j4 = j3 + jink;
      j5 = j4 + jink;
      j6 = j5 + jink;
      j7 = j6 + jink;
      z  = 1.0      / n;
      y  = SQ2 / n;

      LOOP
        a0p4 = a[i0+i] + a[i4+i];
        a0m4 = a[i0+i] - a[i4+i];
        a1p5 = a[i1+i] + a[i5+i];
        a1m5 = a[i1+i] - a[i5+i];
        a2p6 = a[i2+i] + a[i6+i];
        a2m6 = a[i2+i] - a[i6+i];
        a3p7 = a[i3+i] + a[i7+i];
        a3m7 = a[i3+i] - a[i7+i];

        c[j0+j] = z * (a0p4 + a1p5 + a2p6 + a3p7);
        c[j7+j] = z * (a0p4 - a1p5 + a2p6 - a3p7);

        c[j3+j] = z * (a0p4 - a2p6);
        c[j4+j] = z * (a3p7 - a1p5);

        c[j1+j] = z * a0m4 + y * (a1m5 - a3m7);
        c[j5+j] = z * a0m4 - y * (a1m5 - a3m7);
        c[j2+j] =-z * a2m6 - y * (a1m5 + a3m7);
        c[j6+j] = z * a2m6 - y * (a1m5 + a3m7);
      ENDL
    }
  }
  return 0;
}

void fc2gp(double *fc, double *gp, int Lat, int Lon, int Lev, int Fou)
{
   int Lot;
   int fou;
   int ia;
   int ifac;
   int j;
   int jump;
   int k;
   int la;
   int lat;
   int lev;
   int lon;
   int nfax;
   int rix;
   int wix;

   double *wpt;

/* fc2gp performs fourier to gridpoint transforms using       */
/* multiple fast fourier transform of length Lon              */
/*                                                            */
/* fc  - real array of fourier coefficients fc[Lev][Fou][Lat] */
/* gp  - real array of gridpoints           gp[Lev][Lat][Lon] */
/* Lat - Number of latitudes                                  */
/* Lon - Number of longitudes                                 */
/* Lev - Number of levels                                     */
/* Fou - Number of fourier coefficients on 1 latitude         */

/* x(j) = sum(k=0,...,n-1)(c(k)*exp(2*i*j*k*pi/Lon))          */
/*        where c(k) = a(k) + i*b(k) and c(n-k) = a(k)-i*b(k) */

   jump = (Lon + 2) | 1;
   Lot  = Lev * Lat;
   nfax = ifax[0];

   for (lev = 0; lev < Lev; ++lev)
   {
      for (lat = 0; lat < Lat; ++lat)
      {
         wix = jump * (lat + lev * Lat);
         rix = lat  +  lev * Lat * Fou;
         for (fou = 0  ; fou < Fou && fou < jump ; ++fou)
            wfc[wix + fou] = fc[rix + fou * Lat];
         for (fou = Fou; fou < jump; ++fou)
            wfc[wix + fou] = 0.0;
         wfc[wix+1] = 0.5 * wfc[wix];
      }
   }

   ia = 1;
   la = 1;

   for (k = 0; k < nfax; ++k) {
      ifac = ifax[k+1];
      if (k&1) rpassc(wgp,wgp+la,wfc+ia,wfc+ia+ifac*la,
                      jump,jump,Lot,Lon,ifac,la);
      else     rpassc(wfc+ia,wfc+ia+la,wgp,wgp+ifac*la,
                      jump,jump,Lot,Lon,ifac,la);
      la *= ifac;
      ia  = 0;
   }

   if (nfax & 1) wpt = wgp;
   else          wpt = wfc;

   for (j = 0; j < Lot ; ++j)
      for (lon = 0; lon < Lon; ++lon)
         gp[lon + j * Lon] = wpt[lon + j * jump];
}

void gp2fc(double *gp, double *fc, int Lat, int Lon, int Lev, int Fou)
{
   int Lot;
   int fou;
   int ia;
   int ifac;
   int jump;
   int k;
   int la;
   int lat;
   int lev;
   int lon;
   int lot;
   int nfax;
   int rix;
   int wix;

   double *wpt;

/* gp2fc performs gridpoint to fourier transforms using       */
/* multiple fast fourier transform of length Lon              */
/*                                                            */
/* fc  - real array of fourier coefficients fc[Lev][Fou][Lat] */
/* gp  - real array of gridpoints           gp[Lev][Lat][Lon] */
/* Lat - Number of latitudes                                  */
/* Lon - Number of longitudes                                 */
/* Lev - Number of levels                                     */
/* Fou - Number of fourier coefficients on 1 latitude         */

/* a(k) =  (1/n) * sum(j=0,...,n-1)(x(j) * cos(2*j*k*pi/n))   */
/* b(k) = -(1/n) * sum(j=0,...,n-1)(x(j) * sin(2*j*k*pi/n))   */

   if (!gp) Abort("gp2fc called with NULL pointer for gp");
   if (!fc) Abort("gp2fc called with NULL pointer for fc");
   jump = (Lon + 2) | 1;
   Lot  = Lev * Lat;
   nfax = ifax[0];

   rix = 0;
   wix = 0;
   for (lot = 0; lot < Lot; ++lot) {
      for (lon = 0; lon < Lon; ++lon) wgp[wix+lon] = gp[rix+lon];
      wgp[wix+Lon  ] = 0.0;
      wgp[wix+Lon+1] = 0.0;
      rix += Lon;
      wix += jump;
   }

   ia = 0;
   la = Lon;

   for (k = 0; k < nfax; ++k) {
      ifac = ifax[nfax-k];
      la  /= ifac;
      if (k & 1) qpassc(wfc,wfc+ifac*la,wgp+ia,wgp+ia+la,
                        jump,jump,Lot,Lon,ifac,la);
      else       qpassc(wgp+ia,wgp+ia+ifac*la,wfc,wfc+la,
                        jump,jump,Lot,Lon,ifac,la);
      ia = 1;
   }

   if (nfax & 1) wpt = wfc;
   else          wpt = wgp+1;

   for (lev = 0; lev < Lev; ++lev) {
      for (lat = 0; lat < Lat; ++lat) {
         rix = jump * (lat + lev * Lat);
         wix = lat + lev * Lat * Fou;
         fc[wix    ] = wpt[rix];
         fc[wix+Lat] = 0.0;
         for (fou = 2  ; fou < Fou  ; ++fou)
            fc[wix + fou * Lat] = wpt[rix + fou - 1];
      }
   }
}


inline void SwapIEEE(char W[4])
{
   char B;

   B = W[0]; W[0] = W[3]; W[3] = B;
   B = W[1]; W[1] = W[2]; W[2] = B;
}

inline void SwapIEEE64(char W[8])
{
   char B;

   B = W[0]; W[0] = W[3]; W[3] = B;
   B = W[1]; W[1] = W[2]; W[2] = B;
   B = W[4]; W[4] = W[7]; W[7] = B;
   B = W[5]; W[5] = W[6]; W[6] = B;
}


int check_fcw(int fcws, int fcwe)
{
   if (fcwe != fcws)
   {
      printf("\n*************** ERROR reading input file **************\n");
      printf("The FORTRAN control words (FCW) of a record don't match\n");
      printf("FCW before record = %d\n",fcws);
      printf("FCW after  record = %d\n",fcwe);
      printf("File position: %ld\n",ftell(fpi));
      printf("Possible causes are:\n");
      printf("1) Model crashed leaving an incomplete output file\n");
      printf("2) Corrupt data file (cache or disk problems)\n");
      return 1;
   }
   return 0;
}


/* =============================== */
/* Read IEEE 32 bit data from file */
/* =============================== */


inline int ReadINT(void)
{
   int k;
   fread(&k,sizeof(k),1,fpi);
   if (Endian) SwapIEEE((char *)&k);
   return k;
}

inline LONG ReadLONG(void)
{
   LONG l;
   fread(&l,sizeof(l),1,fpi);
   if (Endian) SwapIEEE64((char *)&l);
   return l;
}

inline int ReadFCW(void)
{
   int k;
   LONG l;
   if (LongFCW)
   {
      fread(&l,sizeof(l),1,fpi);
      if (Endian) SwapIEEE64((char *)&l);
      k = (int)l;
   }
   else
   {
      fread(&k,sizeof(k),1,fpi);
      if (Endian) SwapIEEE((char *)&k);
   }
   return k;
}

inline float ReadFLOAT(void)
{
   int i;
   float f;
   i = fread(&f,sizeof(f),1,fpi);
   if (i != 1) Abort("Unexpected EOF in ReadFLOAT");
   if (Endian) SwapIEEE((char *)&f);
   return f;
}

inline double ReadDOUBLE(void)
{
   int i;
   double f;
   i = fread(&f,sizeof(f),1,fpi);
   if (i != 1) Abort("Unexpected EOF in ReadDOUBLE");
   if (Endian) SwapIEEE64((char *)&f);
   return f;
}

inline int ReadINTRecord(void)
{
   int k,b,e;
   b = ReadFCW();
   fread(&k,sizeof(k),1,fpi); /* IEEE 32-bit integer */
   e = ReadFCW();
   if (check_fcw(b,e)) Abort("record control word mismatch in ReadINTRecord");
   if (Endian) SwapIEEE((char *)&k);
   return k;
}

inline int Skip_FORTRAN_Record(void)
{
    int fcws,fcwe;

    fcws = ReadFCW();
    if (feof(fpi)) return 0;
    fseek(fpi,fcws,SEEK_CUR);
    fcwe = ReadFCW();
    if (check_fcw(fcws,fcwe))
       Abort("record control word mismatch in Skip_FORTRAN_Record");
    return fcws;
}

inline void Swap_FORTRAN_Record(char *c, int n)
{
   char b;
   int  i;

   for (i=0 ; i < n ; i+=4)
   {
      b = c[i  ]; c[i  ] = c[i+3]; c[i+3] = b;
      b = c[i+1]; c[i+1] = c[i+2]; c[i+2] = b;
   }
}

inline void Swap_FORTRAN_Double(char *c, int n)
{
   char b;
   int  i;

   for (i=0 ; i < n ; i+=8)
   {
      b = c[i  ]; c[i  ] = c[i+7]; c[i+7] = b;
      b = c[i+1]; c[i+1] = c[i+6]; c[i+6] = b;
      b = c[i+2]; c[i+2] = c[i+5]; c[i+5] = b;
      b = c[i+3]; c[i+3] = c[i+4]; c[i+4] = b;
   }
}

inline int Read_FORTRAN_Record(void)
{
    int fcws,fcwe;

    fcws = ReadFCW();
    if (feof(fpi)) return 0;
    fread(Record_char,1,fcws,fpi);
    fcwe = ReadFCW();
    if (check_fcw(fcws,fcwe)) Abort("record control word mismatch in Read_FORTRAN_Record");
    if (Endian) Swap_FORTRAN_Record(Record_char,fcws);
    return fcws;
}

inline int Read_FORTRAN_Double_Record(void)
{
    int fcws,fcwe;

    fcws = ReadFCW();
    if (feof(fpi)) return 0;
    fread(Record_char,1,fcws,fpi);
    fcwe = ReadFCW();
    if (check_fcw(fcws,fcwe)) Abort("record control word mismatch in Read_FORTRAN_Double_Record");
    if (Endian) Swap_FORTRAN_Double(Record_char,fcws);
    return fcws;
}

int ReadHeaderRecord(void)
{
   int h,i,fcws,fcwe;
   
   fcws = ReadFCW();
   if (feof(fpi)) return 1;
   if (fcws == 8) /* Skip PUMA header */
   {
      if (Debug)
      {
         printf("Skipping %d header words\n",HeaderWords);
         for (i=0 ; i < HeaderWords ; ++i)
         {
            h = ReadINT();
            printf("HW[%2d] = %8x %d\n",i,h,h);
         }
      }
      else
      for (i=0 ; i < HeaderWords ; ++i) ReadINT();
      fcws = ReadFCW();
      if (feof(fpi))
      {
         printf("\n#### Empty data file #####\n");
         printf("Mark [Write Ouput] in MoSt\n");
         printf("or set NOUTPUT=1 in file <puma_namelist>\n");
         Abort("Empty data file");
      }
   }
   if (fcws != HeadSize)
   {
       printf("FCW = %d (should be %d)\n",fcws,HeadSize);
       Abort("Wrong FCW found in ReadHeaderRecord");
   }
   for (i=0 ; i < 8 ; ++i)
   {
      if (HeadSize == 32) HeadIn[i] = ReadINT();
      else                HeadIn[i] = ReadLONG();
   }
   fcwe = ReadFCW();
   if (check_fcw(fcws,fcwe)) Abort("FCW mismatch in ReadHeaderRecord");
   return 0;
}
 
 
void DecodePumaHeader(void)
{
   PumaCode        = HeadIn[0];
   PumaLevel       = HeadIn[1];
   NewDate.tm_year = HeadIn[2] / 10000;
   NewDate.tm_mon  = HeadIn[2] / 100 % 100;
   NewDate.tm_mday = HeadIn[2] % 100;
   NewDate.tm_hour = HeadIn[3] / 100;
   NewDate.tm_min  = HeadIn[3] % 100;
   if (DayDivisor > 1) NewMonth = (TermCount / DPM) % 12 + 1;
   else                NewMonth = NewDate.tm_mon;
   if (HeadIn[4] * HeadIn[5] == DimSP) RepGrib = REP_SPECTRAL;
   else                                RepGrib = REP_GAUSS;

   if (PumaCode < 0 || PumaCode >= CODES)
   {
      printf("Illegal Code %d in header\n",PumaCode);
      Abort("Code < 0 or Code > CODES found");
   }
   All[PumaCode].detected = TRUE;
}

// =============
// ReadPumaArray
// =============

void ReadPumaArray(double *Array)
{
   int i,j,fcws;
   double  zfac;
   double  zmin;

   if (RealSize == sizeof(float)) fcws = Read_FORTRAN_Record();
   else                           fcws = Read_FORTRAN_Double_Record();
   if (fcws == sizeof(float) * (Truncation + 2))  //   Packed spectral data
   {
      for (i=0 ; i <= Truncation ; i++)
      {
         Array[2*i  ] = Record_float[i];
         Array[2*i+1] = 0.0;    // Imaginary parts of zonal modes
      }
      zfac   = 1.0 / Record_float[Truncation+1];
      fcws = Read_FORTRAN_Record();
      if (CoreBigEndian)
      {
         for (i=2*Truncation+2,j=0 ; i < DimSP ; ++i,++j)
         {
            Array[i] = (Record_short[j] - 32768) * zfac;
         }
      }
      else
      {
         for (i=2*Truncation+2,j=0 ; i < DimSP ; i+=2,j+=2)
         {
            Array[i  ] = (Record_short[j+1] - 32768) * zfac;
            Array[i+1] = (Record_short[j  ] - 32768) * zfac;
         }
      }
   }
   else if (fcws == sizeof(float) * DimSP) // Unpacked spectral float data
   {
      for (i=0 ; i < DimSP ; ++i) Array[i] = Record_float[i];
   }
   else if (fcws == sizeof(double) * DimSP) // Unpacked spectral double data
   {
      memcpy(Array,Record_double,fcws);
   }
   else if (fcws == sizeof(float) * DimGG)  // Unpacked grid float data
   {
      for (i=0 ; i < DimGG ; ++i) Array[i] = Record_float[i];
   }
   else if (fcws == sizeof(double) * DimGG) // Unpacked grid double data
   {
      memcpy(Array,Record_double,fcws);
   }
   else if (fcws == 8)               /* Packed grid data */
   {
      zmin = Record_float[0];
      zfac = 1.0 / Record_float[1];
      fcws = Read_FORTRAN_Record();
      if (CoreBigEndian)
      {
         for (i=0 ; i < DimGG ; ++i)
         {
            Array[i] = Record_short[i] * zfac + zmin;
         }
      }
      else
      {
         for (i=0 ; i < DimGG ; i+=2)
         {
            Array[i  ] = Record_short[i+1] * zfac + zmin;
            Array[i+1] = Record_short[i  ] * zfac + zmin;
         }
      }
   }
   else Abort("fcws error in ReadPumaArray");
}


/* ============= */
/* SkipPumaArray */
/* ============= */

void SkipPumaArray(void)
{
   int fcws;
   fcws = Skip_FORTRAN_Record();
   if (fcws == 8 || fcws == 4 * (Truncation + 2))
      fcws = Skip_FORTRAN_Record();
}


/* ============================================= */
/* phcs - Compute values of Legendre polynomials */
/*        and their meridional derivatives       */
/* ============================================= */

void phcs(double *PNM, double *HNM, int Waves, double PMU,
     double *ZTEMP1, double *ZTEMP2)
{
   long TwoWaves;

   long JK;
   long JN;
   long JM;

   double JNmJK;
   double ZCOS2;
   double Lat;
   double ZAN;
   double ZSINPAR;
   double ZCOSPAR;
   double ZSQP;
   double ZCOSFAK;
   double ZSINFAK;
   double ZQ;
   double ZWM2;
   double ZW;
   double ZWQ;
   double ZQ2M1;
   double ZWM2Q2;
   double Z2Q2;
   double ZCNM;
   double ZDNM;
   double ZENM;

   TwoWaves  = Waves << 1;
   ZCOS2     = sqrt(1.0 - PMU * PMU);
   Lat       = acos(PMU);
   ZAN       = 1.0;
   ZTEMP1[0] = 0.5;

   for (JN = 1; JN < TwoWaves; JN++) {
      ZSQP    = 1.0 / sqrt((double)(JN+JN*JN));
      ZAN    *= sqrt(1.0 - 1.0 / (4 * JN * JN));

      ZCOSPAR = cos(Lat * JN);
      ZSINPAR = sin(Lat * JN) * JN * ZSQP;
      ZCOSFAK = 1.0;

      for (JK = 2; JK < JN; JK += 2) {
         JNmJK = JN - JK;
         ZCOSFAK *= (JK-1) * (JN+JNmJK+2) / (JK * (JN+JNmJK+1));
         ZSINFAK  = ZCOSFAK * (JNmJK) * ZSQP;
         ZCOSPAR += ZCOSFAK * cos(Lat * JNmJK);
         ZSINPAR += ZSINFAK * sin(Lat * JNmJK);
      }

      /*  Code for JK == JN */

      if ((JN & 1) == 0) {
      ZCOSFAK *= (double)((JN-1) * (JN+2)) / (double)(JN * (JN+1));
      ZCOSPAR += ZCOSFAK * 0.5;
      }
      ZTEMP1[JN  ] = ZAN * ZCOSPAR;
      ZTEMP2[JN-1] = ZAN * ZSINPAR;
   }

   memcpy(PNM,ZTEMP1,Waves * sizeof(double));
   PNM += Waves;
   memcpy(PNM,ZTEMP2,Waves * sizeof(double));
   PNM += Waves;

   HNM[0] = 0.0;
   for (JN = 1; JN < Waves; JN++) HNM[JN] =
      JN * (PMU * ZTEMP1[JN] - sqrt((JN+JN+1.0) / (JN+JN-1.0)) * ZTEMP1[JN-1]);
   HNM += Waves;

   HNM[0] = PMU * ZTEMP2[0];
   for (JN = 1; JN < Waves; JN++)
      HNM[JN] = (JN+1) * PMU * ZTEMP2[JN]
              - sqrt((double)((JN+JN+3) * ((JN+1) * (JN+1) - 1))
              / (double)(JN+JN+1)) * ZTEMP2[JN-1];
   HNM += Waves;

   for (JM = 2; JM < Waves; JM++) {
      PNM[0] = sqrt(1.0 + 1.0 / (JM+JM)) * ZCOS2 * ZTEMP2[0];
      HNM[0] = JM * PMU * PNM[0];
      for (JN = 1; JN < TwoWaves-JM; JN++) {
          ZQ      = JM + JM + JN - 1;
          ZWM2    = ZQ+JN;
          ZW      = ZWM2+2;
          ZWQ     = ZW*ZQ;
          ZQ2M1   = ZQ*ZQ-1.;
          ZWM2Q2  = ZWM2*ZQ2M1;
          Z2Q2    = ZQ2M1*2;
          ZCNM    = sqrt((ZWQ*(ZQ-2.))/(ZWM2Q2-Z2Q2));
          ZDNM    = sqrt((ZWQ*(JN+1.))/ZWM2Q2);
          ZENM    = sqrt(ZW * JN /((ZQ+1.0) * ZWM2));
          PNM[JN] = ZCNM * ZTEMP1[JN] - PMU
                  * (ZDNM * ZTEMP1[JN+1] - ZENM * PNM[JN-1]);
          HNM[JN] = (JM + JN) * PMU * PNM[JN]
                  - sqrt(ZW * JN * (ZQ+1) / ZWM2) * PNM[JN-1];
      }
      memcpy(ZTEMP1,ZTEMP2,TwoWaves * sizeof(double));
      memcpy(ZTEMP2,PNM   ,TwoWaves * sizeof(double));
      PNM += Waves;
      HNM += Waves;
   }
}

void legini(void)
{
   int jlat,jm,jn,jz;
   int jsp,pdim,hdim;
   double *hnm,*pnm,*ZTEMP1,*ZTEMP2,gmusq;
   double znn1,zgmu;
   char tb[COLS+2];
   double poliv,pol2v,pliuv,plivv;

   hdim = 2 * Waves * Waves;

   if (PolyCreate) /* Generate filenames for polynomials */
   {
      sprintf(polin,"b6poli.T%d",Truncation);
      sprintf(pol2n,"b6pol2.T%d",Truncation);
      sprintf(pliun,"b6pliu.T%d",Truncation);
      sprintf(plivn,"b6pliv.T%d",Truncation);
   }
   else if (PolyDisk) /* Generate filenames for polynomials */
   {
      sprintf(polin,"/comm/T1365/b6poli.t%d",Truncation);
      sprintf(pol2n,"/comm/T1365/b6pol2.t%d",Truncation);
      sprintf(pliun,"/comm/T1365/b6pliu.t%d",Truncation);
      sprintf(plivn,"/comm/T1365/b6pliv.t%d",Truncation);
      polif = fopen(polin,"r");
      pol2f = fopen(pol2n,"r");
      pliuf = fopen(pliun,"r");
      plivf = fopen(plivn,"r");

      if (polif && pol2f && pliuf && plivf)
         sprintf(tb,"Legendre Polynomials read from disk");
      else
      {
         sprintf(tb,"Legendre Polynomials calculated");
         PolyDisk = 0;
      }
      CenterText(tb);
   }

   if (PolyDisk) pdim = Lats;
   else          pdim = Lats * DimSP_half;
   poli   = new double[pdim];
   pol2   = new double[pdim];
   pliu   = new double[pdim];
   pliv   = new double[pdim];
   pnm    = new double[hdim];
   hnm    = new double[hdim];
   ZTEMP1 = new double[Fouriers];
   ZTEMP2 = new double[Fouriers];

   // if gridtype for output is not selected, choose Gauss grid
   // for matching resolution and regular grid else

   if (GaussGrid < 0) GaussGrid = (Lats == Gats && Lons == Gons);

   Gaulat = new GauLat(Gats,"Gaulat"); // Gaussian latitudes of input grid
   Outlon = new RegLon(Lons,"Outlon"); // Regular longitudes of output grid
   if (GaussGrid) Outlat = new GauLat(Lats,"Outlat");
   else           Outlat = new RegLat(Lats,"Outlat");

   if (Debug)
   {
      Gaulat->PrintArray();
      if (Lats != Gats || !GaussGrid) Outlat->PrintArray();
      Outlon->PrintArray();
   }

   if (PolyCreate)
   {
      polif = fopen(polin,"w");
      pol2f = fopen(pol2n,"w");
      pliuf = fopen(pliun,"w");
      plivf = fopen(plivn,"w");
      for (jlat = 0 ; jlat < Lats ; ++jlat)
      {
         gmusq = 1.0 - Outlat->gmu[jlat] * Outlat->gmu[jlat];
         zgmu  = sqrt(gmusq);
         phcs(pnm,hnm,Waves,Outlat->gmu[jlat],ZTEMP1,ZTEMP2);
         for (jm = 0; jm < Waves; jm++)
         {
            for (jn = 0; jn < Waves - jm ; jn++)
            {
               jz = jm+jn;
               znn1 = 0.0;
               if (jz > 0) znn1 = 1.0 / (jz * (jz+1)); 
               poliv = pnm[jm*Waves+jn] * 2.0;
               fwrite(&poliv,sizeof(double),1,polif);
               pol2v = hnm[jm*Waves+jn] / PlanetRadius;
               fwrite(&pol2v,sizeof(double),1,pol2f);
               pliuv = pnm[jm*Waves+jn] * 2.0 * PlanetRadius * znn1 * jm / zgmu;
               fwrite(&pliuv,sizeof(double),1,pliuf);
               plivv = hnm[jm*Waves+jn] * 2.0 * PlanetRadius * znn1 / zgmu;
               fwrite(&plivv,sizeof(double),1,plivf);
            }
         }
      }
   }
   else if (PolyDisk)
   {
      for (jlat = 0 ; jlat < Lats ; ++jlat)
      {
         gmusq = 1.0 - Outlat->gmu[jlat] * Outlat->gmu[jlat];
         CosPhi[jlat] =  sqrt(gmusq);
         RevCosPhi[jlat] = 1.0 / CosPhi[jlat];
         DerivationFactor[jlat] = RevCosPhi[jlat] / PlanetRadius;
      }
   }
   else   /* Normal computation of polynomials */
   {
      for (jlat = 0 ; jlat < Lats ; ++jlat)
      {
         gmusq = 1.0 - Outlat->gmu[jlat] * Outlat->gmu[jlat];
         CosPhi[jlat] =  sqrt(gmusq);
         RevCosPhi[jlat] = 1.0 / CosPhi[jlat];
         DerivationFactor[jlat] = RevCosPhi[jlat] / PlanetRadius;
         phcs(pnm,hnm,Waves,Outlat->gmu[jlat],ZTEMP1,ZTEMP2);
         jsp = jlat;
         for (jm = 0; jm < Waves; jm++)
         {
            for (jn = 0; jn < Waves - jm ; jn++)
            {
               jz = jm+jn;
               znn1 = 0.0;
               if (jz > 0) znn1 = 1.0 / (jz * (jz+1)); 
               poli[jsp] = pnm[jm*Waves+jn] * 2.0;
               pol2[jsp] = hnm[jm*Waves+jn] / PlanetRadius;
               pliu[jsp] = pnm[jm*Waves+jn] * 2.0 * PlanetRadius * znn1 * jm / sqrt(gmusq);
               pliv[jsp] = hnm[jm*Waves+jn] * 2.0 * PlanetRadius * znn1 / sqrt(gmusq);
               jsp += Lats;
            }
         }
      }
   }

   delete [] pnm;
   delete [] hnm;
   delete [] ZTEMP1;
   delete [] ZTEMP2;
}


void spvfc(double *sd, double *sz, double *fu, double *fv, int klev,int nlat,int nfc,int nt)
{
   int lev,jmm,jfc,lat;
   double sdr,sdi,szr,szi;
   double *fur,*fui,*fvr,*fvi;
   double *poa,*pob;

   DoubleZero(fu,klev*nlat*nfc);
   DoubleZero(fv,klev*nlat*nfc);

   for (lev = 0; lev < klev; lev++)
   {
      if (PolyDisk)
      {
         rewind(pliuf);
         rewind(plivf);
      }
      poa = pliu;
      pob = pliv;
      for (jmm = 0; jmm <= nt; jmm++)
      {
         for (jfc = jmm; jfc <= nt; jfc++)
         {
            sdr = *sd++;
            sdi = *sd++;
            szr = *sz++;
            szi = *sz++;
            fur = fu       ;
            fui = fu + nlat;
            fvr = fv       ;
            fvi = fv + nlat;
            if (PolyDisk)
            {
               fread(poa=pliu,sizeof(double),Lats,pliuf);
               fread(pob=pliv,sizeof(double),Lats,plivf);
            }
            for (lat = 0; lat < nlat; lat++)
            {
               *fur += -*pob * szr + *poa * sdi;
               *fui += -*pob * szi - *poa * sdr;
               *fvr +=  *poa * szi + *pob * sdr;
               *fvi += -*poa * szr + *pob * sdi;
               fur++;
               fui++;
               fvr++;
               fvi++;
               poa++;
               pob++;
            }
         }
         fu += 2 * nlat;
         fv += 2 * nlat;
      }
   }
}


// ========================================
// sp2fci - Inverse Legendre Transformation
// ========================================

void sp2fci(double *sa,double *fa,int klev)
{
   int lev,m,n;
   double sar,sai;
   double *Far,*Fai,*pol;
   DoubleZero(fa,klev*DimFC);
   for (lev = 0; lev < klev; ++lev)
   {
      if (PolyDisk) rewind(polif);
      pol = poli;
      for (n = 0; n <= Truncation; ++n)
      {
         if (PolyDisk) fread(pol=poli,sizeof(double),Lats,polif);
         sar = *sa;
         sa += 2;
         for (Far=fa; Far < fa+Lats;++Far,++pol)
         {
            *Far += *pol * sar;
         }
      }
      fa += 2 * Lats;

      for (m = 1; m <= Truncation; ++m)
      {
         for (n = m; n <= Truncation; ++n)
         {
            if (PolyDisk) fread(pol=poli,sizeof(double),Lats,polif);
            sar = *sa++     ;
            sai = *sa++     ;
            for (Far=fa,Fai=fa+Lats; Far<fa+Lats; ++Far,++Fai,++pol)
            {
               *Far += *pol * sar;
               *Fai += *pol * sai;
            }
         }
         fa += 2 * Lats;
      }
   }
}

void sp2fcd(double *sa,double *fa,int klev,int nlat,int nfc,int nt)
{
   int lev,jmm,jfc,lat;
   double sar,sai;
   double *Far,*fai,*pol;
   double zpo;
   DoubleZero(fa,klev*nlat*nfc);
   for (lev = 0; lev < klev; lev++)
   {
      pol = pol2;
      if (PolyDisk) rewind(pol2f);
      for (jmm = 0; jmm <= nt; jmm++)
      {
         for (jfc = jmm; jfc <= nt; jfc++)
         {
            sar = *sa++     ;
            sai = *sa++     ;
            Far = fa        ;
            fai = fa + nlat ;
            if (PolyDisk) fread(pol=pol2,sizeof(double),Lats,pol2f);
            for (lat = 0; lat < nlat; lat++)
            {
               zpo = -2.0 * *pol * RevCosPhi[lat];
               *Far += zpo * sar;
               *fai += zpo * sai;
               Far++;
               fai++;
               pol++;
            }
         }
         fa += 2 * nlat;
      }
   }
}

void fc2sp(double *fa, double *sa, int klev, int nlat, int nt)
{
   int lev,jmm,jfc,lat;
   double sar,sai,*Far,*fai,*pol;
   double zpo;
   for (lev = 0; lev < klev; lev++)
   {
      pol = poli;
      if (PolyDisk) rewind(polif);
      for (jmm = 0; jmm <= nt; jmm++)
      {
         for (jfc = jmm; jfc <= nt; jfc++)
         {
            Far = fa        ;
            fai = fa + nlat ;
            sar = 0.0       ;
            sai = 0.0       ;
            if (PolyDisk) fread(pol=poli,sizeof(double),Lats,polif);
            for (lat = 0; lat < nlat; lat++)
            {
               zpo  = *pol * 0.5 * Outlat->gwt[lat];
               sar += zpo * *Far;
               sai += zpo * *fai;
               Far++;
               fai++;
               pol++;
            }
            *sa++ = sar;
            *sa++ = sai;
         }
         fa += 2 * nlat;
      }
   }
}

void OMEGA(void)
{
   int i,j;
   double DeltaHybrid;
   double Cterm;
   double Pterm;
   double *omega = &Omega->hgp[0];
   double *diver = &Divergence->hgp[0];
   double *halfp = &HalfPress->hgp[0];
   double *fullp = &FullPress->hgp[0];
   double *uwind = &u_wind->hgp[0];
   double *vwind = &v_wind->hgp[0];

/* Compute half level part of vertical velocity */

   for (i = 0; i < DimGP ; i++) omega[i] = 0.0;
   for (j = 0; j < SigLevs; j++) {
      DeltaHybrid = vct[SigLevs+j+2] - vct[SigLevs+j+1];
      for (i = 0; i < DimGP; i++) {
        omega[DimGP] = *omega
                      - *diver * (halfp[DimGP] - *halfp) - DeltaHybrid
                      * (*uwind * dpsdx->hgp[i] 
                      +  *vwind * dpsdy->hgp[i]);
        omega++;
        halfp++;
        diver++;
        uwind++;
        vwind++;
      }
   }

/* interpolate to full levels  */

   omega = &Omega->hgp[0];
   for (i = 0; i < Dim3GP; i++)
      omega[i] = 0.5 * (omega[i] + omega[i+DimGP]);

/* compute full level part of vertical velocity */

   omega = &Omega->hgp[0];
   halfp = &HalfPress->hgp[0];
   fullp = &FullPress->hgp[0];
   uwind = &u_wind->hgp[0];
   vwind = &v_wind->hgp[0];

   for (j = 0; j < SigLevs; j++) {
      DeltaHybrid = vct[SigLevs+j+2] - vct[SigLevs+j+1];
      if (DeltaHybrid) {
         Cterm = vct[j+1] * vct[SigLevs+j+1] - vct[j] * vct[SigLevs+j+2];
         for (i = 0; i < DimGP; i++) {
            if (Cterm != 0.0) Pterm = Cterm /
               (halfp[DimGP] - *halfp) * log(halfp[DimGP] / *halfp);
            else Pterm = 0.0;

           *omega += *fullp *
              (*uwind * dpsdx->hgp[i] + *vwind * dpsdy->hgp[i])
              / (halfp[DimGP] - *halfp) * (DeltaHybrid + Pterm);
           omega++;
           halfp++;
           fullp++;
           uwind++;
           vwind++;
         }
      }
      else {
         omega += DimGP;
         halfp += DimGP;
         fullp += DimGP;
         uwind += DimGP;
         vwind += DimGP;
      }
   }
}

void Omega_w(double w[], double om[], double T[], double P[])
{
   int i;

   for (i=0 ; i < Dim3GP ; ++i)
   {
      w[i] = -om[i] * RD * T[i] / (Grav * P[i]);
   }

}

void Extrap(double *slp, double *aph, double *apf,
            double *Geopotential, double *t, int nhor)
{
   double zrg    = 1.0 / Grav;
   double alpha,tstar,tmsl,ZPRT,ZPRTAL;

   int j;

   for (j = 0; j < nhor; ++j) {
      if (Geopotential[j] < 0.0001 && Geopotential[j] > -0.0001) slp[j] = aph[j];
      else {
         alpha = RD * RLAPSE * zrg;
         tstar = (1.0 + alpha * (aph[j]/apf[j] - 1.0)) * t[j];
         if (tstar < 255.0) tstar = 0.5 * (255.0 + tstar);
         tmsl = tstar + RLAPSE * zrg * Geopotential[j];
         if (tmsl > 290.5 && tstar > 290.5) {
            tstar = 0.5 * (290.5 + tstar);
            tmsl  = tstar;
         }
         if (tmsl-tstar < 0.000001 && tstar-tmsl < 0.000001) alpha = 0.0;
         else if (Geopotential[j] > 0.0001 || Geopotential[j] < -0.0001)
            alpha = RD * (tmsl-tstar) / Geopotential[j];
         ZPRT   = Geopotential[j] / (RD * tstar);
         ZPRTAL = ZPRT * alpha;
         slp[j] = aph[j] * exp(ZPRT*(1.0-ZPRTAL*(0.5-ZPRTAL/3.0)));
      }
   }
}

double ExtraT(double PRES, double APH, double APF, double GEOS, double T)
{
   double zrg    = 1.0 / Grav;
   double tstar,ztsz,Z1,ZTMSL,ZALPH,PEVAL,ZHTS,ZALP;

   tstar = (1.0 + RLAPSE * RD * zrg * (APH/APF - 1.0)) * T;
   ztsz   = tstar;
   Z1     = tstar + RLAPSE * zrg * GEOS;
   if (tstar < 255.0) tstar = 0.5 * (255.0 + tstar);
   ZTMSL = tstar + RLAPSE * zrg * GEOS;
   if (ZTMSL > 290.5 && tstar > 290.5)
   {
      tstar = 0.5 * (290.5 + tstar);
      ZTMSL  = tstar;
   }
   if (ZTMSL > 290.5 && tstar <= 290.5) ZTMSL=290.5;
   ZALPH=RD*RLAPSE*zrg;
   if ( ZTMSL-tstar < 0.000001 && tstar-ZTMSL < 0.000001) ZALPH=0.0;
   if ((ZTMSL-tstar > 0.000001 || tstar-ZTMSL > 0.000001) &&
       (GEOS > 0.0001 || GEOS < -0.0001))
      ZALPH=RD*(ZTMSL-tstar)/GEOS;
   if (PRES <= APH) PEVAL = ((APH-PRES)*T+ (PRES-APF)*tstar)/ (APH-APF);
   else
   {
      ZTMSL=Z1;
      tstar=ztsz;
      ZHTS=GEOS*zrg;
      if (ZHTS > 2000. && Z1 > 298.)
      {
         ZTMSL=298.;
         if (ZHTS < 2500.) ZTMSL=0.002*((2500.-ZHTS)*Z1+(ZHTS-2000.)*ZTMSL);
      }
      if ((ZTMSL-tstar) < 0.000001) ZALPH=0.;
      else if (GEOS > 0.0001 || GEOS < -0.0001)  ZALPH=RD*(ZTMSL-tstar)/GEOS;
      else                          ZALPH=RD*RLAPSE*zrg;
      ZALP=ZALPH*log(PRES/APH);
      PEVAL=tstar*(1.0+ZALP*(1.0+ZALP*(0.5+0.16666666667*ZALP)));
   }
   return PEVAL;
}

double ExtraZ(double pres, double aph, double apf, double sg, double t)
{
   double zrg    = 1.0 / Grav;
   double alpha,tstar,tmsl,ZALP,ZALPAL;

   alpha = RD * RLAPSE * zrg;
   tstar = (1.0 + alpha * (aph/apf - 1.0)) * t;
   if (tstar < 255.0) tstar = 0.5 * (255.0 + tstar);
   tmsl = tstar + RLAPSE * zrg * sg;
   if (tmsl > 290.5 && tstar > 290.5)
   {
      tstar = 0.5 * (290.5 + tstar);
      tmsl  = tstar;
   }
   if (tmsl > 290.5 && tstar <= 290.5) tmsl = 290.5;
   if (tmsl-tstar < 0.000001 && tstar-tmsl < 0.000001) alpha = 0.0;
   else if (sg > 0.0001 || sg < -0.0001) alpha = RD * (tmsl-tstar) / sg;
   ZALP   = log(pres/aph);
   ZALPAL = ZALP * alpha;
   return (sg - RD * tstar * ZALP * (1.0 + ZALPAL * (0.5 + ZALPAL/6.0))) * zrg;
}


void Interpolate_T(int Code)
{
   int lp, i;
   int nl,nh;
   double pres;

   int    *nx = vert_index;

   double *gt = &All[Code].hgp[0];
   double *pt = &All[Code].pgp[0];
   double *pf = &FullPress->hgp[0];
   double *ph = &HalfPress->hgp[0];
   
   for (lp = 0; lp < OutLevs; lp++)
   {
      pres = level[lp];
      for (i = 0; i < DimGP; i++)
      {
         nl = i  + DimGP * *nx;
         nh = nl + DimGP;
         if (nl < 0)              // Above top level
         {
            *pt = gt[i];
         }
         else if (nh >= Dim3GP)   // Below bottom level
         {
                 if (!mars && Code == TCODE) *pt = ExtraT(pres,ph[nh],pf[nl],Orography[i],gt[nl]);
            else                             *pt = gt[nl];
         }
         else                     // Inside
         {
            *pt = gt[nl] + (pres-pf[nl]) * (gt[nh]-gt[nl]) / (pf[nh]-pf[nl]);
         }
         ++nx;
         ++pt;
      }
   }
}

void Interpolate_Z(void)
{
   int lp, i;
   int nl,nh;
   double pres;

   int    *nx = vert_index;

   double *gz = &GeopotHeight->hgp[0];
   double *pz = &GeopotHeight->pgp[0];
   double *gt = &Temperature->hgp[0];
   double *pf = &FullPress->hgp[0];
   double *ph = &HalfPress->hgp[0];
   
   for (lp = 0; lp < OutLevs; lp++)
   {
      pres = level[lp];
      for (i = 0; i < DimGP; i++)
      {
         nl = i  + DimGP * *nx;
         if (pres > ph[nl+DimGP]) nl += DimGP;
         nh = nl + DimGP;
         if (nl < 0)              // Above top level
         {
            *pz = gz[i];
         }
         else if (nl >= Dim3GP)   // Below bottom level
         {
            if (mars) *pz = gz[nl-DimGP];
            else      *pz = ExtraZ(pres,ph[nl], pf[nl-DimGP],Orography[i],gt[nl-DimGP]);
         }
         else                     // Inside
         {
            *pz = gz[nl] + (pres-ph[nl]) * (gz[nh]-gz[nl]) / (ph[nh] - ph[nl]);
         }
         ++nx;
         ++pz;
      }
   }
}


void CheckDependencies(void)
{

           u_wind->needed = (  Divergence->needed &&   !Divergence->detected) ||
                            (   Vorticity->needed &&    !Vorticity->detected) ||
                            (     VeloPot->needed &&      !VeloPot->detected) ||
                            (     StreamF->needed &&      !StreamF->detected) ||
                            (       Omega->needed &&        !Omega->detected) ||
                            (       speed->needed &&        !speed->detected) ||
                            (      v_wind->needed &&       !v_wind->detected) ||
                                   u_wind->selected;
           v_wind->needed = (  Divergence->needed &&   !Divergence->detected) ||
                            (   Vorticity->needed &&    !Vorticity->detected) ||
                            (     VeloPot->needed &&      !VeloPot->detected) ||
                            (     StreamF->needed &&      !StreamF->detected) ||
                            (       Omega->needed &&        !Omega->detected) ||
                            (       speed->needed &&        !speed->detected) ||
                            (      u_wind->needed &&       !u_wind->detected) ||
                                   v_wind->selected;
       Divergence->needed = (      u_wind->needed &&       !u_wind->detected) ||
                            (      v_wind->needed &&       !v_wind->detected) ||
                            (       Omega->needed &&        !Omega->detected) ||
                            (     VeloPot->needed &&      !VeloPot->detected) ||
                               Divergence->selected;
        Vorticity->needed = (      u_wind->needed &&       !u_wind->detected) ||
                            (      v_wind->needed &&       !v_wind->detected) ||
                            (     StreamF->needed &&      !StreamF->detected) ||
                                Vorticity->selected;
         Humidity->needed = (GeopotHeight->needed && !GeopotHeight->detected) ||
                            (   Rhumidity->needed &&    !Rhumidity->detected) ||
                                 Humidity->selected;
              Ps->needed |=         dpsdx->needed ||         dpsdy->needed    ||
                                Rhumidity->needed ||         Omega->needed;
            LnPs->needed |=            Ps->needed;
      Temperature->needed = (GeopotHeight->needed && !GeopotHeight->detected) ||
                            (   Rhumidity->needed &&    !Rhumidity->detected) ||
                            (         SLP->needed &&          !SLP->detected) ||
                            (      w_wind->needed &&       !w_wind->detected) ||
                              Temperature->selected;

          All[176].needed = (    net_heat->needed &&     !net_heat->detected) ||
                            (     net_bot->needed &&      !net_bot->detected) ||
                            (     net_atm->needed &&      !net_atm->detected) ||
                            (      sw_atm->needed &&       !sw_atm->detected) ||
                                  All[176].selected;
          All[177].needed = (    net_heat->needed &&     !net_heat->detected) ||
                            (     net_bot->needed &&      !net_bot->detected) ||
                            (     net_atm->needed &&      !net_atm->detected) ||
                            (      lw_atm->needed &&       !lw_atm->detected) ||
                                  All[177].selected;
          All[178].needed = (     net_top->needed &&      !net_top->detected) ||
                            (     net_atm->needed &&      !net_atm->detected) ||
                            (      sw_atm->needed &&       !sw_atm->detected) ||
                                  All[178].selected;
          All[179].needed = (     net_top->needed &&      !net_top->detected) ||
                            (     net_atm->needed &&      !net_atm->detected) ||
                            (      lw_atm->needed &&       !lw_atm->detected) ||
                                  All[179].selected;
}

void CheckContent(void)
{
   int code;

   for (code = 0; code < 256; code++)
   {
      if (code == GEOSCODE) continue;
      if (code ==  SLPCODE) continue;
      if (code ==    ZCODE) continue;
      if (code ==  STRCODE) continue;
      if (code ==  VELCODE) continue;
      if (code ==    UCODE) continue;
      if (code ==    VCODE) continue;
      if (code ==    WCODE) continue;
      if (code ==   RHCODE) continue;
      if (code ==   PSCODE) continue;
      if (code ==   WZCODE) continue;
      if (code ==   SHCODE)
      {
         if (All[code].needed && !All[code].selected &&
             All[code].hsp.size() == 0 &&
             All[code].hgp.size() == 0 &&
             HumInfo == 0)
         {
            printf("\n ********* I N F O **********\n");
            printf(  " * No humidity in data file *\n");
            printf(  " *  Humidity set to zero !  *\n");
            printf(  " ****************************\n");
            All[code].needed = FALSE;
            HumInfo = 1;
         }
      }
      else
      {
         if (All[code].needed &&
             All[code].hsp.size() == 0 &&
             All[code].hgp.size() == 0)
         {
            printf("\n ****** E R R O R ******\n");
            printf(" * Code  %3d not found *\n",code);
            printf(" ***********************\n");
            exit(1);
         }
      }
   }
}

void Dependencies(void)
{
   int code;

   for (code = 0; code < CODES; code++)
   All[code].needed = All[code].selected;

   CheckDependencies();

   if (OutRep >= PRE_GRID)
   {
      u_wind->needed |= Divergence->needed;
      v_wind->needed |= Divergence->needed;
      u_wind->needed |=  Vorticity->needed;
      v_wind->needed |=  Vorticity->needed;
      u_wind->needed |=    VeloPot->needed;
      v_wind->needed |=    VeloPot->needed;
      u_wind->needed |=    StreamF->needed;
      v_wind->needed |=    StreamF->needed;
   }

   Omega->needed |= w_wind->needed;
   dpsdx->needed |=  Omega->needed;
   dpsdy->needed |=  Omega->needed;

   if (VerType == 'p')
   {
             LnPs->needed  = TRUE;
       Divergence->needed |=        Omega->needed;
           u_wind->needed |=        Omega->needed;
           v_wind->needed |=        Omega->needed;
           u_wind->needed |=        speed->needed;
           v_wind->needed |=        speed->needed;
      Temperature->needed |=        Omega->needed || SLP->needed;
         Humidity->needed |= GeopotHeight->needed || Rhumidity->needed;
      Temperature->needed |= GeopotHeight->needed || Rhumidity->needed ||
                                   ThetaF->needed;
   }

   Divergence->needed |=    u_wind->needed || v_wind->needed;
    Vorticity->needed |=    u_wind->needed || v_wind->needed;
   Divergence->needed |=   VeloPot->needed;
    Vorticity->needed |=   StreamF->needed;
         LnPs->needed |= HalfPress->needed || dpsdx->needed ||
                                Ps->needed || Rhumidity->needed;

 All[139].needed |=      ThetaF->selected;
 All[142].needed |=      precip->selected ||   net_water->selected ||
                    fresh_water->selected || surf_runoff->selected;
 All[143].needed |=      precip->selected ||   net_water->selected ||
                    fresh_water->selected || surf_runoff->selected;
 All[146].needed |=    net_heat->selected;     /* sensible heat */
 All[147].needed |=    net_heat->selected;     /* latent   heat */
 All[160].needed |=   net_water->selected;     /* Runoff        */
 All[176].needed |=    net_heat->selected ||
                        net_bot->selected ||
                        net_atm->selected ||      sw_atm->selected;
 All[177].needed |=    net_heat->selected ||
                        net_bot->selected ||
                        net_atm->selected ||      lw_atm->selected;
 All[178].needed |=     net_top->selected ||
                        net_atm->selected ||      sw_atm->selected;
 All[179].needed |=     net_top->selected ||
                        net_atm->selected ||      lw_atm->selected;
 All[182].needed |=   net_water->selected ||
                    fresh_water->selected || surf_runoff->selected;
 All[218].needed |=    net_heat->selected;     /* snow melt         */
 All[221].needed |= surf_runoff->selected;     /* snow depth change */

}

void Speed(double *speed, double *u, double *v)
{
   int i;

   for (i = 0; i < Dim3GP; i++)
   speed[i] = sqrt(u[i] * u[i] + v[i] * v[i]);
}


// ======================================================
// Compute derivation d(LnPs)/d(sin phi) in fourier space
// ======================================================

void Deriva(double field[], double derilam[])
{
   int l,n;
   int i;

   i = 0;
   for (n = 0; n < Waves    ; n++)
   {
     for (l = 0; l < Lats; l++,i++) // cosine coefficients
       derilam[i] = -n * field[i+Lats] * DerivationFactor[l];
     for (l = 0; l < Lats; l++,i++) //   sine coefficients
       derilam[i] =  n * field[i-Lats] * DerivationFactor[l];
   }
}

void scaluv(double *fu, double Factor[], int nlat, int nlot)
{
   for (int lot = 0; lot < nlot; lot++)
   for (int lat = 0; lat < nlat; lat++)
   {
      *fu++ *= Factor[lat];
   }
}

void uv2dv(double *fu, double *fv, double *sd, double *sv,
           double *pol2, double *poli, int klev, int nlat, int nt)
{
   int lev,jmm,jfc,lat;
   double dir,dii,vor,voi;
   double *ufr,*ufi,*vfr,*vfi;
   double *po2,*pod;
   double  zo2, zod;
   double  gmuq;
   for (lev = 0; lev < klev; lev++)
   {
      if (PolyDisk)
      {
         rewind(pol2f);
         rewind(polif);
      }
      po2 = pol2;
      pod = poli;
      for (jmm = 0; jmm <= nt; jmm++)
      {
         for (jfc = jmm; jfc <= nt; jfc++)
         {
            ufr = fu        ;
            ufi = fu + nlat ;
            vfr = fv        ;
            vfi = fv + nlat ;
            dir = 0.0       ;
            dii = 0.0       ;
            vor = 0.0       ;
            voi = 0.0       ;
            if (PolyDisk)
            {
               fread(po2=pol2,sizeof(double),Lats,pol2f);
               fread(pod=poli,sizeof(double),Lats,polif);
            }
            for (lat = 0; lat < nlat; lat++)
            {
               gmuq = 1.0 - Outlat->gmu[lat] * Outlat->gmu[lat];
               zod  = *pod * 0.5 * jmm * Outlat->gwt[lat] / (PlanetRadius * gmuq);
               zo2  = *po2 * Outlat->gwt[lat] / gmuq;
               dir += *vfr * zo2 - *ufi * zod;
               dii += *vfi * zo2 + *ufr * zod;
               vor -= *ufr * zo2 + *vfi * zod;
               voi -= *ufi * zo2 - *vfr * zod;
               ufr++;
               ufi++;
               vfr++;
               vfi++;
               po2++;
               pod++;
            }
            *sd++ = dir;
            *sd++ = dii;
            *sv++ = vor;
            *sv++ = voi;
         }
         fu += 2 * nlat;
         fv += 2 * nlat;
      }
   }
}


void genind(int *Interpolation_Index, double lv[],
            double *Full_Level_Pressure, int DimGP, int OutLevs)
{
   int h,k,l;
   int   *nx;
   double Pressure,*pf;
   nx = Interpolation_Index;
   for (h=0 ; h < DimGP * OutLevs ; ++h) nx[h] = -1;
   for (k = 0; k<OutLevs; k++)
   {
      Pressure = lv[k];
      pf       = Full_Level_Pressure;
      for (l = 0; l<SigLevs; l++)
      for (h = 0; h<DimGP ; h++)
      {
         if (Pressure > *pf) nx[h] = l;
         pf++;
      }
      nx += DimGP;
   }
}


void theta(double *PThetaF, double *PThetaH, double *PH, double *PS,
           double *TF, double *TS, int DimGP, int Dim3GP)
{
   int h,l;
   double  Kappa = RD / RCPD;
   double *ThetaH = PThetaH;
   double *ThetaF = PThetaF;
 
   for (h = 0; h < DimGP; h++) ThetaH[h] = 0.0;
   ThetaH += DimGP;
   for (l = 0; l < SigLevs - 1; l++)
   {
      for (h = 0; h < DimGP; h++)
      {
         ThetaH[h] = 0.5 * (TF[h] + TF[h+DimGP]) * pow((PS[h]/PH[h]),Kappa);
      }
      PH += DimGP;
      TF += DimGP;
      ThetaH += DimGP;
   }
   memcpy(ThetaH,TS,DimGP * sizeof(double));
   ThetaH = PThetaH;
   for (h = 0; h < Dim3GP; h++)
   {
      ThetaF[h] = 0.5 * (ThetaH[h] + ThetaH[h+DimGP]);
   }
}


void presh(double *pf, double *php, double *vct, double *ps)
{
   int h,l;
   double zp,ze;
   double *ph = php;

   for (l = 0; l<SigLevs; l++)
   {
      zp = vct[l];
      ze = vct[l+SigLevs+1];
      for (h = 0; h<DimGP; h++) ph[h] = zp + ze * ps[h];
      ph += DimGP;
   }
   memcpy(ph,ps,DimGP * sizeof(double));
   ph = php;
   for (h = 0; h<Dim3GP; h++) pf[h] = 0.5 * (ph[h] + ph[h+DimGP]);
}


/*****************************/
/* Compute relative humidity */
/*****************************/

double relhum(double q, double t, double p)
{
   double rh;
   double gascon;
   double rv;
   double TMELT;
   double ra1;
   double ra2;
   double ra4;
   double rdbrv;
   double zqsat;

   rv     = 461.51;
   TMELT  = 273.16;
   gascon = 287.0 ;
   ra1    = 610.78;
   ra2    =  17.2693882;
   ra4    =  35.86;
   rdbrv  = gascon / rv;

   zqsat  = rdbrv * ra1 * exp(ra2 * (t-TMELT) / (t-ra4)) / p;
   zqsat *= 1.0 / (1.0 - (1.0 / rdbrv - 1.0) * zqsat);
   rh     = q * 100.0 / zqsat;
   if (rh <   0.0) rh =   0.0;
   if (rh > 100.0) rh = 100.0;

   return rh;
}

/*****************************/
/* Compute relative humidity */
/*****************************/

void sh2rh(double *sphum, double *rhum, double *t, int lev)
{
   int jhor,jlev;
   double *pp;    // pressure

   pp = &FullPress->hgp[0];

   for (jlev = 0; jlev < lev; jlev++)
      for (jhor = 0; jhor < DimGP; jhor++)
         *rhum++ = relhum(*sphum++,*t++,*pp++);
}


void dv2ps(double *div, double *pot, int lev)
{
   for (int l = 0; l <  lev       ; l++)
   for (int m = 0; m <= Truncation; m++)
   for (int n = m; n <= Truncation; n++)
   {
      if (n)
      {
         *pot++ = *div++ * SQUARE_RADIUS / (n * n + n);
         *pot++ = *div++ * SQUARE_RADIUS / (n * n + n);
      }
      else
      {
         *pot++ = 0.0;
         *pot++ = 0.0;
         div   += 2;
      }
   }
}


void MakeGeopotHeight(double *geop, double* gt, double *gq, double *ph, int nhor, int nlev)
{
   int i;

   double VTMP = (RV / RD) - 1.0;
   double zrg  = 1.0 / Grav;

   if (gq) /* Humidity is present */ {
      for (i = nhor * nlev - 1; i >= nhor; i--)
         geop[i] = geop[i+nhor] + RD * gt[i] * (1.0 + VTMP * gq[i])
                 * log(ph[i+nhor] / ph[i]);

      for (i = 0; i < nhor; i++)
         geop[i] = geop[i+nhor] + RD * gt[i] * (1.0 + VTMP * gq[i])
                 * 2.0 * log(2.0);
   }
   else    /* No humidity */ {
      for (i = nhor * nlev - 1; i >= nhor; i--)
         geop[i] = geop[i+nhor] + RD * gt[i] * log(ph[i+nhor] / ph[i]);

      for (i = 0; i < nhor; i++)
         geop[i] = geop[i+nhor] + RD * gt[i] * 2.0 * log(2.0);
   }

   for (i = 0; i < nhor * (nlev+1); i++) geop[i] *= zrg;
}


void gp2fc_uv(void)
{
   u_wind->SetPFour();
   v_wind->SetPFour();
   gp2fc(&u_wind->pgp[0],&u_wind->pfc[0],Lats,Lons,OutLevs,Fouriers);
   gp2fc(&v_wind->pgp[0],&v_wind->pfc[0],Lats,Lons,OutLevs,Fouriers);
}


void fc2sp_uv(void)
{
   scaluv(&u_wind->pfc[0],CosPhi,Lats,Fouriers*OutLevs);
   scaluv(&v_wind->pfc[0],CosPhi,Lats,Fouriers*OutLevs);

   u_wind->SetPSpec();
   v_wind->SetPSpec();
   Divergence->SetPSpec();
   Vorticity->SetPSpec();

   fc2sp(&u_wind->pfc[0],&u_wind->psp[0],OutLevs,Lats,Truncation);
   fc2sp(&v_wind->pfc[0],&v_wind->psp[0],OutLevs,Lats,Truncation);

   uv2dv(&u_wind->pfc[0],&v_wind->pfc[0],
      &Divergence->psp[0],&Vorticity->psp[0],pol2,poli,OutLevs,Lats,Truncation);

   if (VeloPot->needed)
   {
      VeloPot->plev = OutLevs;
      VeloPot->SetPSpec();
      dv2ps(&Divergence->psp[0],&VeloPot->psp[0],OutLevs);
   }

   if (StreamF->needed)
   {
      StreamF->plev = OutLevs;
      StreamF->SetPSpec();
      dv2ps(&Vorticity->psp[0],&StreamF->psp[0],OutLevs);
   }
}


void sp2fc_uv(void)
{
   for (int si = 0 ; si < 4 ; ++si)
   {
      int code = SpecialCodes[si];
      if (All[code].selected && All[code].psp.size())
      {
         All[code].SetPFour();
         sp2fci(&All[code].psp[0],&All[code].pfc[0],OutLevs);
      }
   }
   if (u_wind->selected && u_wind->pfc.size())
      scaluv(&u_wind->pfc[0],RevCosPhi,Lats,Fouriers*OutLevs);
   if (v_wind->selected && v_wind->pfc.size())
      scaluv(&v_wind->pfc[0],RevCosPhi,Lats,Fouriers*OutLevs);
}


void fc2gp_uv(void)
{
   for (int si = 0 ; si < 4 ; ++si)
   {
      int code = SpecialCodes[si];
      if (All[code].selected && All[code].pfc.size())
      {
         All[code].SetPGrid();
         fc2gp(&All[code].pfc[0],&All[code].pgp[0],Lats,Lons,OutLevs,Fouriers);
      }
   }
}


void PumaProcess(void)
{
   int code;

   MeanCount++; // Count term inside month
   TermCount++; // Count term 

#ifdef NETCDF_OUTPUT
   if (TermCount == 1 && NetCDF) NetVarDefine(); // Define NetCDF variables
#endif
   if (MeanCount == 1) CheckContent();           // Everything OK ?
   if (TermCount > 60) Debug = 0;                // Limit debug output

   // Reset level offset for all variables

   for (code = 0; code < CODES; code++) All[code].loff = 0;

   // Derive velocity potential and streamfunction from divergence and vorticity

   if (VeloPot->needed && !VeloPot->detected && VerType == 's')
   {
      VeloPot->SetHSpec(SigLevs,OutLevs,FALSE);
      dv2ps(&Divergence->hsp[0],&VeloPot->hsp[0],SigLevs);
   }

   if (StreamF->needed && !StreamF->detected && VerType == 's')
   {
      StreamF->SetHSpec(SigLevs,OutLevs,FALSE);
      dv2ps(&Vorticity->hsp[0],&StreamF->hsp[0],SigLevs);
   }

   // -------------------------
   // Output of spectral fields
   // -------------------------

   if (OutRep == HYB_SPEC)
   {
      HybSpec->Write_hspec();
      return;
   }

   // =====================================================
   // Transformation from spectral domain to fourier domain
   // =====================================================

   // Derive wind components u*cos(phi) and v*cos(phi)

   if ((u_wind->needed   ||  v_wind->needed) &&
      (!u_wind->detected || !v_wind->detected))
   {
      u_wind->SetHFour(SigLevs,OutLevs,FALSE);
      v_wind->SetHFour(SigLevs,OutLevs,FALSE);
      spvfc(&Divergence->hsp[0],&Vorticity->hsp[0],
            &u_wind->hfc[0]    ,&v_wind->hfc[0]   ,
            Divergence->hlev ,Lats,Fouriers,Truncation);
   }

   // If divergence and vorticity were needed for u and v computation
   // only, they are now released

   Vorticity->needed  = Vorticity->selected;
   Divergence->needed = Divergence->selected || Omega->needed;

   // Perform inverse Legendre transformation for all needed variables

   for (code = 0 ; code < CODES; code++)
   if (All[code].needed && All[code].hsp.size())
   {
      All[code].SetHFour(All[code].hlev,All[code].plev,All[code].twod);
      sp2fci(&All[code].hsp[0],&All[code].hfc[0],All[code].hlev);
   }

   // Compute d(Lnps)/dx and d(Lnps)/dy if needed

   if (dpsdx->needed || dpsdy->needed)
   {
      dpsdx->SetHFour(1,1,TRUE);
      dpsdy->SetHFour(1,1,TRUE);
      Deriva(&LnPs->hfc[0],&dpsdx->hfc[0]);
      sp2fcd(&LnPs->hsp[0],&dpsdy->hfc[0],1,Lats,Fouriers,Truncation);
   }

   /* ------------------------ */
   /* Output of fourier fields */
   /* ------------------------ */

   if (OutRep == HYB_FOUR)
   {
      HybFour->Write_hfour();
      return;
   }

   /* --------------------- */
   /* Output of zonal means */
   /* --------------------- */

   if (OutRep == HYB_ZONM)
   {
      HybSect->Write_hfour();
      return;
   }

   /* ============================ */
   /* Transformation to gridpoints */
   /* ============================ */

   if (SaveMemory) HybSpec->Clear_hspec();

   for (code = 0; code < CODES; code++)
   if (All[code].needed && All[code].hfc.size())
   {
      All[code].SetHGrid(All[code].hlev,All[code].plev,All[code].twod);
      fc2gp(&All[code].hfc[0],&All[code].hgp[0],Lats,Lons,All[code].hlev,Fouriers);
   }

   if (SaveMemory) HybSpec->Clear_hfour();

   if (Humidity->hgp.size()) // Remove spurious negative humidity
   {
      for (int i=0 ; i < Dim3GP ; ++i)
         if (Humidity->hgp[i] < 0.0) Humidity->hgp[i] = 0.0;
   }


   if (LnPs->hgp.size())
   {
      Ps->SetHGrid(1,1,TRUE);
      Ps->hgp = exp(LnPs->hgp);
   }

   LnPs->needed = LnPs->selected;

   if (Orography.size() != DimGP)
   {
      Orography.resize(DimGP);
      if (Geopotential->hgp.size()) Orography = Geopotential->hgp;
      else
      {
         if (Geopotential->selected || VerType == 'p')
         {
           CenterText("Orography not found - using zero orography");
	   Orography = 0.0;
         }
      }
   }

   Geopotential->needed = Geopotential->selected;

   if (Geopotential->needed && !Geopotential->hgp.size())
   {
      Geopotential->SetHGrid(1,1,TRUE);
      Geopotential->hgp = Orography;
   }

   // This section is implemented for pressure level fields only

   if (VerType == 'p' || Omega->needed)
   {
      FullPress->SetHGrid(SigLevs  ,OutLevs,FALSE);
      HalfPress->SetHGrid(SigLevs+1,OutLevs,FALSE);
      presh(&FullPress->hgp[0],&HalfPress->hgp[0],vct,&Ps->hgp[0]);
   
      if (ThetaF->needed)
      {
         ThetaF->SetHGrid(SigLevs,OutLevs,FALSE);
         ThetaH->SetHGrid(SigLevs,OutLevs,FALSE);
         theta(&ThetaF->hgp[0], &ThetaH->hgp[0], &HalfPress->hgp[0], &Ps->hgp[0],
               &Temperature->hgp[0], &Ts->hgp[0], DimGP, Dim3GP);
      }
   
      if (GeopotHeight->needed)
      {
         GeopotHeight->SetHGrid(SigLevs+1,OutLevs,FALSE);
         memcpy(&GeopotHeight->hgp[Dim3GP],&Orography[0],DimGP * sizeof(double));
         MakeGeopotHeight(&GeopotHeight->hgp[0],&Temperature->hgp[0],
                          &Humidity->hgp[0],&HalfPress->hgp[0],DimGP,SigLevs);
         Humidity->needed = Humidity->selected;
      }
   
      if (dpsdx->needed) dpsdx->hgp *= Ps->hgp;
      if (dpsdy->needed) dpsdy->hgp *= Ps->hgp;
   
      if (Omega->needed)
      {
         Omega->SetHGrid(SigLevs+1,OutLevs,FALSE);
         OMEGA();
         dpsdx->needed = dpsdx->selected;
         dpsdy->needed = dpsdy->selected;
      }
   
      if (w_wind->needed)
      {
         w_wind->SetHGrid(SigLevs,OutLevs,FALSE);
         Omega_w(&w_wind->hgp[0],&Omega->hgp[0],&Temperature->hgp[0],&FullPress->hgp[0]);
      }
   
      if (Rhumidity->needed)
      {
         Rhumidity->SetHGrid(SigLevs,OutLevs,FALSE);
         sh2rh(&Humidity->hgp[0],&Rhumidity->hgp[0],
               &Temperature->hgp[0],SigLevs);
   
         Temperature->needed = Temperature->selected;
            Humidity->needed =    Humidity->selected;
      }
   
      if (SLP->needed)
      {
         SLP->SetHGrid(1,1,TRUE);
         Extrap(&SLP->hgp[0],&HalfPress->hgp[0] + Dim3GP,
                &FullPress->hgp[0] + Dim3GP - DimGP , &Orography[0],
                &Temperature->hgp[0] + Dim3GP - DimGP , DimGP);
         Temperature->needed = Temperature->selected || GeopotHeight->selected;
      }

   } // endif (VerType == 'p') 

   if (speed->needed)
   {
      speed->SetHGrid(SigLevs,OutLevs,FALSE);
      Speed(&speed->hgp[0],&u_wind->hgp[0],&v_wind->hgp[0]);
   }
   
   if (precip->needed)
   {
      precip->SetHGrid(1,1,TRUE);
      precip->hgp = All[142].hgp + All[143].hgp;
   }

   if (net_top->needed)
   {
      net_top->SetHGrid(1,1,TRUE);
      net_top->hgp = All[178].hgp + All[179].hgp;
   }

   if (net_bot->needed)
   {
      net_bot->SetHGrid(1,1,TRUE);
      net_bot->hgp = All[176].hgp + All[177].hgp;
   }

   if (net_heat->needed)
   {
      net_heat->SetHGrid(1,1,TRUE);
      net_heat->hgp = All[218].hgp * L_times_rhoH2O
      + All[176].hgp + All[177].hgp + All[146].hgp + All[147].hgp;
   }

   if (net_water->needed)
   {
      net_water->SetHGrid(1,1,TRUE);
      net_water->hgp = All[182].hgp - All[160].hgp + All[142].hgp + All[143].hgp;
   }

   if (sw_atm->needed)
   {
      sw_atm->SetHGrid(1,1,TRUE);
      sw_atm->hgp = All[178].hgp - All[176].hgp;
   }

   if (lw_atm->needed)
   {
      lw_atm->SetHGrid(1,1,TRUE);
      lw_atm->hgp = All[179].hgp - All[177].hgp;
   }

   if (net_atm->needed)
   {
      net_atm->SetHGrid(1,1,TRUE);
      net_atm->hgp = All[178].hgp + All[179].hgp - All[176].hgp - All[177].hgp;
   }

   if (surf_runoff->needed)
   {
      surf_runoff->SetHGrid(1,1,TRUE);
      surf_runoff->hgp = All[182].hgp - All[221].hgp + All[142].hgp + All[143].hgp;
   }

   if (fresh_water->needed)
   {
      fresh_water->SetHGrid(1,1,TRUE);
      fresh_water->hgp = All[142].hgp + All[143].hgp + All[182].hgp;
   }

   // =============================
   // Monthly means on hybrid grids
   // =============================

   if (Mean && OutRep == HYB_GRID)
   for (code = 0; code < CODES; code++)
   if (All[code].selected && All[code].hgp.size())
   {
      if (MeanCount == 1)
      {
         All[code].mgp.resize(All[code].hgp.size());
         All[code].mgp = All[code].hgp;
         All[code].hgp.resize(0);
      }
      else All[code].mgp += All[code].hgp;
      if (EndOfMonth)
      {
         double rmc = 1.0 / MeanCount;
         All[code].hgp = All[code].mgp * rmc;
         All[code].mgp.resize(0);
      }
   }

   // ----------------------------
   // Output of hybrid level grids
   // ----------------------------

   if (OutRep == HYB_GRID)
   {
      if (!Mean || EndOfMonth) HybGrid->Write_hgrid();
      if (SaveMemory) HybGrid->Clear_hgrid();
      return;
   }

   // ======================================
   // Vertical interpolation / extrapolation
   // ======================================

   if (vert_index == NULL) vert_index = new int[OutLevs*DimGP];
   genind(vert_index,level,&FullPress->hgp[0],DimGP,OutLevs);

   for (code = 0; code < CODES; code++)
   if (All[code].needed && All[code].hgp.size())
   {
      All[code].SetPGrid();
      if (!All[code].twod)
      {
         if (code == ZCODE) Interpolate_Z();
         else               Interpolate_T(code);
      }
   }
   Temperature->needed = Temperature->selected;
   if (SaveMemory) HybGrid->Clear_hgrid();

   // ===========================
   // Computation of Montly Means
   // ===========================

   if (Mean)
   for (code = 0; code < CODES; code++)
   if (All[code].needed && All[code].pgp.size())
   {
      if (MeanCount == 1)
      {
         All[code].mgp.resize(All[code].pgp.size());
         All[code].mgp = All[code].pgp;
         All[code].pgp.resize(0);
      }
      else All[code].mgp += All[code].pgp;
      if (EndOfMonth)
      {
         double rmc = 1.0 / MeanCount;
         All[code].pgp = All[code].mgp * rmc;
         All[code].mgp.resize(0);
      }
   }

   if (Mean && !EndOfMonth)
   {
      if (SaveMemory) HybGrid->Clear_pgrid();
      return;
   }

   // ------------------------------
   // Output of pressure level grids
   // ------------------------------

   if (OutRep == PRE_GRID)
   {
      if (SpecialUV)
      {
         gp2fc_uv();
         fc2sp_uv();
         sp2fc_uv();
         fc2gp_uv();
      }

      HybGrid->Write_pgrid();
      if (SaveMemory) HybGrid->Clear_pgrid();
      return;
   }

   // ===============================
   // Transformation to fourier space
   // ===============================

   for (code = 0; code < CODES; code++)
   if (All[code].needed && All[code].pgp.size())
   {
      if (!All[code].pfc.size()) All[code].SetPFour();
      gp2fc(&All[code].pgp[0],&All[code].pfc[0],
            Lats,Lons,All[code].plev,Fouriers);
   }

   if (SaveMemory) HybGrid->Clear_pgrid();

   // ---------------------------------------
   // Output of fourier fields or zonal means
   // ---------------------------------------

   if (OutRep == PRE_FOUR || OutRep == PRE_ZONM)
   {
      if (SpecialUV)
      {
         fc2sp_uv();
         sp2fc_uv();
      }
      if (OutRep == PRE_FOUR) HybFour->Write_pfour();
      else                    HybSect->Write_pfour();
      if (SaveMemory) HybFour->Clear_pfour();
      return;
   }

   // ================================
   // Transformation to spectral space
   // ================================

   if (!SpecialUV && u_wind->pfc.size() && v_wind->pfc.size())
   {
      scaluv(&u_wind->pfc[0],CosPhi,Lats,Fouriers*OutLevs);
      scaluv(&v_wind->pfc[0],CosPhi,Lats,Fouriers*OutLevs);
   }

   for (code = 0; code < CODES; code++)
   if (All[code].needed && All[code].pfc.size() && !All[code].psp.size())
   {
      All[code].SetPSpec();
      fc2sp(&All[code].pfc[0],&All[code].psp[0],All[code].plev,Lats,Truncation);
   }

   if (SpecialUV) fc2sp_uv();

   if (SaveMemory) HybFour->Clear_pfour();

   // -------------------------
   // Output of spectral fields
   // -------------------------

   if (OutRep == PRE_SPEC)
   {
      HybSpec->Write_pspec();
      if (SaveMemory) HybSpec->Clear_pspec();
      return;
   }
}


void PostProcess(void)
{
   int l;
   char tb[COLS+2];
   if (EndOfMonth)
   {
      sprintf(tb,"Processed Month %2d   Year %4d", OldMonth,OldDate.tm_year);
      l = strlen(tb);
      if (MeanCount > 1)
      {
         if (Mean) sprintf(tb+l,"  (Mean of %3d Terms)",MeanCount);
         else      sprintf(tb+l,"   Terms %3d",MeanCount);
      }
      LeftText(tb);
      EndOfMonth = FALSE;
      MeanCount  =     0;
      MonthCount++ ;
   }
}

/* ================= */
/* switch input file */
/* ================= */

void SwitchFile(void)
{
   int l,YY,MM;

   fclose(fpi);
   l = strlen(ifile);
   if (l > 4 && ifile[l-4] == '.' && atoi(ifile+l-3) != 0) // .YYY
   {
      YY = atoi(ifile+l-3) + 1;
      sprintf(ifile+l-3,"%03d",YY);
   }
   else if (l > 5 && ifile[l-5] == '_' && atoi(ifile+l-4) != 0) // _YYYY
   {
      YY = atoi(ifile+l-4) + 1;
      sprintf(ifile+l-4,"%04d",YY);
   }
   else if (l > 7 && ifile[l-7] == '_' && atoi(ifile+l-6) != 0) // _YYYYMM
   {
      MM = atoi(ifile+l-2);
      YY = atoi(ifile+l-6);
      if (MM == 12) YY += 88;
      if (YY > 1) sprintf(ifile+l-6,"%06d",++YY);
   }
   else if (l > 5 && atoi(ifile+l-5) != 0) // "YYYMM" at end
   {
      MM = atoi(ifile+l-2);
      YY = atoi(ifile+l-5);
      if (MM == 12) YY += 88;
      if (YY > 1) sprintf(ifile+l-5,"%05d",++YY);
   }
   else if (l > 4 && atoi(ifile+l-4) != 0) // "YYMM" at end
   {
      MM = atoi(ifile+l-2);
      YY = atoi(ifile+l-4);
      if (MM == 12) YY += 88;
      if (YY > 1) sprintf(ifile+l-4,"%04d",++YY);
   }

   Multi--;
   printf("Continuation File: %s\n",ifile);
   fpi  = fopen(ifile,"rb");
}

/* ====================================== */
/* Interpolate gauss grid to regular grid */
/* ====================================== */

void Regauss(double *r, double *g, double *Ghi)
{
   int j,jlon,jlat;
   double np,sp,fn,fs,fe,fw;
   double rphi,rlam,gdx,rdx;
   double *Gam;

   Gam = new double[Gons];

   gdx = 360.0 / Gons;
   rdx = 360.0 / Lons;

   np = sp = 0.0;
   for (jlon = 0 ; jlon < Gons ; ++jlon)
   {
       np += g[jlon];
       sp += g[jlon + DimGG - Gons];
   }

   np /= Gons;
   sp /= Gons;

   for (jlat = 0 ; jlat < Lats ; ++jlat)
   {
      rphi = Outlat->Phi[jlat];
      if (rphi > Ghi[0]) // far north
      {
         fn = (rphi - Ghi[0]) / (90.0 - Ghi[0]);
         fs = 1.0 - fn;
         for (jlon = 0 ; jlon < Gons ; ++jlon)
            Gam[jlon] = fn * np + fs * g[jlon];
      }
      else if (rphi < Ghi[Gats-1]) // far south
      {
         fs = (Ghi[Gats-1] - rphi) / (Ghi[Gats-1] + 90.0);
         fn = 1.0 - fs;
         for (jlon = 0 ; jlon < Gons ; ++jlon)
            Gam[jlon] = fn * g[jlon + DimGG - Gons] + fs * sp;
      }
      else // inside
      {
         j = 0; // search neighboured gauss latitudes
         while (j < Gats-1 && rphi < Ghi[j]) ++j;
         fn = (rphi - Ghi[j]) / (Ghi[j-1] - Ghi[j]);
         fs = 1.0 - fn;
         for (jlon = 0 ; jlon < Gons ; ++jlon)
            Gam[jlon] = fn * g[jlon+(j-1)*Gons] + fs * g[jlon+j*Gons];
      }

      for (jlon = 0 ; jlon < Lons ; ++jlon)
      {
         rlam = jlon * rdx;
         j = (int)floor(rlam / gdx);
         fe = (rlam - j * gdx) / gdx;
         fw = 1.0 - fe;
         if (j >= Gons-1) r[jlon + jlat * Lons] = fw * Gam[j] + fe * Gam[0];
         else           r[jlon + jlat * Lons] = fw * Gam[j] + fe * Gam[j+1];
      }
   }
   delete Gam;
}


/*******************/
/* SetOutputHeader */
/*******************/

void SetOutputHeader(void)
{
   int MM,DD,HH;

   if (DPM > 99) // Workaround for months with more than 99 days 
   {
      MM = (TermCount / DPM) % 12 + 1;
      HeadOu[2] = OldDate.tm_year * 10000 + MM * 100;
      HeadOu[3] =  0;
      if (!Mean)  // Add day & hour info
      {
         HeadOu[2] += (TermCount % DPM) / DayDivisor + 1;
         HeadOu[3]  = 100 * (24 / DayDivisor) * ((TermCount % DPM) % DayDivisor);
      }
   }
   else
   {
      HeadOu[2] = OldDate.tm_year * 10000 + OldDate.tm_mon * 100;
      HeadOu[3] = 0;
      if (!Mean)
      {
         HeadOu[2] += OldDate.tm_mday;
         HeadOu[3]  = OldDate.tm_hour * 100 + OldDate.tm_min;
      }
   }
}


/*****************/
/* Puma Control */
/*****************/

void PumaControl(void)
{
   int i;
   int LevelOffset;
   int Eof;
   char tb[COLS+2];
   struct tm D1;
   struct tm D2;

   while (1)
   {
      Eof = ReadHeaderRecord();
      if (Eof && Multi)
      {
         SwitchFile();
         if (fpi) Eof = ReadHeaderRecord();
      }

      if (DataStep < 0.01) // Compute time interval
      {
         if (HeadIn[2] != HeadSt[2] || HeadIn[3] != HeadSt[3])
	     {
    	    HeadToDate(HeadSt,&D1);
            HeadToDate(HeadIn,&D2);
            DeltaDy  = D2.tm_mday - D1.tm_mday;
            DeltaHr  = D2.tm_hour - D1.tm_hour;
            DeltaMn  = D2.tm_min  - D1.tm_min ;
            if (DeltaDy < 0) DeltaDy = 1; // month changed after 1.st term
    	    DataStep = DeltaDy + DeltaHr / 24.0 + DeltaMn / 1440.0;
         }
      }

      if ((HeadIn[2] / 100 % 100) > LastMonth && DayDivisor == 0) // Ignore rest of file
      {
         Eof = 1;
         if (Multi)
         {
            SwitchFile();
            if (fpi) Eof = ReadHeaderRecord();
         }
      }

      if (Eof) // Process last read term and finish
      {
         EndOfMonth = TRUE;
         SetOutputHeader();
         PumaProcess();
         PostProcess();
         Dependencies();
         return;
      }

      DecodePumaHeader();

      if (NewMonth < FirstMonth) /* Skip months before FirstMonth */
      {
         SkipPumaArray();
         if (Debug)
         {
            if (RepGrib == REP_SPECTRAL) sprintf(tb,"T%04d",Truncation);
            else                         sprintf(tb,"N%04d",Lons);
            sprintf(tb+5," SKIPPED Code %3d   Level%6d   %2.2d.%2.2d.%2.2d   %2.2d:%2.2d",
            PumaCode,PumaLevel,NewDate.tm_mday,NewDate.tm_mon,NewDate.tm_year,
            NewDate.tm_hour,NewDate.tm_min);
            LeftText(tb);
         }
         continue;
      }

      if (Debug)
      {
         if (RepGrib == REP_SPECTRAL) sprintf(tb,"T%04d",Truncation);
         else                         sprintf(tb,"N%04d",Lons);
         sprintf(tb+5," Code %3d   Level%6d   %2.2d.%2.2d.%2.2d   %2.2d:%2.2d",
         PumaCode,PumaLevel,NewDate.tm_mday,NewDate.tm_mon,NewDate.tm_year,
         NewDate.tm_hour,NewDate.tm_min);
         LeftText(tb);
      }

      if (OldMonth > 0)
      {
         EndOfMonth = NewMonth != OldMonth;
         EndOfTerm = memcmp(&NewDate,&OldDate,sizeof(struct tm));
         if (EndOfTerm && MeanCount == DPM-1) EndOfMonth = 1;
         if (EndOfTerm)
         {
            SetOutputHeader();
            PumaProcess();
            PostProcess();
            Dependencies();
         }
      }
      OldDate  = NewDate;
      OldMonth = NewMonth;

      if (All[PumaCode].needed)
      {
         if (RepGrib == REP_SPECTRAL) // Spectral array
         { 
            if (PumaLevel) All[PumaCode].SetHSpec(SigLevs,OutLevs,FALSE);
            else           All[PumaCode].SetHSpec(      1,      1,TRUE );
            if (VerType != 's' || mom[PumaLevel])
            {
               ReadPumaArray(&All[PumaCode].hsp[0]+All[PumaCode].loff);
               All[PumaCode].loff += DimSP;
            }
            else SkipPumaArray();
         }
         else  // Gridpoint array
         {
            if (PumaLevel) All[PumaCode].SetHGrid(SigLevs,OutLevs,FALSE);
            else           All[PumaCode].SetHGrid(      1,      1,TRUE );
            if (VerType != 's' || mom[PumaLevel])
            {
               if (DimGP == DimGG)
                  ReadPumaArray(&All[PumaCode].hgp[0]+All[PumaCode].loff);
               else
               {
                  double *ArrayBuffer;
                  ArrayBuffer = new double[DimGG];
                  ReadPumaArray(ArrayBuffer);
                  Regauss(&All[PumaCode].hgp[0]+All[PumaCode].loff,ArrayBuffer,Gaulat->Phi);
                  delete ArrayBuffer;
               }
               All[PumaCode].loff += DimGP;
            }
            else SkipPumaArray();
         }
      }
      else SkipPumaArray();
   }
}

char *amatch(char *msr, const char *sub)
{
   int i,nm,ns;
   nm = strlen(msr);
   ns = strlen(sub);

   for (i = 0; i < nm-ns; i++)
   if (strncmp(msr+i,sub,ns) == 0) return (msr+i+ns);

   return NULL;
}

int scanpar(const char *name, int def)
{
   char *cp;
   int value;
   char tb[COLS+2];

   cp = amatch(namelist,name);
   if (cp == NULL)
   {
       value = def;
       sprintf(tb,"%10.10s = %8d  (default)",name,value);
   }
   else
   {
       value = atoi(cp);
       sprintf(tb,"%10.10s = %8d           ",name,value);
   }
   LeftText(tb);
   return value;
}
   
double scanreal(const char *name, double def)
{
   char *cp;
   double value;
   char tb[COLS+2];

   cp = amatch(namelist,name);
   if (cp == NULL)
   {
       value = def;
       sprintf(tb,"%10.10s = %8.3f  (default)",name,value);
   }
   else
   {
       value = strtod(cp,NULL);
       sprintf(tb,"%10.10s = %8.3f           ",name,value);
   }
   LeftText(tb);
   return value;
}
   
char scantex(const char *name, const char choice[])
{
   char *cp;
   char value;
   int i;
   char tb[COLS+2];

   value = choice[0];
   cp = amatch(namelist,name);
   if (cp)
   {
      while (isspace(*cp)) ++cp;
      for (i=1 ; i < strlen(choice) ; ++i)
      {
         if (*cp == choice[i]) value = *cp;
      }
   }
   sprintf(tb,"%10.10s = %c                  ",name,value);
   LeftText(tb);
   return value;
}
   
void scantime(void)
{
   char *cp,*icp;
   int time,i;
   char tb[512];

   nrqh = 0;

   cp = amatch(namelist,"timesel");
   if (cp == NULL)
   {
      hours[nrqh++] = -1;
      sprintf(tb,"%10.10s = all                ","timesel");
      LeftText(tb);
      return;
   }

   time = strtol(cp,&icp,10);
   while ((char *)icp != (char *)cp && nrqh < MAX_HOURS) {
      hours[nrqh++] = time;
      cp = icp;
      time = strtol(cp,&icp,10);
   }
   sprintf(tb,"%10.10s = ","timesel");
   for (time = 0 ; time < nrqh ; ++time)
   {
       i = strlen(tb);
       sprintf(tb+i," %02d",hours[time]);
   }
   LeftText(tb);
}

void PrintCodes(void)
{
   int code;
   char tb[COLS+2];

   DashLine();
   CenterText("Code Id       Name                             Units            ");
   DashLine();
   for (code=0 ; code < MAXCODES ; ++code) if (strncmp(All[code].Na,"Code",4))
   {
      sprintf(tb,"%4d %-8.8s %-32.32s %-16.16s",code,All[code].Id,All[code].Na,All[code].Un);
      CenterText(tb);
   }
}

int CodeOrName(char *a, char **b)
{
   int l,code;
   while (*a == ' ') ++a;
   if (*a == '+' || *a == '-' || (*a >= '0' && *a <='9'))
      return strtol(a,b,10);
   for (code = 0 ; code < CODES ; ++code)
   {
      l = strlen(All[code].Id);
      if (!strncmp(All[code].Id,a,l) && *(a+l) == ' ')
      {
        *b = a+l;
        return code;
      }
   }
   return 0;
}

      
void scancode(void)
{
   char *cp,*icp;
   int code;
   char tb[COLS+2];

   cp = amatch(namelist,"code");
   if (cp == NULL) Abort(" *** No code selected for output ***");

   code = CodeOrName(cp,&icp);
   DashLine();
   while (code > 0 && code < CODES)
   {
      sprintf(tb,"Code %5d = %-6.6s %-24.24s",code,All[code].Id,All[code].Na);
      LeftText(tb);
      All[code].selected = 1;
      cp = icp;
      code = CodeOrName(cp,&icp);
   }
}

void scanmol(void)
{
   char *cp,*icp;
   int lev;

   nrml = 0;
   cp = amatch(namelist,"modlev");
   if (cp == NULL) return;
   lev = strtol(cp,&icp,10);
   while (lev > 0 && nrml < MAX_LEVELS)
   {
      mol[nrml++] = lev;
      mom[lev] = 1;
      cp = icp;
      lev = strtol(cp,&icp,10);
   }
}

void scanhPa(void)
{
   char *cp,*icp;
   double lev;

   nrpl = 0;
   cp = amatch(namelist,"hpa");
   if (!cp) return;
   lev = strtod(cp,&icp);
   while (lev > 0 && nrpl < MAX_LEVELS)
   {
      hPa[nrpl++] = lev;
      cp = icp;
      lev = strtod(cp,&icp);
   }
}

void scanattributes(void)
{
   char *cp;
   int i;

   nattr = 0;
   cp = amatch(namelist,"attributes");
   if (!cp) return;

   while (nattr < ATTR_MAX)
   {
      i = 0;
      while (*cp == ':' || isspace(*cp)) ++cp;
      while ((isalnum(*cp) || *cp == '_') && i < 80) AttrNam[nattr][i++] = *cp++;
      while (isspace(*cp)) ++cp;
      if (*cp != '=') break;
      else ++cp;
      while (isspace(*cp)) ++cp;
      if (*cp != '"') break;
      else ++cp;
      i = 0;
      while (*cp != '"'  && i < 80) AttrVal[nattr][i++] = *cp++;
      ++cp;
      while (isspace(*cp)) ++cp;
      if (*cp != ';') break;
      ++cp;
      ++nattr;
   }
}

void InitAll(void)
{
   char Id[MAX_ID_LEN];
   char Na[MAX_NA_LEN];
 
   for (int code = 0 ; code < MAXCODES ; ++code)
   {
      sprintf(Id,"var%3.3d",code);
      sprintf(Na,"Code[%d]",code);
      All[code].Init(Id,Na,"1",0);
      All[code].code = code;
   }

   All[110].Init("mld"  ,"mixed_layer_depth"               ,"m"        ,1); // Not standard
   All[129].Init("sg"   ,"surface_geopotential"            ,"m2 s-2"   ,1);
   All[130].Init("ta"   ,"air_temperature"                 ,"K"        ,0);
   All[131].Init("ua"   ,"eastward_wind"                   ,"m s-1"    ,0);
   All[132].Init("va"   ,"northward_wind"                  ,"m s-1"    ,0);
   All[133].Init("hus"  ,"specific_humidity"               ,"1"        ,0);
   All[134].Init("ps"   ,"surface_air_pressure"            ,"hPa"      ,1);
   All[135].Init("wap"  ,"vertical_air_velocity"           ,"Pa s-1"   ,0); // shortened
   All[137].Init("wa"   ,"upward_wind"                     ,"m s-1"    ,0); // Not standard
   All[138].Init("zeta" ,"atm_relative_vorticity"          ,"s-1"      ,0);
   All[139].Init("ts"   ,"surface_temperature"             ,"K"        ,1);
   All[140].Init("mrso" ,"lwe_of_soil_moisture_content"    ,"m"        ,1); // shortened
   All[141].Init("snd"  ,"surface_snow_thickness"          ,"m"        ,1);
   All[142].Init("prl"  ,"lwe_of_large_scale_precipitation","m s-1"    ,1); // rate !!
   All[143].Init("prc"  ,"convective_precipitation_rate"   ,"m s-1"    ,1);
   All[144].Init("prsn" ,"lwe_of_snowfall_amount"          ,"m s-1"    ,1); // rate !!
   All[145].Init("bld"  ,"dissipation_in_atmosphere_bl"    ,"W m-2"    ,1); // shortened
   All[146].Init("hfss" ,"surface_sensible_heat_flux"      ,"W m-2"    ,1); // shortened
   All[147].Init("hfls" ,"surface_latent_heat_flux"        ,"W m-2"    ,1); // shortened
   All[148].Init("stf"  ,"streamfunction"                  ,"m2 s-2"   ,0); // Not standard
   All[149].Init("psi"  ,"velocity_potential"              ,"m2 s-2"   ,0); // Not standard
   All[151].Init("psl"  ,"air_pressure_at_sea_level"       ,"hPa"      ,1);
   All[152].Init("pl"   ,"log_surface_pressure"            ,"1"        ,1);
   All[155].Init("d"    ,"divergence_of_wind"              ,"s-1"      ,0);
   All[156].Init("zg"   ,"geopotential_height"             ,"m"        ,0);
   All[157].Init("hur"  ,"relative_humidity"               ,"1"        ,0);
   All[158].Init("tps"  ,"tendency_of_surface_air_pressure","Pa s-1"   ,1);
   All[159].Init("u3"   ,"ustar"                           ,"m3 s-3"   ,1); // Not standard
   All[160].Init("mrro" ,"surface_runoff"                  ,"m s-1"    ,1); // Not standard
   All[161].Init("clw"  ,"liquid_water_content"            ,"1"        ,0); // Not standard
   All[162].Init("cl"   ,"cloud_area_fraction_in_layer"    ,"1"        ,0); // Not standard
   All[163].Init("tcc"  ,"total_cloud_cover"               ,"1"        ,1); // Not standard
   All[164].Init("clt"  ,"cloud_area_fraction"             ,"1"        ,1);
   All[165].Init("uas"  ,"eastward_wind_10m"               ,"m s-1"    ,1); // shortened
   All[166].Init("vas"  ,"northward_wind_10m"              ,"m s-1"    ,1); // shortened
   All[167].Init("tas"  ,"air_temperature_2m"              ,"K"        ,1); // shortened
   All[168].Init("td2m" ,"dew_point_temperature_2m"        ,"K"        ,1); // shortened
   All[169].Init("tsa"  ,"surface_temperature_accumulated" ,"K"        ,1); // Not standard
   All[170].Init("tsod" ,"deep_soil_temperature"           ,"K"        ,1);
   All[171].Init("dsw"  ,"deep_soil_wetness"               ,"1"        ,1);
   All[172].Init("lsm"  ,"land_binary_mask"                ,"1"        ,1);
   All[173].Init("z0"   ,"surface_roughness_length"        ,"m"        ,1);
   All[174].Init("alb"  ,"surface_albedo"                  ,"1"        ,1); // Not standard
   All[175].Init("as"   ,"surface_albedo"                  ,"1"        ,1); // Not standard
   All[176].Init("rss"  ,"surface_net_shortwave_flux"      ,"W m-2"    ,1); // shortened
   All[177].Init("rls"  ,"surface_net_longwave_flux"       ,"W m-2"    ,1); // shortened
   All[178].Init("rst"  ,"toa_net_shortwave_flux"          ,"W m-2"    ,1); // shortened
   All[179].Init("rlut" ,"toa_net_longwave_flux"           ,"W m-2"    ,1); // shortened
   All[180].Init("tauu" ,"surface_eastward_stress"         ,"Pa"       ,1); // shortened
   All[181].Init("tauv" ,"surface_northward_stress"        ,"Pa"       ,1); // shortened
   All[182].Init("evap" ,"lwe_of_water_evaporation"        ,"m s-1"    ,1); // rate !!
   All[183].Init("tso"  ,"climate_deep_soil_temperature"   ,"K"        ,1); // Not standard
   All[184].Init("wsoi" ,"climate_deep_soil_wetness"       ,"1"        ,1);
   All[199].Init("vegc" ,"vegetation_cover"                ,"1"        ,1); // Not standard
   All[203].Init("rsut" ,"toa_outgoing_shortwave_flux"     ,"W m-2"    ,1); // Not standard
   All[204].Init("ssru" ,"surface_solar_radiation_upward"  ,"W m-2"    ,1); // Not standard
   All[205].Init("stru" ,"surface_thermal_radiation_upward","W m-2"    ,1); // Not standard
   All[207].Init("tso2" ,"soil_temperature_level_2"        ,"K"        ,1); // Not standard
   All[208].Init("tso3" ,"soil_temperature_level_3"        ,"K"        ,1); // Not standard
   All[209].Init("tso4" ,"soil_temperature_level_4"        ,"K"        ,1); // Not standard
   All[210].Init("sic"  ,"sea_ice_cover"                   ,"1"        ,1); // Not standard
   All[211].Init("sit"  ,"sea_ice_thickness"               ,"m"        ,1); // Not standard
   All[212].Init("vegf" ,"forest_cover"                    ,"1"        ,1); // Not standard
   All[218].Init("snm"  ,"snow_melt"                       ,"m s-1"    ,1); // Not standard
   All[221].Init("sndc" ,"snow_depth_change"               ,"m s-1"    ,1); // Not standard
   All[230].Init("prw"  ,"atmosphere_water_vapor_content"  ,"kg m-2"   ,1); // Not standard
   All[232].Init("glac" ,"glacier_cover"                   ,"1"        ,1); // Not standard
   All[238].Init("tsn"  ,"snow_temperature"                ,"K"        ,1);
   All[259].Init("spd"  ,"wind_speed"                      ,"m s-1"    ,0); // Not standard
   All[260].Init("pr"   ,"total_precipitation"             ,"m s-1"    ,1); // Not standard
   All[261].Init("ntr"  ,"net_top_radiation"               ,"W m-2"    ,1); // Not standard
   All[262].Init("nbr"  ,"net_bottom_radiation"            ,"W m-2"    ,1); // Not standard
   All[263].Init("hfns" ,"surface_downward_heat_flux"      ,"W m-2"    ,1); // shortened
   All[264].Init("wfn"  ,"net_water_flux"                  ,"m s-1"    ,1); // Not standard
   All[273].Init("dpdx" ,"d(ps)/dx"                        ,"Pa m-1"   ,1); // Not standard
   All[274].Init("dpdy" ,"d(ps)/dy"                        ,"Pa m-1"   ,1); // Not standard
   All[277].Init("hlpr" ,"half_level_pressure"             ,"Pa"       ,0); // Not standard
   All[278].Init("flpr" ,"full_level_pressure"             ,"Pa"       ,0); // Not standard
}

void Usage(void)
{
   char  Line[132];
   fpi = fopen("/usr/local/doc/burn7.txt","r");
   if (fpi)
   do
   {
      fgets(Line,130,fpi);
      printf("%s",Line);
   }
   while (!feof(fpi) && Line[0] != '#');
   if (fpi) fclose(fpi);

   printf("\nburn7 [options] InputFile OutputFile <namelist >printout\n");
   printf("     option -h : help (this output)\n");
   printf("     option -c : print available codes and names\n");
   printf("     option -d : debug mode (verbose output)\n");
   printf("     option -m : Mean=1 output (override namelist option)\n");
   printf("     option -n : NetCDF output (override namelist option)\n");
   printf("     option -s : Save memory (increases CPU time)\n");
   printf("    InputFile  : PUMA or PlaSim data file\n");
#ifdef NETCDF_OUTPUT
   printf("  <OutputFile> : SERVICE, or NetCDF format file\n");
#else
   printf("  <OutputFile> : SERVICE format file\n");
#endif
   printf("  namelist is read from  <stdin>\n");
   printf("  printout is written to <stdout>\n\n");
   exit(1);
}

void parini(void)
{
   char c;
   unsigned int i;
   int jind;
   int lc;
   char tb[COLS+2];

   i = 1;
   namelist[0] = ' ';
   lc = 1;
   c = getchar();
   while (!feof(stdin) && i < MAX_NL)
   {
      if (c == ':' ) lc = 0; // No conversion to lower case 
      if (c == '\n') lc = 1;
      if (c == '#') // Skip comment
      {
         while (!feof(stdin) && c != '\n') c = getchar();
      }
      if (lc)
      {
              if ((c >= '0' && c <= '9') ||
                  c == '-' || c == '.')  namelist[i++] = c;
         else if (c >= 'a' && c <= 'z')  namelist[i++] = c;
         else if (c >= 'A' && c <= 'Z')  namelist[i++] = tolower(c);
         else c = ' ';
         if (c == ' ' && namelist[i-1] != ' ') namelist[i++] = c;
      }
      else namelist[i++] = c;
      c = getchar();
   }
   namelist[i++] = ' ';
   namelist[i] = 0;

   if (Debug)
   {
      sprintf(tb,"Length of namelist: %d bytes",(int)strlen(namelist));
      LeftText(tb);
      for (i = 0; i<strlen(namelist); i+=40)
      {
         sprintf(tb,"namelist[%02d]=%-40.40s",i,namelist+i);
         LeftText(tb);
      }
      StarLine();
   }

   Lats       = scanpar("lats",Lats);
   Lons       = scanpar("lons",Lons);
   DPM        = scanpar("dpm",0);
   if (DPM > 99 && DPM < 2400) /* Days Per Month */
   {
      DayDivisor = DPM / 100 + 1;
      while (24 % DayDivisor && DayDivisor < 24) ++DayDivisor;
      sprintf(tb,"%10.10s = %8d           ","daydivisor",DayDivisor);
      LeftText(tb);
   }
   DPY        = scanpar("dpy",0);
   if (DPY > 0) DaysPerYear = DPY;
   Cyclical   = scanpar("cyclical",Cyclical);
   if (Cyclical) Cyclical = 1;
   if (VerType == 0 || HorType == 0)
   {
      VerType  = scantex("vtype" ,"ps");   // 1. char is default value (p)
      HorType  = scantex("htype" ,"gsfz"); // 1. char is default value (g)
   }
   Multi    = scanpar("multi" ,0);
   LevelFactor = scanreal("levelfactor",1.0);
   if (NetCDF == 0) NetCDF = scanpar("netcdf",0);
   HeadOu[6]    = scanpar("head7",0);
   mars         = scanpar("mars" ,0);
   FirstMonth   = scanpar("first",1);
   LastMonth    = scanpar("last",12);
   PlanetRadius = scanreal("radius",EARTH_RADIUS);
   Grav         = scanreal("gravity",EARTH_GRAV);
   SigmaTop     = scanreal("sigmatop",0.0);
   vct[SigLevs] = SigmaTop;
   if (FirstMonth < 1) FirstMonth = 1;
   if (LastMonth > 12) LastMonth = 12;
   if (LastMonth < FirstMonth) LastMonth = FirstMonth;

   if (VerType == 's')
   {
      switch (HorType)
      {
         case 's': OutRep = HYB_SPEC; break;
         case 'f': OutRep = HYB_FOUR; break;
         case 'z': OutRep = HYB_ZONM; break;
         case 'g': OutRep = HYB_GRID; break;
      }
   }

   if (VerType == 'p')
   {
      switch (HorType)
      {
         case 's': OutRep = PRE_SPEC; break;
         case 'f': OutRep = PRE_FOUR; break;
         case 'z': OutRep = PRE_ZONM; break;
         case 'g': OutRep = PRE_GRID; break;
      }
   }

   if (Mean   == 0) Mean   = scanpar("mean"  ,1);
   if (Mean && OutRep < HYB_GRID)
   {
      sprintf(tb,"      mean =        1  ignored for spectral data on sigma levels");
      LeftText(tb);
      Mean = 0;
   }

   if (Multi) --Multi;

   if (mars) {
      Grav         = MARS_GRAV;
      PlanetRadius = MARS_RADIUS;
      RD           = MARS_RD;
   }

   scantime();
   scancode();
   DashLine();
   if (VerType == 's')   /* model levels */
   {
      scanmol();
      mom[0] = 1; // surface arrays are selected always
      if (nrml) /* Sigma levels explicitely given */
      {
         nrml = 0;
         for (i=1 ; i <= SigLevs ; ++i)
         {
            if (mom[i])
            {
               level[nrml] = mol[nrml] = i;
               nrml++;
            }
         }
         SigLevs = OutLevs = nrml;
      }
      else   /* No sigma levels specified -> select all */
      {
         OutLevs = nrml = SigLevs;
         for (i=0 ; i < OutLevs ; ++i)
         {
            level[i] = mol[i] = mom[i+1] = i+1;
         }
      }
      for (i=0 ; i < OutLevs ; ++i)
      {
         jind = mol[i] - 1;
         LevelUnits[i] = mol[i];
         SigmaF[i]     = (int)(500000.0 * (vct[jind+1] - vct[jind])); // sigma * 1E6
         sprintf(tb,"Level %4d = %10.6f",i+1,
         0.5 * (vct[SigLevs+mol[i]]+vct[SigLevs+mol[i]+1]));
         LeftText(tb);
      }
   }
   else             /* pressure levels */
   {
      scanhPa();
      if (nrpl)
      {
         OutLevs = nrpl;
         for (i=0 ; i < OutLevs ; ++i)
         {
            level[i] = 100.0 * hPa[i];
            if (hPa[i] <    0.0) Abort("pressure level < 0.0 is illegal");
//          if (hPa[i] > 2000.0) Abort("pressure level > 2000 hPa is illegal");
         }
      }
      else
      {
         OutLevs = nrpl = 10;
         for (i=0 ; i < OutLevs ; ++i)
         {
            hPa[i] = (i+1) * 100.0;
            level[i] = 100.0 * hPa[i];
         }
      }
      for (i=0 ; i < OutLevs ; ++i)
      {
         LevelUnits[i] = (int)(LevelFactor * 100.0 * hPa[i] + 0.5);
         sprintf(tb,"Level %4d = %14.4f hPa",i+1,hPa[i]);
         LeftText(tb);
      }
   }
   DashLine();
   scanattributes();
   for (i=0 ; i < nattr ; ++i)
   {
      sprintf(tb,"NetCDF attribute[%2d] :%s = \"%s\" ;",
              i,AttrNam[i],AttrVal[i]);
      LeftText(tb);
   }
}

void dimcalc(void)
{
   PolyDisk   = Lats == 2048;  // Currently T1365 only
   Waves      = Truncation + 1;
   Fouriers   = Waves * 2;
   DimSP      = (Truncation + 1) * (Truncation + 2);
   DimFC      = Lats * Fouriers;
   DimGP      = Lats * Lons;
   DimGG      = Gats * Gons;
   Dim3GP     = SigLevs * DimGP;
   Dim3GG     = SigLevs * DimGG;
   Dim3FC     = SigLevs * DimFC;
   Dim3SP     = SigLevs * DimSP;
   DimSP_half = DimSP / 2;
   DimAB      = MAX(DimGG,DimGP) + MAX(Lons,Gons);

   DashLine();
   if (VerType == 's') LeftText("Vertical   Type = Sigma                [S]");
   if (VerType == 'p') LeftText("Vertical   Type = Pressure             [P]");
   if (HorType == 's') LeftText("Horizontal Type = Spherical Harmonics  [S]");
   if (HorType == 'f') LeftText("Horizontal Type = Fourier Coefficients [F]");
   if (HorType == 'z') LeftText("Horizontal Type = Zonal Means          [Z]");
   if (HorType == 'g')
   {
      if (Lons == Gons && Lats == Gats)
      {
         LeftText("Horizontal Type = Gaussian Grid        [G]");
         GaussianOutput = 1;
      }
      else
      {
         LeftText("Horizontal Type = Grid (Lons x Lats)   [G]");
         GaussianOutput = 0;
      }
   }
   DashLine();

   Record_double    = new double[DimAB];
   Record_float     = (float *)Record_double;
   Record_int       = (int   *)Record_double;
   Record_short     = (unsigned short *)Record_double;
   Record_char      = (char  *)Record_double;

   CosPhi           = new double[Lats];
   RevCosPhi        = new double[Lats];
   DerivationFactor = new double[Lats];
}

/* ------------------------------------------------------ */
/* Check file OutRep and decode level number and truncation */
/* ------------------------------------------------------ */

void AnalyzeFile(void)
{
   int  i;
   LONG fcb,fce;       /* Fortran Record Control Words */
   char Id[8];
   char tb[COLS+2];
   
   union EndianCheck
   {
      char b[sizeof(int)];
      int  i;
   } ec;

   ec.i = 8;
   CoreBigEndian = ec.b[0] == 0;

   fread(Id,1,8,fpi);      // Read first 8 bytes
   FileBigEndian =   Id[0] == 0;
   Endian = CoreBigEndian != FileBigEndian;

   if (FileBigEndian)
   {
      if (Id[3] == 0) LongFCW = 1;
   }
   else
   {
      if (Id[4] == 0) LongFCW = 1;
   }

   rewind(fpi);
   HeadSize = fcb = ReadFCW();
   if (fcb != 8 && fcb != 32 && fcb != 64) Abort("Not a PUMA/PLASIM file");

   CenterText("Found PUMA/Planet Simulator data set");
   if (CoreBigEndian) CenterText("Running on  BIG   endian machine");
   else               CenterText("Running on little endian machine");
   if (FileBigEndian) CenterText("File is in  BIG   endian  format");
   else               CenterText("File is in little endian  format");
   if (Endian) CenterText("Endian swap activated");
   if (LongFCW) CenterText("Record control words have 64 bits");
   else         CenterText("Record control words have 32 bits");
   if (fcb ==  8) CenterText("Header format: PUMA-II       ");
   if (fcb == 32) CenterText("Header format: Service 32 bit");
   if (fcb == 64) CenterText("Header format: Service 64 bit");
   BlankLine();

   if (fcb == 8)
   {
      HeadSize = 32;
      i = fread(Id,1,8,fpi);
      if (strncmp(Id,"PUMA-II ",8)) Abort("PUMA-II header missing");
   
      fce = ReadFCW();
      if (check_fcw(fcb,fce)) Abort("Wrong FORTRAN control word after PUMA header");
   
      Truncation = ReadINTRecord();
      Gats = ReadINTRecord();
      AllLevs = SigLevs = ReadINTRecord();
      if (SigLevs < 1 || SigLevs > 1000) Abort("Illegal value for Level");
      SingleLevel = (SigLevs == 1);
      nvct  = 2 * SigLevs + 2;
   
      /* Check length of sigma vector and determine precision */
   
      fcb = ReadFCW();
      RealSize = fcb / SigLevs; // Should be float (4) or double (8)
      if (RealSize != sizeof(float) && RealSize != sizeof(double))
         Abort("FCW error on level record");
      for (i = 0 ; i < SigLevs ; i++)
      {
         if (RealSize == sizeof(float)) vct[i+SigLevs+2] = ReadFLOAT();
         else                           vct[i+SigLevs+2] = ReadDOUBLE();
      }
      fce = ReadFCW();
      if (check_fcw(fcb,fce)) Abort("FCW mismatch on level record");
   
   
      /* Header = "PUMA-II " Truncation   Lats      SigLevs      sigmah     */
      /*             2-FCW   FCW-1-FCW  FCW-1-FCW FCW-1-FCW FCW-SigLevs-FCW */
   
      HeaderWords = SigLevs * RealSize / 4 + 5 + 9 * (1 + LongFCW);
   
      ReadHeaderRecord();
      if (HeadIn[7] > 100) DaysPerYear = HeadIn[7];
   }
   else
   {
      rewind(fpi);
      ReadHeaderRecord();
      if (HeadIn[0] != 333) Abort("Initial code 333 was not found");
      Gats = HeadIn[5];
      AllLevs = SigLevs = HeadIn[6];
      SingleLevel = (SigLevs == 1);
      nvct  = 2 * SigLevs + 2;
      Truncation = HeadIn[7];

      /* Check length of array and determine precision */
   
      fcb = ReadFCW();
      RealSize = fcb / (Gats * Gats * 2); // Should be float (4) or double (8)
      fseek(fpi,-4*LongFCW,SEEK_CUR);
      if (RealSize != sizeof(float) && RealSize != sizeof(double))
         Abort("FCW error on first array");
      for (i = 0 ; i < SigLevs ; i++)
      {
         if (RealSize == sizeof(float)) vct[i+SigLevs+2] = ReadFLOAT();
         else                           vct[i+SigLevs+2] = ReadDOUBLE();
      }
      if (RealSize == sizeof(float)) DaysPerYear = ReadFLOAT();
      else                           DaysPerYear = ReadDOUBLE();
   }
   HeadSt = HeadIn;
   sprintf(tb,"Truncation                        = %6d",Truncation);
   CenterText(tb);
   sprintf(tb,"Latitudes                         = %6d",Gats);
   CenterText(tb);
   sprintf(tb,"Longitudes                        = %6d",Gats*2);
   CenterText(tb);
   sprintf(tb,"Sigma Levels                      = %6d",SigLevs);
   CenterText(tb);
   if (RealSize == 8)
      sprintf(tb,"Double precision data:      Bytes = %6d",RealSize);
   else if (RealSize == 4)
      sprintf(tb,"Single precision data:      Bytes = %6d",RealSize);
   else
      sprintf(tb,"Size of real data                 = %6d",RealSize);
   CenterText(tb);
   BlankLine();
   sprintf(tb,"Half Level             [p]         [sigma]");
   CenterText(tb);
   for (i = 0; i<nvct/2; i++)
   {
      sprintf(tb,"%10.1f  %14.4f  %14.4f",(float)i,vct[i],vct[i+nvct/2]);
      CenterText(tb);
   }
   StarLine();
   
   rewind(fpi);
}

const char *MoName[13] = {"Nix","Jan","Feb","Mar","Apr","May","Jun",
                                "Jul","Aug","Sep","Oct","Nov","Dec"};
void WriteGradsCtl()
{
   int i,j,yy,mm,dd,code,varcodes;
   FILE *fp;
   double DelLon;

   fp = fopen(gfile,"w");
   if (HorType == 'z') 
      fprintf(fp,"DSET ^%s\n",rfile);
   else
      fprintf(fp,"DSET ^%s\n",ofile);
   fprintf(fp,"UNDEF 9e09i\n");
   if (HorType == 'z')
      fprintf(fp,"XDEF 1 LINEAR 0 1\n");
   else
   {
      DelLon = 360.0 / Gons;
      fprintf(fp,"XDEF %5d LINEAR 0.0 %14.8f\n",Gons,DelLon);
   }
   fprintf(fp,"OPTIONS YREV\n");
   fprintf(fp,"YDEF %5d LEVELS\n",Gats);

   for (j=Gats-1 ; j >= 0 ; j-=8)
   {
      for (i=j ; i >= 0 && i> j-8 ; --i)
         fprintf(fp,"%14.8f",Outlat->Phi[i]);
      fprintf(fp,"\n");
   }
   if (HorType != 'z')
   {
      fprintf(fp,"OPTIONS SEQUENTIAL\n");
      fprintf(fp,"XYHEADER %d\n",40 + 8 * LongFCW);
   }
   if (OutLevs < 2)
   {
      if (VerType == 'p') fprintf(fp,"ZDEF 1 LINEAR %14.8f 1\n",hPa[0]);
      else                fprintf(fp,"ZDEF 1 LINEAR %d 1\n",mol[0]);
   }
   else
   {
      fprintf(fp,"ZDEF %d LEVELS\n",OutLevs);
      if (VerType == 'p')
      {
         
         for (j=0 ; j < OutLevs ; j+=8)
         {
            for (i=j ; i < OutLevs && i < j+8 ; ++i)
            if (HorType == 'z' && hPa[0] < hPa[OutLevs-1])
               fprintf(fp,"%14.8f",hPa[OutLevs-1-i]);
            else
               fprintf(fp,"%14.8f",hPa[i]);
            fprintf(fp,"\n");
         }
         if (HorType == 'z' && hPa[0] < hPa[OutLevs-1])
            fprintf(fp,"OPTIONS ZREV\n");
      }
      else
      {
         for (j=0 ; j < OutLevs ; j+=8)
         {
            for (i=j ; i < OutLevs && i < j+8 ; ++i)
            if (HorType == 'z' && mol[0] < mol[OutLevs-1])
               fprintf(fp,"%10d",mol[OutLevs-1-i]);
            else
               fprintf(fp,"%10d",mol[i]);
            fprintf(fp,"\n");
         }
         if (HorType == 'z' && mol[0] < mol[OutLevs-1])
            fprintf(fp,"OPTIONS ZREV\n");
      }
   }
   yy = HeadSt[2] / 10000;
   mm = HeadSt[2] / 100 % 100;
   dd = HeadSt[2] % 100;

   if (Mean)
      fprintf(fp,"TDEF %d LINEAR 00:00z%d%s%d 1mo\n",
         MonthCount,dd,MoName[mm],yy);
   else if (DeltaMn > 0)
      fprintf(fp,"TDEF %d LINEAR 00:00z%d%s%d %dmn\n",
         TermCount,dd,MoName[mm],yy,DeltaMn);
   else if (DeltaHr > 0)
      fprintf(fp,"TDEF %d LINEAR 00:00z%d%s%d %dhr\n",
         TermCount,dd,MoName[mm],yy,DeltaHr);
   else
      fprintf(fp,"TDEF %d LINEAR 00:00z%d%s%d %ddy\n",
         TermCount,dd,MoName[mm],yy,DeltaDy);


   varcodes = 0;
   for (code = 0 ; code < CODES ; ++code)
      if (All[code].selected) ++varcodes;
   fprintf(fp,"VARS %d\n",varcodes);
   for (code = 0 ; code < CODES ; ++code)
      if (All[code].selected)
         fprintf(fp,"%s %d 99 %s\n",
            All[code].Id,All[code].plev,All[code].Na);
   fprintf(fp,"ENDVARS\n");
   fclose(fp);
}



int main(int argc, char *argv[])
{
   int i,l;
   int wdim;
   char tb[COLS+2];
   time_t StartTime;

   StartTime = time(NULL);

   /*********************/
   /* print information */
   /*********************/
   
   StarLine();
   CenterText(V0);
   CenterText(V1);
   DashLine();
   CenterText(V2);
   CenterText(V3);
   CenterText(V4);
   StarLine();
   sprintf(tb,"Run started at %s",ctime(&StartTime));
   strcpy(tb+strlen(tb)-1," local time");
   CenterText(tb);
   StarLine();
   
   InitAll();

   /***********************/
   /* options & filenames */
   /***********************/

   for (i = 1 ; i < argc ; ++i) {
      if (argv[i][0] == '-') {
         if      (argv[i][1] == 'c') {PrintCodes(); exit(1);}
         else if (argv[i][1] == 'd') Debug      =  1;
         else if (argv[i][1] == 'v') Debug      =  1;
         else if (argv[i][1] == 'n') NetCDF     =  1;
         else if (argv[i][1] == 'g') Grads      =  1;
         else if (argv[i][1] == 'm') Mean       =  1;
         else if (argv[i][1] == 'p') PolyCreate =  1;
         else if (argv[i][1] == 'i') GaussGrid  =  1;
         else if (argv[i][1] == 'r') GaussGrid  =  0;
         else if (argv[i][1] == 's') SaveMemory =  1;
         else Usage();
      }
      else if (ifile[0] == '\0') strcpy(ifile,argv[i]);
      else if (ofile[0] == '\0') strcpy(ofile,argv[i]);
      else if (strcmp("Debug",argv[i]) == 0) Debug = 1;
      else Usage();
   }

   if (NetCDF) Grads = 0;

   if (ifile[0] == '\0' || ofile[0] == '\0') {
      printf("*** Missing filename ***\n");
      Usage();
   }

   /*******************/
   /* open input file */
   /*******************/
   
   fpi = fopen(ifile,"rb");
   if (fpi == 0) {
      printf("could not open input file %s\n",ifile);
      exit(1);
   }

   /******************/
   /* pre-processing */
   /******************/

   AnalyzeFile();

   Gons = Gats * 2;
   Lats = Gats;
   Lons = Gons;

   /*******************/
   /* initializations */
   /*******************/

   parini();

   l = strlen(ofile);
   if (NetCDF) // Add ".nc" extension for NetCDF filename if not specified 
   {
      if (l < 4 || strcmp(ofile+l-3,".nc")) strcat(ofile,".nc");
   }
   else        // Add ".srv" extension for Service format if not specified
   {
      if (l > 4 && !strcmp(ofile+l-4,".srv")) l -= 4;
      ofile[l] = 0;
      strcpy(gfile,ofile);
      strcpy(rfile,ofile);
      strcat(ofile,".srv");
      strcat(gfile,".ctl");
      strcat(rfile,".gra");
   }

   sprintf(tb," Input File: %s",ifile); LeftText(tb);
   sprintf(tb,"Output File: %s",ofile); LeftText(tb);
   if (Grads) {sprintf(tb,"Grads  File: %s",gfile); LeftText(tb);}
   if (Grads && HorType == 'z')
      {sprintf(tb,"Grads  Data: %s",rfile); LeftText(tb);}
   if (Debug) LeftText("Debug on");

#ifndef NETCDF_OUTPUT
   if (NetCDF) Abort("This executable was not compiled for NetCDF");
#endif

   /********************/
   /* open output file */
   /********************/

   if (!NetCDF)
   {
      gp = fopen(ofile,"wb");
      if (gp == 0)
      {
         printf("could not open output file <%s>\n",ofile);
         exit(1);
      }
   }

   dimcalc();
   
   HybSpec = new ServiceGrid(gp,DimSP,1);
   HybFour = new ServiceGrid(gp,Lats,Fouriers);
   HybGrid = new ServiceGrid(gp,Lons,Lats);
   HybSect = new ServiceSect(gp,Lats,OutLevs);

   /* FFT initialization */

   if (OutLevs <= SigLevs) wdim = (Lons + 3) * Lats * SigLevs;
   else                    wdim = (Lons + 3) * Lats * OutLevs  ;
   wfc  = new double[wdim];
   wgp  = new double[wdim];
   fft_set(Lons);

   filename = strrchr(ifile,'/');
   if (filename == 0) filename = ifile;
   else               filename++ ;
   HeadOu[7] = atol(filename);

   if (VerType == 'p'       &&
      (Divergence->selected || VeloPot->selected ||
        Vorticity->selected || StreamF->selected))
      SpecialUV = 1;

   Dependencies();

   // Check correct vertical coordinate

   if (GeopotHeight->selected && VerType != 'p')
   {
      printf("\n ****************** E R R O R ************************\n");
      printf(" * Geopotential height (156) requires pressure level *\n");
      printf(" *****************************************************\n");
      exit(1);
   }

   Geopotential->needed |= OutRep >= PRE_GRID
                        || SLP->needed || GeopotHeight->needed;
   if (OutRep) legini();

   if (PolyCreate) exit(0); /* Called for Legendre Polynomials only */

#ifdef NETCDF_OUTPUT
   if (NetCDF) NetOpen(ofile);
#endif

   PumaControl();
   if (gp) fclose(gp);

#ifdef NETCDF_OUTPUT
   if (NetCDF) NetClose();
#endif

   if (Grads) WriteGradsCtl();

   StarLine();
   if (DaysPerYear < 1) DaysPerYear = 360;
   if (DPM > 99) DaysPerYear = DPM * 12;
   sprintf(tb,"Using %d day calendar",DaysPerYear);
   CenterText(tb);
   sprintf(tb,"Data interval = %6.2f hours",DataStep * 24.0);
   CenterText(tb);
   StarLine();
   {
      double ut,st;
      struct rusage ru;
      getrusage(RUSAGE_SELF,&ru);
      ut = ru.ru_utime.tv_sec + 0.000001 * ru.ru_utime.tv_usec;
      st = ru.ru_stime.tv_sec + 0.000001 * ru.ru_stime.tv_usec;

      sprintf(tb,"User     time: %10.3lf seconds",ut);
      LeftText(tb);
      sprintf(tb,"System   time: %10.3lf seconds",st);
      LeftText(tb);
      sprintf(tb,"Total    time: %10.3lf seconds",ut+st);
      LeftText(tb);
      if (ru.ru_maxrss > 0)
      {
         sprintf(tb,"Memory  usage: %10.3lf MBytes",0.000001 * ru.ru_maxrss);
         LeftText(tb);
      }
      if (ru.ru_minflt > 0)
      {
         sprintf(tb,"Page reclaims: %10ld page",ru.ru_minflt);
         if (ru.ru_minflt != 1) strcat(tb,"s");
         LeftText(tb);
      }
      if (ru.ru_majflt > 0)
      {
         sprintf(tb,"Page   faults: %10ld page",ru.ru_majflt);
         if (ru.ru_majflt != 1) strcat(tb,"s");
         LeftText(tb);
      }
      if (ru.ru_nswap > 0)
      {
         sprintf(tb,"Swaps        : %10ld",ru.ru_nswap);
         LeftText(tb);
      }
      if (ru.ru_inblock > 0)
      {
         sprintf(tb,"Disk     read: %10ld block",ru.ru_inblock);
         if (ru.ru_inblock != 1) strcat(tb,"s");
         LeftText(tb);
      }
      if (ru.ru_oublock > 0)
      {
         sprintf(tb,"Disk    Write: %10ld block",ru.ru_oublock);
         if (ru.ru_oublock != 1) strcat(tb,"s");
         LeftText(tb);
      }
      StarLine();
   }
   return 0;
}
