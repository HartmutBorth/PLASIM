/*
    pumax - a GUI for PUMA & PLASIM
    (c) 2007 Edilbert Kirk (E.Kirk@gmx.de)

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

int Debug = 0; // set in initgui

#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>
#include <X11/keysym.h>

#define INT int
#define REAL float
#define INTXU unsigned int
#define INTXS int

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#define ROTATE_LEFT  (XK_Left)
#define ROTATE_RIGHT (XK_Right)

#define TOLELO 0x0001
#define TOLEHI 0x0002
#define TORILO 0x0004
#define TORIHI 0x0008
#define BOLELO 0x0010
#define BOLEHI 0x0020
#define BORILO 0x0040
#define BORIHI 0x0080
#define TOLEIN 0x0100
#define TOREIN 0x0200
#define BOLEIN 0x0400
#define BOREIN 0x0800

#define IPX(l,m,h) (VGAX * (x + (l-m) / (l-h)))
#define IPY(l,m,h) (VGAY * (y + (l-m) / (l-h)))

#define NUMWIN 9
#define NUMPAL 13

/* Picture types */

#define ISOHOR  0
#define ISOCS   1
#define ISOHOV  2
#define ISOTS   3
#define ISOTAB  4
#define ISOSH   5
#define ISOLON  6
#define ISOTRA  7
#define ISOCOL  8
#define ISOROT  9
#define MAPHOR 10
#define MAPTRA 11

/* Models */

#define PUMA   0

char *IsoNames[] =
{
   "ISOHOR",
   "ISOCS",
   "ISOHOV",
   "ISOTS",
   "ISOTAB",
   "ISOSH",
   "ISOLON",
   "ISOTRA",
   "ISOCOL",
   "ISOROT",
   "MAPHOR",
   "MAPTRA"
};

int IsoTypes = sizeof(IsoNames) / sizeof(char *);

char *PalNames[] =
{
   "AUTO",
   "U",
   "V",
   "T",
   "P",
   "Q",
   "MARST",
   "AMPLI",
   "VEG",
   "OCET",
   "OCES",
   "DCC",
   "DTDT"
};

int PalTypes = sizeof(PalNames) / sizeof(char *);

/* Map projections */

#define CYLINDRICAL 0
#define POLAR       1
#define AZIMUTHAL   2

#define MAXMAPS     3

char *ProNames[] =
{
   "CYLINDRICAL",
   "POLAR",
   "AZIMUTHAL"
};


/* Coordinates of button areas */

int ShowQueue = 0;
int ShowOnce  = 1;
int OffsetY = 20;

int Stop_XL;
int Stop_XH;
int Stop_YL;
int Stop_YH;

int Start_XL;
int Start_XH;
int Start_YL;
int Start_YH;

int Pause_XL;
int Pause_XH;
int Pause_YL;
int Pause_YH;

int Help_XL;
int Help_XH;
int Help_YL;
int Help_YH;

int Minus_XL;
int Minus_XH;
int Minus_YL;
int Minus_YH;

int FBWD_XL;
int FBWD_XH;
int FBWD_YL;
int FBWD_YH;

int Plus_XL;
int Plus_XH;
int Plus_YL;
int Plus_YH;

int FFWD_XL;
int FFWD_XH;
int FFWD_YL;
int FFWD_YH;

int Parc_XL;
int Parc_XH;
int Parc_XD;
int Parc_YL;
int Parc_YH;
int Parc_YD;

int Wbox_XL;
int Wbox_XH;
int Wbox_YL;
int Wbox_YH;

int Grid_XL;
int Grid_XH;
int Grid_YL;
int Grid_YH;

XPoint StopButton[4] = {{10,10},{40,10},{40,40},{10,40}};
XPoint StartButton[3] = {{50,10},{80,25},{50,40}};
XPoint PauseButton1[4] = {{95,10},{102,10},{102,40},{95,40}};
XPoint PauseButton2[4] = {{108,10},{115,10},{115,40},{108,40}};

#define PARCS 100

struct ParcStruct
{
   char Name[8];
   float Inc;
   float Min;
   float Max;
   float Val;
};

int Parcs;
int cpi = 1;

struct ParcStruct Parc[PARCS];

int screen_num;

Display *display;
unsigned int ScreenW,ScreenH,ScreenBot;
int SmallScreen;

static char *progname = "PUMA";

int BigEndian;
int ScreenOffset;
int Model;
int CowX = -1;
int CowY = -1;
int CowW;
int CowH;
int NumWin = NUMWIN;
int WinCols = 3;
int WinRows = 3;
int DimT = 100; // Default length of time axis
int DimX;
int DimY;
int DimZ;
int DimXY;
int Uindex = -1;
int Vindex = -1;
int ScreenD;
int OffX;
int OffY = 0;
INTXU WinXSize;
INTXU WinYSize;
int InXSize;
int InYSize;
int Latitudes;
int Lines;
int SizeChanged;
int Shutdown;
int Grid = 1;
int GridLabel = 1;
int WhitePix;
int BlackPix;
int *Flag;
int FlagSize;
int nstep;
int MRpid;
int MRnum;
int wto;
int FirstKey;
int LastKey;
int SymsPerKey;

KeySym *KeyMap;

char wtrun[8];

unsigned Seed;

float RotInc = 0.2;
int LineCo[NUMPAL];
int MapPro[NUMWIN];
int HovInx[NUMWIN];
int Indez[NUMWIN];
int MaxZ[NUMWIN];
int RedrawFlag[NUMWIN];
int LevelChanged[NUMWIN];
int RotLon[NUMWIN];
double AutoDelta[NUMWIN];
double AutoXdelt[NUMWIN];
double AutoLo[NUMWIN];
REAL *TSdata[NUMWIN];
REAL *Dmin[NUMWIN];
REAL *Dmax[NUMWIN];
REAL *SpeedScale;
XPoint *TSxp[NUMWIN];


#define MAXPAR 16384
float pax[MAXPAR];
float pay[MAXPAR];
long  pal[MAXPAR];

int Delpar = 8; /* Interval for particle injection */
int ParInd = 0;
int ColInd = 0;
long TracerColor;

#define NUMSCALARS 9

char TSName[NUMSCALARS][80];
char TSubsc[NUMSCALARS][80];
char TSUnit[NUMSCALARS][80];
char TScale[NUMSCALARS][80];

struct PixStruct
{
   int DimX;
   int DimY;
   Pixmap Pix;
} WinPixMap[NUMWIN];

Pixmap pix;

#define NUMARRAYS 100

struct ArrayStruct
{
   char Name[80];
   float *Data;
   int DimX;
   int DimY;
   int DimZ;
   int Flag;
} Array[NUMARRAYS];

int NumArrays;

struct WinAttStruct
{
   unsigned int x;
   unsigned int y;
   unsigned int w;
   unsigned int h;
   int active;
   char array_name[80];
   int Plot;
   int Palette;
} WinAtt[NUMWIN];

int BorderWidth =  0;
int WM_top_area =  0; // Window manager top    area
int WinLM       =  0;
int WinRM       =  0;
int WinTM       = 24;
int WinBM       =  0;

struct ColorStrip
{
   double  Lo;
   double  Hi;
   char *Name;
   unsigned long pixel;
};

#define AUTOCOLORS 14
struct ColorStrip Autostrip[AUTOCOLORS+1] =
{
   {0.0,0.0,"red"},
   {0.0,0.0,"OrangeRed"},
   {0.0,0.0,"Orange"},
   {0.0,0.0,"yellow"},
   {0.0,0.0,"GreenYellow"},
   {0.0,0.0,"green"},
   {0.0,0.0,"aquamarine"},
   {0.0,0.0,"cyan"},
   {0.0,0.0,"SkyBlue"},
   {0.0,0.0,"blue"},
   {0.0,0.0,"orchid"},
   {0.0,0.0,"magenta"},
   {0.0,0.0,"violet"},
   {0.0,0.0,"DarkViolet"},
   {0.0,0.0,NULL}
};

struct ColorStrip Ustrip[] =
{
   {-99.0,-10.0,"red"},
   {-10.0, -5.0,"OrangeRed"},
   { -5.0,  0.0,"Orange"},
   {  0.0,  5.0,"yellow"},
   {  5.0, 10.0,"GreenYellow"},
   { 10.0, 15.0,"green"},
   { 15.0, 20.0,"aquamarine"},
   { 20.0, 25.0,"cyan"},
   { 25.0, 30.0,"SkyBlue"},
   { 30.0, 35.0,"blue"},
   { 35.0, 40.0,"orchid"},
   { 40.0, 45.0,"magenta"},
   { 45.0, 50.0,"violet"},
   { 50.0, 99.0,"DarkViolet"},
   {  0.0,  0.0,NULL}
};

struct ColorStrip Qstrip[] =
{
   {-99.0,  0.0,"black"},
   {  0.0,  5.0,"OrangeRed"},
   {  5.0, 10.0,"Orange"},
   { 10.0, 15.0,"yellow"},
   { 15.0, 20.0,"GreenYellow"},
   { 20.0, 25.0,"green"},
   { 25.0, 30.0,"aquamarine"},
   { 30.0, 35.0,"cyan"},
   { 35.0, 40.0,"SkyBlue"},
   { 40.0, 45.0,"blue"},
   { 45.0, 50.0,"orchid"},
   { 50.0, 55.0,"magenta"},
   { 55.0, 60.0,"violet"},
   { 60.0, 99.0,"DarkViolet"},
   {  0.0,  0.0,NULL}
};

struct ColorStrip Vstrip[] =
{
   {-99.0,-10.0,"brown"},
   {-10.0, -8.0,"gold"},
   { -8.0, -6.0,"purple"},
   { -6.0, -4.0,"red"},
   { -4.0, -2.0,"OrangeRed"},
   { -2.0,  0.0,"Orange"},
   {  0.0,  2.0,"yellow"},
   {  2.0,  4.0,"GreenYellow"},
   {  4.0,  6.0,"green"},
   {  6.0,  8.0,"aquamarine"},
   {  8.0, 10.0,"cyan"},
   { 10.0, 99.0,"blue"},
   {  0.0,  0.0,NULL}
};

struct ColorStrip Tstrip[] =
{
   {-99.0,-50.0,"MidnightBlue"},
   {-50.0,-40.0,"RoyalBlue4"},
   {-40.0,-30.0,"RoyalBlue3"},
   {-30.0,-20.0,"RoyalBlue2"},
   {-20.0,-10.0,"RoyalBlue1"},
   {-10.0,  0.0,"violet"},
   {  0.0, 10.0,"IndianRed1"},
   { 10.0, 20.0,"IndianRed2"},
   { 20.0, 30.0,"IndianRed3"},
   { 30.0, 40.0,"IndianRed4"},
   { 40.0, 50.0,"red"},
   {  0.0,  0.0,NULL}
};

struct ColorStrip Kstrip[] =
{
   {100.0,210.0,"RoyalBlue4"},
   {210.0,220.0,"RoyalBlue3"},
   {220.0,230.0,"RoyalBlue2"},
   {230.0,240.0,"RoyalBlue1"},
   {240.0,250.0,"violet"},
   {250.0,260.0,"IndianRed1"},
   {260.0,270.0,"IndianRed2"},
   {270.0,280.0,"IndianRed3"},
   {280.0,290.0,"IndianRed4"},
   {290.0,400.0,"red"},
   {  0.0,  0.0,NULL}
};

struct ColorStrip Pstrip[] =
{
   {  00.0, 985.0,"RoyalBlue4"},
   { 985.0, 990.0,"RoyalBlue3"},
   { 990.0, 995.0,"RoyalBlue2"},
   { 995.0,1000.0,"RoyalBlue1"},
   {1000.0,1005.0,"violet"},
   {1005.0,1010.0,"IndianRed1"},
   {1010.0,1015.0,"IndianRed2"},
   {1015.0,1020.0,"IndianRed3"},
   {1020.0,1025.0,"IndianRed4"},
   {1025.0,9990.0,"red"},
   {  0.0,  0.0,NULL}
};

struct ColorStrip MarsTStrip[] =
{
   {-200.0,-90.0,"RoyalBlue4"},
   { -90.0,-80.0,"RoyalBlue3"},
   { -80.0,-70.0,"RoyalBlue2"},
   { -70.0,-60.0,"RoyalBlue1"},
   { -60.0,-50.0,"violet"},
   { -50.0,-40.0,"IndianRed1"},
   { -40.0,-30.0,"IndianRed2"},
   { -30.0,-20.0,"IndianRed3"},
   { -20.0,  0.0,"IndianRed4"},
   {   0.0,100.0,"red"},
   {   0.0,  0.0,NULL}
};

struct ColorStrip AmpliStrip[] =
{
   {   0.0,  1.0,"DarkBlue"},
   {   8.0,  9.0,"blue"},
   {   8.0,  9.0,"DarkGreen"},
   {   8.0,  9.0,"green"},
   {   8.0,  9.0,"yellow"},
   {   9.0, 10.0,"Orange"},
   {  10.0, 11.0,"OrangeRed"},
   {  11.0, 12.0,"red"},
   {   0.0,  0.0,NULL}
};

struct ColorStrip Vegstrip[] =
{
   {-999.9,0.001,"DarkBlue"},
   { 0.001, 10.0,"brown"},
   {  10.0, 20.0,"SeaGreen4"},
   {  20.0, 30.0,"SeaGreen3"},
   {  30.0, 40.0,"SeaGreen2"},
   {  40.0, 50.0,"SeaGreen1"},
   {  50.0, 60.0,"YellowGreen"},
   {  60.0, 70.0,"green4"},
   {  70.0, 80.0,"green3"},
   {  80.0, 90.0,"green2"},
   {  90.0,999.0,"green1"},
   {   0.0,  0.0,NULL}
};

struct ColorStrip Tstripoce[] =
{
   {-99.0, 0.00001,"Black"},
   { 0.00001, 4.0 ,"RoyalBlue4"},
   { 4.0 ,8.0,"RoyalBlue3"},
   { 8.0, 12.0,"RoyalBlue2"},
   { 12.0, 16.0,"RoyalBlue1"},
   { 16.0, 20.0,"violet"},
   { 20.0, 24.0,"IndianRed1"},
   { 24.0, 28.0,"IndianRed2"},
   { 28.0, 32.0,"IndianRed3"},
   { 32.0, 36.0,"IndianRed4"},
   { 36.0, 40.0,"red"},
   {  0.0,  0.0,NULL}
};

struct ColorStrip Sstripoce[] =
{
   {-99.0, 30.0,"Black"},
   { 30.0, 32.0 ,"RoyalBlue4"},
   { 32.0 ,32.5,"RoyalBlue3"},
   { 32.5, 33.0,"RoyalBlue2"},
   { 33.0, 34.0,"SeaGreen"},
   { 34.0, 34.5,"YellowGreen"},
   { 34.5, 35.0,"LightGreen"},
   { 35.0, 35.5,"Orange"},
   { 35.5, 36.0,"DarkOrange"},
   { 36.0, 36.5,"OrangeRed"},
   { 36.5, 40.0,"red"},
   {  0.0,  0.0,NULL}
};

struct ColorStrip DCCstrip[] =
{
   {-99.0, 0.0, "Black"},
   {  0.0, 10.0,"Yellow"},
   { 10.0, 20.0,"YellowGreen"},
   { 20.0, 30.0,"Green"},
   { 30.0, 40.0,"Aquamarine"},
   { 40.0, 50.0,"Cyan"},
   { 50.0, 60.0,"SkyBlue"},
   { 60.0, 70.0,"Blue"},
   { 70.0, 80.0,"Orchid"},
   { 80.0, 90.0,"Magenta"},
   { 90.0,100.0,"DarkViolet"},
   { 0.0, 0.0,NULL}
};

struct ColorStrip DTDTstrip[] =
{
   {-99.0,  -5.0,"brown"},
   {- 5.0,  -3.0,"gold"},
   { -3.0,  -1.5,"purple"},
   { -1.5,  -0.5,"red"},
   { -0.5,  -0.25,"OrangeRed"},
   { -0.25,  0.0,"Orange"},
   {  0.0,   0.25,"yellow"},
   {  0.25,  0.5,"GreenYellow"},
   {  0.5,   1.5,"green"},
   {  1.5,   3.0,"aquamarine"},
   {  3.0,   5.0,"cyan"},
   {  5.0, 99.0,"blue"},
   {  0.0,  0.0,NULL}
};

#define AMPLI_COLS 8

struct ColorStrip *Cstrip;
struct ColorStrip *Pallet[NUMPAL] =
{ Autostrip, Ustrip, Vstrip, Tstrip, Pstrip, Qstrip,MarsTStrip,AmpliStrip,
  Vegstrip, Tstripoce, Sstripoce, DCCstrip, DTDTstrip};

REAL VGAX;
REAL VGAY;
REAL CurVal;

REAL *Field;  // Pointer for ISO data
REAL *Ampli;  // Amplitudes of spherical harmonics
XColor *Acol; // Amplitude colors

Colormap colormap;

XColor xcolor1,xcolor2;
XColor Red,Green,Blue,Yellow,Grey,LightRed,DarkRed,LightBlue,DarkBlue;
XColor LightGreen,DarkGreen,Dummy;

unsigned long TSColor[10];

Window Win[NUMWIN];
int    win; // Current window
Window Cow; // Control bar
Window HelpWindow;

XTextProperty WinconName1;
XTextProperty WinconName3;
XTextProperty WinconPause;
XTextProperty WinconReady;
char *PauseTitle = "Run Paused";
char *ReadyTitle = "Press Start Button";

XEvent CowEvent;

XSizeHints WinSizeHints;
XSizeHints CowSizeHints;
int OutXSize,OutYSize;
int count;
XEvent report;
GC gc;

/*
char BigFontName[80] = "10x20";
char FixFontName[80] = "9x15bold";
char SubFontName[80] = "6x10";
*/
char BigFontName[80] = "9x15bold";
char FixFontName[80] = "8x13";
char SubFontName[80] = "6x10";

XFontStruct *FixFont;
XFontStruct *SubFont;
XFontStruct *BigFont;
int FixFontHeight;
int SubFontHeight;
int BigFontHeight;
int FixFontWidth;
int SubFontWidth;
int BigFontWidth;
int Paused = 0;

char *display_name = NULL;
XWMHints wm_hints;
XClassHint class_hints;
Atom Delwin;

char *WindowTitle[NUMWIN];
XTextProperty WindowName[NUMWIN];
XTextProperty HelpName;

char *mona[12] =
{
   "Jan","Feb","Mar",
   "Apr","May","Jun",
   "Jul","Aug","Sep",
   "Oct","Nov","Dec"
};

char *wena[7] = {"Mon","Tue","Wed","Thu","Fri","Sat","Sun"};

char datch[64];
char PlanetName[32] = "Earth";

char GUI_default[80] = "GUI.cfg";
char GUI_config[80] = "GUI_last_used.cfg";


#define PSDIM 30
char *pixelstar[PSDIM]  = {
"..............................",
"..............................",
"..............................",
"...**....**....**....**.......",
"...**....**....**....**.......",
"...**....**....**....**.......",
"...**....**....**....**.......",
"...**....**************.......",
"...**....**************.......",
"...**....**....**....**.......",
"...***...**....**....**.......",
"....******.....**....**.......",
".....****......**....**.......",
"..............................",
"..............................",
"..............................",
"..............................",
"..............................",
"..............................",
"..............................",
"..............................",
"..............................",
".......................*......",
"......................***.....",
"...................**.***.**..",
"...................**.***.**..",
"...................**.***.**..",
"...................****.****..",
"...................***...***..",
"...................***...***.."};

// Variables used for calculating "frames per second"

int fps;
int rmuf = 20;  // Rotating map update frequency [1/sec]
int rmui =  1;  // Rotating map update interval  [steps]
int SkipFreq;
int ThisSecond;
int LastSecond;
int LastStep;
int SecEvent;

int ndatim[6];  /* year, month, day, hour, minute, weekday */
int LastMinute; /* Last value of ndatim(4)                 */
int DeltaTime ; /* Computed timestep interval [min]        */

struct BMIstruct
{
   int   Size;
   int   Width;
   int   Height;
   short Planes;
   short Count;
   int   Compression;
   int   SizeImage;
   int   XPelsPerMeter;
   int   YPelsPerMeter;
   int   ClrUsed;
   int   ClrImportant;
};

struct BMIstruct ImageBMI;

struct MapImageStruct
{
   char   *d; // Bitmap data
   int     w; // Image width
   int     h; // Image height
   int     f; // Rotation factor
   float   l; // Reference longitude
   float   r; // Rotation speed [deg/step]
   XImage *X; // XImage structure
};

// Convert an RGB value to an X11 Pixel

unsigned long create_pixel(long red, long green, long blue)
{
   if (ScreenD == 24)       // 24 bit true color
   {
       return red | green << 8 | blue << 16;
   }
   else if (ScreenD == 16)  // 16 bit RGB 565
   {
      return red >> 3 | (green >> 2) << 5 | (blue >> 3) << 11;
   }
   else return 0;
}

struct MapImageStruct MapHR; // Hires (2560x1280) Earth image
struct MapImageStruct MapLR[NUMWIN];

void ScaleImage(struct MapImageStruct *s, struct MapImageStruct *d)
{
   int a,b,x,y;
   unsigned long px;
   unsigned long h,w;
   double f,g;

   if (d->X) XDestroyImage(d->X);
   w = (d->w + 7) & 0xFFF8;
   h =  d->h;
   d->d = malloc(4 * w * h);
   d->X = XCreateImage(display,CopyFromParent,ScreenD,ZPixmap,0,
                  d->d,w,h,8,0);

   f = (double)s->w / (double)d->w;
   g = (double)s->h / (double)d->h;

   for (y = 0 ; y < d->h ; ++y)
   for (x = 0 ; x < d->w ; ++x)
   {
       a = x * f;
       b = y * g;
       XPutPixel(d->X,x,y,XGetPixel(s->X,a,b));
   }
}


void PolarImage(struct MapImageStruct *s, struct MapImageStruct *d)
{
   int a,b,x,y;
   unsigned long h,w;
   double f,lfa,xnp,ynp,xsp,ysp,dx,dy;

   if (d->X) XDestroyImage(d->X);
   w = (d->w + 7) & 0xFFF8;
   h =  d->h;
   d->d = malloc(4 * w * h);
   d->X = XCreateImage(display,CopyFromParent,ScreenD,ZPixmap,0,
                  d->d,w,h,8,0);

   f   = (double)s->h / (double)d->h;
   lfa = (s->w-1) / (2.0 * M_PI);
   xnp = (d->w) * 0.25;
   ynp = (d->h) * 0.50;
   xsp = (d->w) * 0.75;
   ysp = ynp;

   for (y = 0 ; y < d->h   ; ++y)
   {
      dy = ynp - y;
      for (x = 0 ; x < d->w/2 ; ++x) /* Northern hemisphere */
      {
         dx = xnp - x;
         a = atan2(dx,dy) * lfa;
         b = sqrt(dx * dx + dy * dy) * f;
         if (a < 0.0) a += (s->w);
         if (a >= 0 && a < s->w && b >= 0 && b <= s->h/2)
              XPutPixel(d->X,x,y,XGetPixel(s->X,a,b));
         else XPutPixel(d->X,x,y,BlackPix);
      }
      dy = ysp - y;
      for (x = d->w/2 ; x < d->w ; ++x) /* Southern hemisphere */
      {
         dx = x - xsp;
         a = atan2(dx,dy) * lfa;
         b = s->h - f * sqrt(dx * dx + dy * dy);
         if (a < 0.0) a += (s->w);
         if (a >= 0 && a < s->w && b >= 0 && b >= s->h/2 && b <= s->h)
              XPutPixel(d->X,x,y,XGetPixel(s->X,a,b));
         else XPutPixel(d->X,x,y,BlackPix);
      }
   }
}


void RevOrtho(double R, double lam0, double phi0, int xx, int yy, double *lam, double *phi)
{
   double rho,c,rhoq,x,y;

   x = xx;
   y = yy;

   rhoq = x*x + y*y;
   rho  = sqrt(rhoq);
   if (rho > R)
   {
      *phi = 999.9;
      *lam = 999.9;
   }
   else if (rho > R / 100000.0)
   {
      c   = asin(rho / R);
      *phi = asin(cos(c) * sin(phi0) + y * sin(c) * cos(phi0) / rho);
      *lam = lam0 + atan2(x * sin(c) , rho * cos(phi0) * cos(c) - y * sin(phi0) * sin(c));
   }
   else
   {
      *phi = phi0;
      *lam = lam0;
   }
   // printf("[%6.1lf / %6.1lf] --> (%6.1lf / %6.1lf)\n",x,y,*lam,*phi);
}

/* ==================================================== */
/* AzimuthalImage - Display map in azimuthal projection */
/* ==================================================== */

void AzimuthalImage(struct MapImageStruct *s, struct MapImageStruct *d)
{
   int lam;           // lambda pixel coordinate in source image
   int phi;           // phi    pixel coordinate in source image
   int x  ;           // x      pixel coordinate in destination image
   int y  ;           // y      pixel coordinate in destination image
   int dxc;           // pixel coordinate of centre     position
   int dxl;           // pixel coordinate of left  most position
   int dxr;           // pixel coordinate of right most position
   int dx ;           // centre relative x position
   int dy ;           // centre relative y position
   int p00;           // euqator
   int l00;           // reference longitude
   
   unsigned int dih;  // destination image height
   unsigned int diw;  // destination image width
   unsigned int dpw;  // destination image padded width

   double rad;        // pixel radius of new image
   double ysc;        // yscale = source height / destination height
   double rho;        // distance from centre
   double xpi;        // x scale factor = 2 * PI / width
   double ypi;        // y scale factor =     PI / height
   double xrf;
   double yrf;

   XImage *sX;        // source      image
   XImage *dX;        // destination image


   // Destroy old image structure inclusive data storage

   if (d->X) XDestroyImage(d->X);

   // Set width of new image

   diw = d->w;

   // Pad width of new image to a multiple of 8

   dpw = (diw + 7) & 0xFFF8;

   // Set height of new image

   dih =  d->h;

   // Allocate space for image data

   d->d = calloc(dpw * dih,4);

   // Create image structure

   dX = d->X = XCreateImage(display,CopyFromParent,ScreenD,ZPixmap,0,d->d,dpw,dih,8,0);
   sX = s->X;

   p00 = s->h >> 1;
   rad = dih  >> 1;
   dxc = diw  >> 1;
   dxl = MAX(dxc - rad,   0);
   dxr = MIN(dxc + rad, diw);
   xpi = s->w / M_PI * 0.5;
   ypi = s->h / M_PI;
   ysc = (double)s->h / (double)dih;
   l00 = (int)((d->l * s->w) / 360 + s->w/2) % s->w;
   xrf = (double)diw / pow(2.0,31),
   yrf = (double)dih / pow(2.0,31);

   // Paint some stars on the sky

   srandom(Seed);
   for (y = 0 ; y < (diw * dih) >> 8  ; ++y)
      XPutPixel(dX,xrf*random(),yrf*random(),WhitePix);

   for (y = 0 ; y < dih ; ++y)
   {
      dy  = y - rad;
      phi = y * ysc;
      for (x = dxl ; x < dxr ; ++x)
      {
         dx = x - dxc;
         rho = sqrt(dx * dx + dy * dy);
         if (rho < rad)
         {
            lam = l00 + xpi * atan2(dx / rad, cos(asin(rho / rad)));
            phi = p00 + ypi *  asin(dy / rad);
            XPutPixel(dX,x,y,XGetPixel(sX,lam,phi));
         }
      }
   }
}


void SwapIEEE16(char W[2])
{
   char B;

   B = W[0]; W[0] = W[1]; W[1] = B; 
}

void SwapIEEE32(char W[4])
{
   char B;

   B = W[0]; W[0] = W[3]; W[3] = B; 
   B = W[1]; W[1] = W[2]; W[2] = B; 
}

int ReadINT(FILE *fpi)
{
   int i,k;
   i = fread(&k,sizeof(k),1,fpi);
   if (BigEndian) SwapIEEE32((char *)&k);
   return k;
}

short ReadSHORT(FILE *fpi)
{
   short i,k;
   i = fread(&k,sizeof(k),1,fpi);
   if (BigEndian) SwapIEEE16((char *)&k);
   return k;
}

void ReadImage(struct MapImageStruct *ei, char *filename)
{
   char ch;
   int i,n,x,y,z;
   long r,g,b;
   int FileSize;
   int Reserved;
   int OffBits;
   int PicBytes;
   int ImgBytes;
   int PadWidth;
   int PadBytes;
   int bpp;
   int byr;
   FILE *fp;
   unsigned char *BuffImageData;

   if (!(fp = fopen(filename,"r"))) return;
   ch = fgetc(fp);
   if (ch != 'B') return;
   ch = fgetc(fp);
   if (ch != 'M') return;
   FileSize = ReadINT(fp);
   Reserved = ReadINT(fp);
   OffBits  = ReadINT(fp);

   if (Debug)
   {
      printf("Properties of %s:\n",filename);
      printf("FileSize   = %d\n",FileSize);
      printf("FileOffset = %d\n",OffBits);
   }

   ImageBMI.Size          = ReadINT(fp);
   ImageBMI.Width         = ReadINT(fp);
   ImageBMI.Height        = ReadINT(fp);
   ImageBMI.Planes        = ReadSHORT(fp);
   ImageBMI.Count         = ReadSHORT(fp);
   ImageBMI.Compression   = ReadINT(fp);
   ImageBMI.SizeImage     = ReadINT(fp);
   ImageBMI.XPelsPerMeter = ReadINT(fp);
   ImageBMI.YPelsPerMeter = ReadINT(fp);
   ImageBMI.ClrUsed       = ReadINT(fp);
   ImageBMI.ClrImportant  = ReadINT(fp);

   PadWidth = (ImageBMI.Width + 7) & 0xFFF8;
   bpp      =  ImageBMI.Count >> 3;
   PadBytes = (4 - ((ImageBMI.Width * bpp) % 4)) % 4 ;

   if (Debug)
   {
      printf("BMI Size     = %d\n",ImageBMI.Size);
      printf("BMI Width    = %d\n",ImageBMI.Width);
      printf("BMI Height   = %d\n",ImageBMI.Height);
      printf("BMI Planes   = %d\n",ImageBMI.Planes);
      printf("BMI Count    = %d\n",ImageBMI.Count);
      printf("BMI ClrUsed  = %d\n",ImageBMI.ClrUsed);
      printf("Pad Bytes    = %d\n",PadBytes);
      printf("ScreenD= %d\n",ScreenD);
   }

   PicBytes = (bpp * ImageBMI.Width + PadBytes) * ImageBMI.Height;
   ImgBytes = 4 * PadWidth * ImageBMI.Height;
   ei->d = calloc(ImgBytes,1);
   BuffImageData = malloc(PicBytes);
   n = fread(BuffImageData,PicBytes,1,fp);
   fclose(fp);  
   
   ei->X = XCreateImage(display,CopyFromParent,ScreenD,ZPixmap,0,
                  ei->d,PadWidth,ImageBMI.Height,8,0);

   for (y = 0 ; y < ImageBMI.Height ; ++y)
   for (x = 0 ; x < ImageBMI.Width  ; ++x)
   {
       i = bpp * x + (ImageBMI.Height - 1 - y) * (ImageBMI.Width*bpp+PadBytes);
       r = BuffImageData[i  ];
       g = BuffImageData[i+1];
       b = BuffImageData[i+2];
       XPutPixel(ei->X, x, y, create_pixel(r,g,b));
   }
   free(BuffImageData);
   ei->w = ImageBMI.Width;
   ei->h = ImageBMI.Height;
}


void ntocdat(void)
{
   if (ndatim[1] > 0)
   {
      if (ndatim[5] >= 0)
         sprintf(datch,"%s %2d-%s-%04d %2d:%02d",wena[ndatim[5]],ndatim[2],mona[ndatim[1]-1],
           ndatim[0],ndatim[3],ndatim[4]);
      else
      sprintf(datch,"%2d-%s-%04d %2d:%02d",ndatim[2],mona[ndatim[1]-1],
           ndatim[0],ndatim[3],ndatim[4]);
   }
}

/* ============================================== */
/* CharAlloc - Allocate space for character array */
/* ============================================== */

char *CharAlloc(int bytes,char *array_name)
{
   char *result = NULL;
   if (bytes > 0)
   {
      result = (char *)calloc(bytes, sizeof(char));
      if (!result || Debug)
      {
         if (!result) printf("*** Out of Memory *** on ");
         printf("CharAlloc:  at%10p %6d char for %-30.30s\n",
                result,bytes,array_name);
      }
   }
   return(result);
}

/* ======================================= */
/* IntAlloc - Allocate space for int array */
/* ======================================= */

int *IntAlloc(int words,char *array_name)
{
   int *result = NULL;
   if (words > 0)
   {
      result = (int *) calloc(words, sizeof(int));
      if (!result || Debug)
      {
         if (!result) printf("*** Out of Memory *** on ");
         printf("IntAlloc:   at%10p %6d int for %s\n",
               result,words,array_name);
      }
   }
   return(result);
}

/* =========================================== */
/* FloatAlloc - Allocate space for float array */
/* =========================================== */

float *FloatAlloc(int words, char *array_name)
{
   float *result = NULL;
   if (words > 0)
   {
      result = (float *) calloc(words, sizeof(float));
      if (!result || Debug)
      {
         if (!result) printf("*** Out of Memory *** on ");
         printf("FloatAlloc: at%10p %6d float for %s\n",
               result,words,array_name);
      }
   }
   return(result);
}

/* ========================================== */
/* SizeAlloc - Allocate space for sized array */
/* ========================================== */

void *SizeAlloc(int words, int Size, char *array_name)
{
   void *result = NULL;
   if (words > 0)
   {
      result = (void *) calloc(words, Size);
      if (!result || Debug)
      {
         if (!result) printf("*** Out of Memory *** on ");
         printf("SizeAlloc:  at%10p %6d void for %s\n",
               result,words,array_name);
      }
   }
   return(result);
}

void LoadFonts(void)
{
   if ((FixFont = XLoadQueryFont(display,FixFontName)) == NULL)
   {
      printf("%s: Cannot open %s font\n",progname,FixFontName);
      exit(-1);
   }
   if ((SubFont = XLoadQueryFont(display,SubFontName)) == NULL)
   {
      SubFont = FixFont;
   }
   if ((BigFont = XLoadQueryFont(display,BigFontName)) == NULL)
   {
      if ((BigFont = XLoadQueryFont(display,"10x20")) == NULL)
      {
         printf("%s: Cannot open %s font\n",progname,BigFontName);
         exit(-1);
      }
   }
   FixFontWidth  = XTextWidth(FixFont,"X",1);
   SubFontWidth  = XTextWidth(SubFont,"X",1);
   BigFontWidth  = XTextWidth(BigFont,"X",1);
   FixFontHeight = FixFont->ascent + FixFont->descent;
   SubFontHeight = SubFont->ascent + SubFont->descent;
   BigFontHeight = BigFont->ascent + BigFont->descent;
}

int AllocateColorCells(struct ColorStrip cs[])
{
   int i = 0;
   XColor xcolor1,xcolor2;

   while (cs[i].Name)
   {
      XAllocNamedColor(display,colormap,cs[i].Name,&xcolor1,&xcolor2);
      cs[i].pixel = xcolor1.pixel;
      ++i;
   }   
   return i;
}

void CalcButtonAreas(void)
{
   int of,ds,bh,fh,fw,j;

   of = OffsetY;
   ds = 10;
   bh = 30;
   fh = BigFontHeight;
   fw = BigFontWidth;

   /* 1. Stop Button */

   StopButton[0].y += OffsetY;
   StopButton[1].y += OffsetY;
   StopButton[2].y += OffsetY;
   StopButton[3].y += OffsetY;

   Stop_XL = ds;
   Stop_XH = Stop_XL + bh;
   Stop_YL = of + ds;
   Stop_YH = Stop_YL + bh;

   /* 2. Start Button */

   Start_XL = Stop_XH  + ds;
   Start_XH = Start_XL + bh;
   Start_YL = Stop_YL;
   Start_YH = Stop_YH;

   StartButton[0].y += OffsetY;
   StartButton[1].y += OffsetY;
   StartButton[2].y += OffsetY;

   /* 3. Pause Button */

   Pause_XL = Start_XH + 15;
   Pause_XH = Start_XH + 35;
   Pause_YL = Start_YL;
   Pause_YH = Start_YH;

   PauseButton1[0].y += OffsetY;
   PauseButton1[1].y += OffsetY;
   PauseButton1[2].y += OffsetY;
   PauseButton1[3].y += OffsetY;

   PauseButton2[0].y += OffsetY;
   PauseButton2[1].y += OffsetY;
   PauseButton2[2].y += OffsetY;
   PauseButton2[3].y += OffsetY;

   /* 4. Help  Button */

   Help_XL = 200;
   Help_XH = Help_XL + 5 * fw;
   Help_YL = Pause_YL;
   Help_YH = Pause_YH;

   /* 5. Minus & FBWD Button */

   Minus_XL = 300;
   Minus_XH = Minus_XL + fh;
   Minus_YL = Pause_YL + (bh - fh) / 2;
   Minus_YH = Minus_YL + fh;

   FBWD_XH = Minus_XL - ds/2;
   FBWD_XL = FBWD_XH - fh;
   FBWD_YL = Minus_YL;
   FBWD_YH = Minus_YH;

   /* 6. Parameter Change Area */

   Parc_XL = Minus_XH + ds;
   Parc_XH = Parc_XL  + fw * 14;
   Parc_YL = Minus_YL;
   Parc_YH = Minus_YH;
   Parc_XD = Parc_XH - Parc_XL;
   Parc_YD = Parc_YH - Parc_YL;

   /* 7. Plus & FFWD Button */

   Plus_XL = Parc_XH + ds;
   Plus_XH = Plus_XL + fh;
   Plus_YL = Minus_YL;
   Plus_YH = Minus_YH;

   FFWD_XL = Plus_XH + ds/2;
   FFWD_XH = FFWD_XL + fh;
   FFWD_YL = Plus_YL;
   FFWD_YH = Plus_YH;

   /* 8. Grid on/off box */

   Grid_XL = Stop_XL;
   Grid_XH = Grid_XL + FixFontHeight;
   Grid_YL = Stop_YH + 3 * BigFontHeight / 2;
   Grid_YH = Grid_YL + FixFontHeight;

   /* Window status and select boxes */

   Wbox_XL = 580;
   Wbox_XH = Wbox_XL + FixFontHeight;
   Wbox_YL =   0;
   Wbox_YH = NumWin  * FixFontHeight;
}


/* Hack for getting information of window decoration  */
/* especially borderwidth and titleheight             */
/* These are properties of the window manager and are */
/* usually adjustable by the user                     */

void TestWindow(int xx, int yy, int wi, int he, int *x, int *y,
                int *X, int *Y, unsigned int *W, unsigned int *H)
{
   unsigned int DepthReturn  = 0;
   unsigned int BorderReturn = 0;
   Window TW,Parent,Child;
   
   TW = XCreateWindow(display,RootWindow(display,screen_num),
        xx,yy,wi,he,                   // x,y,w,h
        0,CopyFromParent,InputOutput,  // border, depth, class
          CopyFromParent,0,NULL);      // visual, valuemask, attributes
   XSelectInput(display,TW,ExposureMask);
   XMapWindow(display,TW);           
   XMoveWindow(display,TW,xx,yy);
   XNextEvent(display,&CowEvent);            // Wait until WM mapped it
   XGetGeometry(display,TW,&Parent,x,y,W,H,&BorderReturn,&DepthReturn);
   XTranslateCoordinates(display,TW,Parent,0,0,X,Y,&Child);
   XDestroyWindow(display,TW);

   if (Debug)
   {
      printf("TW[%x][%x]: %4d / %4d   %4d x %4d x %2d [%2d:%2d]\n",
              (int)TW,(int)Parent,*x,*y,*W,*H,DepthReturn,*X,*Y);
   }
}


void CreateTestWindow(void)
{
   int X,Y,xp,yp;
   unsigned int W,H;
   
   // Open a test window with displacement x = 100 and y = 100
   // Test the returned coordinates for upper left corner
   // and compute left margin (border) and top margin (title bar)

   TestWindow(100,100,100,100,&xp,&yp,&X,&Y,&W,&H);

   // Left margin

   if (X >= 100 && X < 110) WinLM = X - 100;

   // Right margin

   WinRM = WinLM;

   // Top margin

        if (Y >= 100 && Y < 150) WinTM = Y - 100;
   else if (Y >    0 && Y <  50) WinTM = Y;

   // Bottom margin

   WinBM = WinLM;

   // Open a full screen test window
   // Test the returned coordinates for upper left corner
   // the y displacement should be WinTM plus Windowmanager menu bar

   TestWindow(0,0,ScreenW,ScreenH,&xp,&yp,&X,&Y,&W,&H);

   // Window manager menu bar on top

   if (yp > WinTM) WM_top_area = yp - WinTM;


   // Last usable scan line

   ScreenBot = WM_top_area + WinTM + H;
   if (ScreenBot > ScreenH) ScreenBot = ScreenH;

   // Are we running on a 27" iMac :-)

   if (ScreenW == 2560 && ScreenH > 1400) ScreenBot -= 100;

   if (Debug)
   {
      printf("Screen       = %4d x %4d x %2d\n",ScreenW,ScreenH,ScreenD);
      printf("ScreenBot    = %4d\n",ScreenBot);
      printf("Top reserved = %4d\n",WM_top_area);
      printf("Bot reserved = %4d\n",ScreenH-ScreenBot);
      printf("Left  margin = %4d\n",WinLM);
      printf("Right margin = %4d\n",WinRM);
      printf("Top   margin = %4d\n",WinTM);
      printf("Bot   margin = %4d\n",WinBM);
   }
}


void TrimCopy(char *Trim, char *Line, int n)
{
   int l;
   char *s;

   l = strlen(Line);
   Line[--l] = 0;                          // Remove linefeed
   while (Line[l-1] == ' ') Line[--l] = 0; // Remove trailing blanks
   s = Line;
   while (*s == ' ') ++s;                  // Remove leading blanks
   strncpy(Trim,s,n);
}


int ReadConfig(char *filename)
{
   FILE *fp;
   int i,w,l,s,p;
   char Line[256];
   char Plot[80];
   char Palette[80];
   char Projection[80];
   char Rotation[80];
   char *t;

   fp = fopen(filename,"r");
   if (!fp) return 1;

   w = -2;
   s = -1;
   p = -1;
   while (!feof(fp))
   {
      t = fgets(Line,256,fp);
      if (Line[0] == '#') continue;
      if (strncmp(Line,"[Window "    , 8) == 0) w = atoi(Line+8);
      if (strncmp(Line,"[Control"    , 8) == 0) w = -1;
      if (strncmp(Line,"[Scalar "    , 8) == 0) s = atoi(Line+8);
      if (strncmp(Line,"WinRows ="   , 9) == 0)
      {
         WinRows = atoi(Line+9);
         if (WinRows < 1) WinRows = 1;
         if (WinRows > 3) WinRows = 3;
         NumWin = WinCols * WinRows;
      }
      if (strncmp(Line,"WinCols ="   , 9) == 0)
      {
         WinCols = atoi(Line+9);
         if (WinCols < 1) WinCols = 1;
         if (WinCols > 3) WinCols = 3;
         NumWin = WinCols * WinRows;
      }
      if (strncmp(Line,"[Parameter " ,11) == 0)
      {
         p = atoi(Line+11);
         if (p >= Parcs) Parcs = p+1;
      }
      if (strncmp(Line,"Geometry:",9) == 0 && w >= 0 && w < NUMWIN)
      {
         sscanf(Line+9,"%d %d %d %d",&WinAtt[w].w,&WinAtt[w].h,
			             &WinAtt[w].x,&WinAtt[w].y);
         WinAtt[w].active = 1;
      }
      if (strncmp(Line,"Inactive:",9) == 0 && w >= 0 && w < NUMWIN)
      {
         sscanf(Line+9,"%d %d %d %d",&WinAtt[w].w,&WinAtt[w].h,
			             &WinAtt[w].x,&WinAtt[w].y);
         WinAtt[w].active = 0;
      }
      if (strncmp(Line,"Geometry:",9) == 0 && w == -1)
      {
         sscanf(Line+9,"%d %d %d %d",&CowW,&CowH,&CowX,&CowY);
      }
      if (strncmp(Line,"Title:",6) == 0 && w >= -1 && w < NUMWIN)
      {
         if (w == -1) TrimCopy(PlanetName,Line+6,sizeof(PlanetName)-1);
         else         TrimCopy(WindowTitle[w]+wto,Line+6,79-wto);
      }
      if (strncmp(Line,"Array:",6) == 0 && w >= 0 && w < NUMWIN)
         TrimCopy(WinAtt[w].array_name,Line+6,79);
      if (strncmp(Line,"Plot:",5) == 0 && w >= 0 && w < NUMWIN)
      {
         TrimCopy(Plot,Line+5,79);
	 for (i=0 ; i < IsoTypes ; ++i)
            if (!strcmp(Plot,IsoNames[i])) WinAtt[w].Plot = i;
      }
      if (strncmp(Line,"Projection:",11) == 0 && w >= 0 && w < NUMWIN)
      {
         TrimCopy(Projection,Line+11,79);
	 for (i=0 ; i < MAXMAPS ; ++i)
            if (!strcmp(Projection,ProNames[i])) MapPro[w] = i;
      }
      if (strncmp(Line,"Rotation factor:",15) == 0 && w >= 0 && w < NUMWIN)
      {
         TrimCopy(Rotation,Line+15,79);
         MapLR[w].f = atoi(Rotation);
      }
      if (strncmp(Line,"Palette:",8) == 0 && w >= 0 && w < NUMWIN)
      {
         TrimCopy(Palette,Line+8,79);
	 for (i=0 ; i < PalTypes ; ++i)
            if (!strcmp(Palette,PalNames[i])) WinAtt[w].Palette = i;
      }

      // Dimension for time axis

      if (strncmp(Line,"DimT:" ,5) == 0) DimT = atoi(Line+5);

      // Attributes of scalars for timeseries windows

      if (strncmp(Line,"Name:" ,5) == 0 && s >= 0 && s < NUMSCALARS)
         TrimCopy(TSName[s],Line+5,79);
      if (strncmp(Line,"Sub:"  ,4) == 0 && s >= 0 && s < NUMSCALARS)
         TrimCopy(TSubsc[s],Line+4,79);
      if (strncmp(Line,"Unit:" ,5) == 0 && s >= 0 && s < NUMSCALARS)
         TrimCopy(TSUnit[s],Line+5,79);
      if (strncmp(Line,"Scale:",6) == 0 && s >= 0 && s < NUMSCALARS)
         TrimCopy(TScale[s],Line+6,79);

      // Attributes of parameter in change menu

      if (strncmp(Line,"ParName:",8) == 0 && p >= 0 && p < PARCS)
         TrimCopy(Parc[p].Name,Line+8,6);
      if (strncmp(Line,"ParInc:" ,7) == 0 && p >= 0 && p < PARCS)
         Parc[p].Inc = atof(Line+7);
      if (strncmp(Line,"ParMin:" ,7) == 0 && p >= 0 && p < PARCS)
         Parc[p].Min = atof(Line+7);
      if (strncmp(Line,"ParMax:" ,7) == 0 && p >= 0 && p < PARCS)
         Parc[p].Max = atof(Line+7);
      if (strcmp(Parc[p].Name,"SELLON") == 0) Parc[p].Inc = 360.0 / (2.0 * Latitudes);
   }
   fclose(fp);
   return 0;
}


void CreateControlWindow(void)
{
   int X,Y,W,H,B,D,x,y;
   char Title1[256];
   char *WinconTitle1;
   char *WinconTitle3 = {"RUN finished - click on red stop button"};
   Window Rootwin;
   Window Child;

   strcpy(Title1,"Unknown model");
   if (Model == 0) // PUMA
   {
      if (MRpid < 0) 
         strcpy(Title1,"PUMA - KlimaCampus Hamburg");
      else
         sprintf(Title1,"Run %d: PUMA - KlimaCampus Hamburg",MRpid);
   }
   else if (Model == 1) sprintf(Title1,"SAM - KlimaCampus");
   else if (Model == 2)
   {
       if (MRpid < 0)
       sprintf(Title1,"Planet Simulator (%s) - KlimaCampus",PlanetName);
       else
       sprintf(Title1,"Run %d: Planet Simulator (%s) - KlimaCampus",MRpid,PlanetName);
   }
   else if (Model == 3) sprintf(Title1,"Planet Simulator / LSG - KlimaCampus");
   else printf("*** unknown model number %d in pumax ***\n",Model);
   WinconTitle1 = Title1;
   if (CowW < 920 || CowW > ScreenW) CowW = 920;
   if (NumWin < 5) CowH = 6 * FixFontHeight + 5;
   else            CowH = NumWin * FixFontHeight + 5;
   if (CowX < 0 || CowX > ScreenW - CowW) CowX = (ScreenW - CowW) / 2;
   if (CowY < 0 || CowY > ScreenBot - CowH - WinTM - WinBM) CowY = ScreenBot - CowH - WinTM - WinBM;
   CowSizeHints.flags      = PPosition | PSize | PMinSize | PMaxSize;
   CowSizeHints.min_width  = CowW;
   CowSizeHints.min_height = CowH;
   CowSizeHints.max_width  = ScreenW - WinLM - WinRM;
   CowSizeHints.max_height = CowH;
   if (MRnum == 2 && MRpid == 1) CowX += ScreenW;
   if (Debug) printf("Control Window %d/%d %dx%d\n",CowX,CowY,CowW,CowH);
   Cow = XCreateWindow(display,RootWindow(display,screen_num), // display, parent
               CowX,CowY,CowW,CowH,
               BorderWidth,CopyFromParent,InputOutput,         // border, depth, class
               CopyFromParent,0,NULL);                         // visual, valuemask, attributes
   XStringListToTextProperty(&WinconTitle1,1,&WinconName1);
   XStringListToTextProperty(&WinconTitle3,1,&WinconName3);
   XStringListToTextProperty(&PauseTitle,1,&WinconPause);
   XStringListToTextProperty(&ReadyTitle,1,&WinconReady);
   XSetWMProtocols(display,Cow,&Delwin,1);
   XSetWMProperties(display,Cow,&WinconName1,NULL,
		NULL,0,&CowSizeHints,&wm_hints,&class_hints);
   XSelectInput(display,Cow,ButtonPressMask | KeyPressMask | ExposureMask);
   XMapWindow(display,Cow);
}


void ShowStep(void)
{
   char Text[80];
   XSetFont(display, gc, BigFont->fid);
   XSetForeground(display,gc,BlackPix);
   XSetBackground(display,gc,Grey.pixel);
   ntocdat();
   if (ndatim[5] == 6) XSetForeground(display,gc,Red.pixel); // sunday
   XDrawImageString(display,Cow,gc,10,BigFontHeight,datch,strlen(datch));
   XSetForeground(display,gc,BlackPix);
   if (ShowQueue)
   {
      sprintf(Text,"%8d Events   ",XPending(display));
      XDrawImageString(display,Cow,gc,10,60 + BigFontHeight,Text,strlen(Text));
   }
   else if (fps > 1)
   {
      if (SkipFreq < 2) sprintf(Text,"%5d fps",fps);
      else sprintf(Text,"%5d fps [%d]",fps,SkipFreq);
      XDrawImageString(display,Cow,gc,10,60 + BigFontHeight,Text,strlen(Text));
   }
}


int FormatVal(char *Name, float V, char *Text)
{
   char Format[80];

   if (V > 9999.9 || V < -999.9)
   {
      sprintf(Text,"%6.6s = *****",Name);
      return strlen(Text);
   }
   else if (V > 999.99  ) strcpy(Format,"%6.6s =%6.1f");
   else if (V >  99.999 ) strcpy(Format,"%6.6s =%6.2f");
   else if (V >   9.9999) strcpy(Format,"%6.6s =%6.3f");
   else if (V >   0.0   ) strcpy(Format,"%6.6s =%6.4f");
   else if (V >  -9.9   ) strcpy(Format,"%6.6s =%6.3f");
   else if (V > -99.9   ) strcpy(Format,"%6.6s =%6.2f");
   else                   strcpy(Format,"%6.6s =%6.1f");

   sprintf(Text,Format,Name,V);
   return strlen(Text);
}


void ShowParcs(void)
{
   char Text[80];
   int len;
   float V;

   XSetFont(display, gc, BigFont->fid);
   len = FormatVal(Parc[cpi].Name,Parc[cpi].Val,Text);

   XSetBackground(display,gc,WhitePix);
   XSetForeground(display,gc,Red.pixel);
   XFillRectangle(display,Cow,gc,Minus_XL,Minus_YL,Minus_XH-Minus_XL,Minus_YH-Minus_YL);
   XFillRectangle(display,Cow,gc,FBWD_XL,FBWD_YL,FBWD_XH-FBWD_XL,FBWD_YH-FBWD_YL);
   XSetForeground(display,gc,WhitePix);
   XFillRectangle(display,Cow,gc,Minus_XL+3,(Minus_YH+Minus_YL)/2-1,Minus_XH-Minus_XL-5,3);
   XFillRectangle(display,Cow,gc,FBWD_XL+3,(FBWD_YH+FBWD_YL)/2-1,FBWD_XH-FBWD_XL-5,3);
   XSetForeground(display,gc,BlackPix);
   XDrawImageString(display,Cow,gc,Parc_XL,Parc_YH-3,Text,len);
   XSetForeground(display,gc,Green.pixel);
   XFillRectangle(display,Cow,gc,Plus_XL,Plus_YL,Plus_XH-Plus_XL,Plus_YH-Plus_YL);
   XFillRectangle(display,Cow,gc,FFWD_XL,FFWD_YL,FFWD_XH-FFWD_XL,FFWD_YH-FFWD_YL);
   XSetForeground(display,gc,WhitePix);
   XFillRectangle(display,Cow,gc,Plus_XL+3,(Plus_YH+Plus_YL)/2-1,Plus_XH-Plus_XL-5,3);
   XFillRectangle(display,Cow,gc,FFWD_XL+3,(FFWD_YH+FFWD_YL)/2-1,FFWD_XH-FFWD_XL-5,3);
   XFillRectangle(display,Cow,gc,(Plus_XH+Plus_XL)/2-1,Plus_YL+3,3,Plus_YH-Plus_YL-5);
   XFillRectangle(display,Cow,gc,(FFWD_XH+FFWD_XL)/2-1,FFWD_YL+3,3,FFWD_YH-FFWD_YL-5);
   XSetForeground(display,gc,BlackPix);
   XDrawRectangle(display,Cow,gc,Minus_XL,Minus_YL,Minus_XH-Minus_XL,Minus_YH-Minus_YL);
   XDrawRectangle(display,Cow,gc,FBWD_XL,FBWD_YL,FBWD_XH-FBWD_XL,FBWD_YH-FBWD_YL);
   XDrawRectangle(display,Cow,gc,Parc_XL-1,Parc_YL-1,Parc_XH-Parc_XL+1,Parc_YH-Parc_YL);
   XDrawRectangle(display,Cow,gc,Plus_XL,Plus_YL,Plus_XH-Plus_XL,Plus_YH-Plus_YL);
   XDrawRectangle(display,Cow,gc,FFWD_XL,FFWD_YL,FFWD_XH-FFWD_XL,FFWD_YH-FFWD_YL);

   /* Show next parameter */

   strcpy(Text,"              ");
   if (cpi < Parcs-1) FormatVal(Parc[cpi+1].Name,Parc[cpi+1].Val,Text);
   len = strlen(Text);

   XSetBackground(display,gc,WhitePix);
   XSetForeground(display,gc,BlackPix);
   XDrawImageString(display,Cow,gc,Parc_XL,2*Parc_YH-Parc_YL-3,Text,len);
   XDrawRectangle(display,Cow,gc,Parc_XL-1,Parc_YH-1,Parc_XH-Parc_XL+1,Parc_YH-Parc_YL);

   /* Show previous parameter */

   strcpy(Text,"              ");
   if (cpi > 0) FormatVal(Parc[cpi-1].Name,Parc[cpi-1].Val,Text);
   len = strlen(Text);

   XSetBackground(display,gc,WhitePix);
   XSetForeground(display,gc,BlackPix);
   XDrawImageString(display,Cow,gc,Parc_XL,Parc_YL-3,Text,len);
   XDrawRectangle(display,Cow,gc,Parc_XL-1,2*Parc_YL-Parc_YH-1,Parc_XH-Parc_XL+1,Parc_YH-Parc_YL);
   XSetFont(display, gc, FixFont->fid);

}


void CheckMark(int x, int y, int d)
{
   XDrawLine(display,Cow,gc,x+2,y+2,x+d-1,y+d-1);
   XDrawLine(display,Cow,gc,x+3,y+2,x+d-1,y+d-2);
   XDrawLine(display,Cow,gc,x+2,y+3,x+d-2,y+d-1);

   XDrawLine(display,Cow,gc,x+d-1,y+1,x+2,y+d-2);
   XDrawLine(display,Cow,gc,x+d-2,y+1,x+2,y+d-3);
   XDrawLine(display,Cow,gc,x+d-1,y+2,x+3,y+d-2);
}


void ShowWindowStatus(void)
{
   char *cp;
   char Text[80];
   int len,w,x,d;

   XSetFont(display, gc, FixFont->fid);
   XSetForeground(display,gc,BlackPix);
   XSetBackground(display,gc,Grey.pixel);
   x = Wbox_XL;
   d = FixFontHeight;
   for (w=0 ; w < NumWin ; ++w)
   {
      strncpy(Text,WindowTitle[w]+wto,40);
      Text[40] = 0;
      cp = strstr(Text,"Level");
      if (cp) *cp = 0;
      cp = strstr(Text,"Latitude");
      if (cp) *cp = 0;
      len = strlen(Text);
      XDrawImageString(display,Cow,gc,x+20,(w+1) * d,Text,len);
   }
   XSetForeground(display,gc,WhitePix);
   for (w=0 ; w < NumWin ; ++w)
   {
      XFillRectangle(display,Cow,gc,x,w*d,d,d);
   }
   XSetForeground(display,gc,BlackPix);
   for (w=0 ; w < NumWin ; ++w)
   {
      XDrawRectangle(display,Cow,gc,x,w*d,d,d);
   }
   for (w=0 ; w < NumWin ; ++w)
   if (Win[w]) CheckMark(x,w*d,d);
}


void ShowGridStatus(void)
{
   int x,y,d;

   d = FixFontHeight;
   x = Grid_XL;
   y = Grid_YL;
   XSetFont(display, gc, FixFont->fid);
   XSetForeground(display,gc,BlackPix);
   XSetBackground(display,gc,Grey.pixel);
   XDrawImageString(display,Cow,gc,Grid_XH+FixFontWidth,Grid_YH,"Grid",4);
   XSetForeground(display,gc,WhitePix);
   XFillRectangle(display,Cow,gc,x,Grid_YL,d,d);
   XSetForeground(display,gc,BlackPix);
   XDrawRectangle(display,Cow,gc,x,Grid_YL,d,d);

   if (Grid) CheckMark(x,y,d);
}

      
int RedrawControlWindow(void)
{
   char Text[80];
   int status,WinXSize,WinYSize,font_height,width,len;
   int i,j,x1,x2,y1,y2,y3;
   XWindowAttributes CurAtt;

   status = XGetWindowAttributes(display,Cow,&CurAtt);
   WinXSize = CurAtt.width;
   WinYSize = CurAtt.height;
   
   XSetForeground(display,gc,BlackPix);
   XSetBackground(display,gc,WhitePix);

   XSetWindowBackground(display,Cow,Grey.pixel);
   XClearWindow(display,Cow);

   /* Red stop button */

   XSetForeground(display,gc,Red.pixel);
   XFillPolygon(display,Cow,gc,StopButton,4,Convex,CoordModeOrigin);
   XSetForeground(display,gc,WhitePix);
   x1 = StopButton[0].x ;
   y1 = StopButton[0].y ;
   for (j=0 ; j < PSDIM ; ++j)
   for (i=0 ; i < PSDIM ; ++i)
   if (pixelstar[j][i] == '*')
      XDrawPoint(display,Cow,gc,i+x1,j+y1);
  
   XSetForeground(display,gc,LightRed.pixel);
   x1 = 10; x2 = 40 ; y1 = OffsetY + 10 ; y2 = OffsetY + 40;
   XDrawLine(display,Cow,gc,x2  ,y1-1,x2  ,y2  );
   XDrawLine(display,Cow,gc,x2+1,y1-2,x2+1,y2+1);
   XDrawLine(display,Cow,gc,x1-1,y1-1,x2  ,y1-1);
   XDrawLine(display,Cow,gc,x2+1,y2-2,x2+1,y1-2);
   XSetForeground(display,gc,DarkRed.pixel);
   XDrawLine(display,Cow,gc,x1-1,y1-1,x1-1,y2  );
   XDrawLine(display,Cow,gc,x1-2,y1-2,x1-2,y2+1);
   XDrawLine(display,Cow,gc,x1-1,y2  ,x2  ,y2  );
   XDrawLine(display,Cow,gc,x1-2,y2+1,x2+1,y2+1);

   // Green <play> button

   XSetForeground(display,gc,Green.pixel);
   XFillPolygon(display,Cow,gc,StartButton,3,Convex,CoordModeOrigin);
   x1 = 50; x2 = 80 ; y1 = OffsetY + 10 ; y2 = OffsetY + 25; y3 = OffsetY + 40;
   XSetForeground(display,gc,LightGreen.pixel);
   XDrawLine(display,Cow,gc,x1  ,y1-1,x2  ,y2  );
   XDrawLine(display,Cow,gc,x1  ,y1-2,x2+1,y2  );
   XSetForeground(display,gc,DarkGreen.pixel);
   XDrawLine(display,Cow,gc,x1-1,y1-1,x1-1,y3+1);
   XDrawLine(display,Cow,gc,x1-2,y1-2,x1-2,y3+2);

   // Blue <pause> button

   XSetForeground(display,gc,Blue.pixel);
   XFillPolygon(display,Cow,gc,PauseButton1,4,Convex,CoordModeOrigin);
   XFillPolygon(display,Cow,gc,PauseButton2,4,Convex,CoordModeOrigin);
   
   x1 = 95; x2 = 102; y1 = OffsetY + 10; y2 = OffsetY + 40;
   XSetForeground(display,gc,DarkBlue.pixel);
   XDrawLine(display,Cow,gc,x1-1,y1-1,x1-1,y2  );
   XDrawLine(display,Cow,gc,x1-2,y1-2,x1-2,y2+1);
   XDrawLine(display,Cow,gc,x1-1,y2  ,x2  ,y2  );
   XDrawLine(display,Cow,gc,x1-2,y2+1,x2+1,y2+1);
   XSetForeground(display,gc,LightBlue.pixel);
   XDrawLine(display,Cow,gc,x2  ,y1-1,x2  ,y2  );
   XDrawLine(display,Cow,gc,x2+1,y1-2,x2+1,y2+1);
   XDrawLine(display,Cow,gc,x1-1,y1-1,x2  ,y1-1);
   XDrawLine(display,Cow,gc,x2+1,y2-2,x2+1,y1-2);

   x1 = 108; x2 = 115; y1 = OffsetY + 10; y2 = OffsetY + 40;
   XSetForeground(display,gc,DarkBlue.pixel);
   XDrawLine(display,Cow,gc,x1-1,y1-1,x1-1,y2  );
   XDrawLine(display,Cow,gc,x1-2,y1-2,x1-2,y2+1);
   XDrawLine(display,Cow,gc,x1-1,y2  ,x2  ,y2  );
   XDrawLine(display,Cow,gc,x1-2,y2+1,x2+1,y2+1);
   XSetForeground(display,gc,LightBlue.pixel);
   XDrawLine(display,Cow,gc,x2  ,y1-1,x2  ,y2  );
   XDrawLine(display,Cow,gc,x2+1,y1-2,x2+1,y2+1);
   XDrawLine(display,Cow,gc,x1-1,y1-1,x2  ,y1-1);
   XDrawLine(display,Cow,gc,x2+1,y2-2,x2+1,y1-2);

   // Help button

   x1 = 200; x2 = x1 + 5 * BigFontWidth; y1 = OffsetY + 10; y2 = y1 + 30;
   XSetForeground(display,gc,Red.pixel);
   XFillRectangle(display,Cow,gc,x1,y1,x2-x1,y2-y1);
   XSetForeground(display,gc,BlackPix);
   XDrawRectangle(display,Cow,gc,x1-1,y1-1,x2-x1+1,y2-y1+1);
   XSetFont(display,gc,BigFont->fid);
   XSetForeground(display,gc,Yellow.pixel);
   XSetBackground(display,gc,Red.pixel);
   XDrawImageString(display,Cow,gc,x1+BigFontWidth/2,(y1+y2)/2+BigFont->ascent/2,"HELP",4);

   ShowStep();
   ShowParcs();
   ShowWindowStatus();
   ShowGridStatus();
   return 0;
}


void ClearTracer(void)
{
   int j;
   for (j=0 ; j < MAXPAR ; ++j)
   {
      pax[j] = -1.0;
      pay[j] = -1.0;
   }
}


int CheckEndianess(void)
{
   union EndianCheck
   {
      char b[sizeof(int)];
      int  i;
   } ec;

   ec.i = 8; 
   return (ec.b[0] == 0);
}


void initgui_(int *model, int *debug, int *lats, int *mrpid, int *mrnum)
{
   int x,y,w,h;
   int argc = 1;
   int i,j;
   unsigned long valuemask = 0; /* ignore XGCvalues and use defaults */
   XGCValues values;
   XEvent Event;
   struct timeval TimeVal;

   // Check for big endian computer
   // The bitmap format ".bmp" ist little endian

   BigEndian = CheckEndianess();

   // Set seed for random number generator from clock

   gettimeofday(&TimeVal,NULL);
   Seed = TimeVal.tv_sec;
   
   Model = *model;
   Debug = *debug;
   Latitudes = *lats;
   MRpid = *mrpid;
   MRnum = *mrnum;

   if (MRpid >= 0)
   {
      sprintf(wtrun,"[%d] ",MRpid);
      wto = strlen(wtrun);
   }

   if (Debug) printf("initgui(%d,%d)\n",Model,Debug);
   if ((display=XOpenDisplay(display_name)) == NULL)
   {
      fprintf(stderr,"%s: cannot connect to X server %s\n", 
              progname, XDisplayName(display_name));
      exit(1);
   }
   screen_num = DefaultScreen(display);
   ScreenD  = XDefaultDepth(display,screen_num);
   ScreenW  = DisplayWidth (display,screen_num);
   ScreenH  = DisplayHeight(display,screen_num);
   BlackPix = BlackPixel(display,screen_num);
   WhitePix = WhitePixel(display,screen_num);
   SmallScreen = ScreenH <= 768;
   if (SmallScreen)
   {
      WinRows = 2;
      NumWin  = WinCols * WinRows;
      Grid = 0;
   }

   if (MRnum == 2) ScreenW /= 2; /* Run with two GUI's */

   /* Fix defaults for multihead displays */

   if (ScreenH > ScreenW && ScreenH > 1200) ScreenH = 1200;
   if (ScreenW > 2 * ScreenH && ScreenW >= 2048) ScreenW /= 2;

   for (i = 0 ; i < NUMWIN ; ++i)
   {
      WindowTitle[i] = CharAlloc(256,"WindowTitle");
      if (wto == 0) sprintf(WindowTitle[i],"Window %d",i+1);
      else          sprintf(WindowTitle[i],"%sWindow %d",wtrun,i+1);
      WinAtt[i].x = -1;
      WinAtt[i].y = -1;
      WinAtt[i].w = -1;
      WinAtt[i].h = -1;
      WinAtt[i].active = 1;
   }

   CreateTestWindow();

   if (MRpid >= 0) sprintf(GUI_default,"GUI_%2.2d.cfg",MRpid);
   ReadConfig(GUI_default);
   if (MRpid >= 0) sprintf(GUI_config,"GUI_last_used_%2.2d.cfg",MRpid);
   ReadConfig(GUI_config);

   LoadFonts();
   CreateControlWindow();

   OutXSize = (ScreenW-ScreenOffset) / WinCols;
   OutYSize = (CowY-ScreenOffset-1-WM_top_area)  / WinRows;
   WinXSize = OutXSize - WinLM - WinRM;
   WinYSize = OutYSize - WinTM - WinBM;
   WinSizeHints.flags      = PPosition | PSize | PMinSize;
   WinSizeHints.min_width  = 200;
   WinSizeHints.min_height = 100;
   if (Debug)
   {
      printf("Outer windowsize = %dx%d\n",OutXSize,OutYSize);
      printf("Inner windowsize = %dx%d\n",WinXSize,WinYSize);
   }

   wm_hints.initial_state = NormalState;
   wm_hints.input = True;
   wm_hints.flags = StateHint | InputHint;
   
   class_hints.res_name = progname;
   class_hints.res_class = "PUMA";

   Delwin = XInternAtom(display,"WM_DELETE_WINDOW",0);

   ReadImage(&MapHR,"map.bmp");

   for (i = 0 ; i < NumWin ; ++i)
   {
      x = WinAtt[i].x;
      y = WinAtt[i].y;
      w = WinAtt[i].w;
      h = WinAtt[i].h;

      if (x < ScreenOffset || x >= ScreenW) x = ScreenOffset+(i%WinCols)*OutXSize;
      if (y < ScreenOffset || y >= ScreenH) y = ScreenOffset+(i/WinCols)*OutYSize+WM_top_area;
      if (w <  WinSizeHints.min_width ) w = WinXSize;
      if (h <  WinSizeHints.min_height) h = WinYSize;
      XStringListToTextProperty(&WindowTitle[i],1,WindowName+i);
      if (WinAtt[i].active)
      {
         if (MRnum == 2 && MRpid == 1) x += ScreenW;
         Win[i] = XCreateWindow(display,RootWindow(display,screen_num), // display, parent
                  x,y,w,h,
                  BorderWidth,CopyFromParent,InputOutput,         // border, depth, class
                  CopyFromParent,0,NULL);                         // visual, valuemask, attributes
         XSetWMProtocols(display,Win[i],&Delwin,1);
         XSetWMProperties(display,Win[i],WindowName+i,NULL,
   		NULL,0,&WinSizeHints,&wm_hints,&class_hints);
         XSelectInput(display,Win[i],ButtonPressMask | KeyPressMask | ExposureMask);
         XMapWindow(display,Win[i]);
      }
   }

   /* Prepare GC */

   gc = XCreateGC(display, Cow, valuemask, &values);
   XSetFont(display, gc, FixFont->fid);
   colormap = XDefaultColormap(display,screen_num);
   XSetForeground(display,gc,BlackPix);
   XSetBackground(display,gc,WhitePix);

   /* Get keyboard information */

   XDisplayKeycodes(display,&FirstKey,&LastKey);
   KeyMap = XGetKeyboardMapping(display,FirstKey,LastKey-FirstKey+1,&SymsPerKey);

   /* Allocate color cells */

   for (i=0 ; i < NUMPAL ; ++i)
      LineCo[i] = AllocateColorCells(Pallet[i]);
  
   /* Color cells for control window */

   XAllocNamedColor(display,colormap,"red"        ,&Red       ,&Dummy);
   XAllocNamedColor(display,colormap,"green"      ,&Green     ,&Dummy);
   XAllocNamedColor(display,colormap,"blue"       ,&Blue      ,&Dummy);
   XAllocNamedColor(display,colormap,"grey"       ,&Grey      ,&Dummy);
   XAllocNamedColor(display,colormap,"yellow"     ,&Yellow    ,&Dummy);
   XAllocNamedColor(display,colormap,"hot pink"   ,&LightRed  ,&Dummy);
   XAllocNamedColor(display,colormap,"dark red"   ,&DarkRed   ,&Dummy);
   XAllocNamedColor(display,colormap,"light blue" ,&LightBlue ,&Dummy);
   XAllocNamedColor(display,colormap,"dark blue"  ,&DarkBlue  ,&Dummy);
   XAllocNamedColor(display,colormap,"light green",&LightGreen,&Dummy);
   XAllocNamedColor(display,colormap,"dark green" ,&DarkGreen ,&Dummy);

   TSColor[0] = Red.pixel;
   TSColor[1] = Green.pixel;
   TSColor[2] = Blue.pixel;
   TSColor[3] = WhitePix;
   TSColor[4] = LightRed.pixel;
   TSColor[5] = Grey.pixel;
   TSColor[6] = Yellow.pixel;
   TSColor[7] = LightBlue.pixel;
   TSColor[8] = LightGreen.pixel;

   ClearTracer();
   CalcButtonAreas();
   RedrawControlWindow();
   XSync(display,0);
}

void FillPoly(int n, REAL Poly[])
{
   int i;
   XPoint xpol[8];
   for (i=0; i < n ; ++i)
   {
      xpol[i].x = OffX + Poly[i+i  ] + 0.5;
      xpol[i].y = OffY + Poly[i+i+1] + 0.5;
   }
   XFillPolygon(display,pix,gc,xpol,n,Convex,CoordModeOrigin);
}

/* ======= */
/* IsoArea */
/* ======= */

void IsoArea(INT y, REAL vl, REAL vh, REAL Top[], REAL Bot[], INT Dim)
{
   INT  f,x,p;
   REAL xl,xr,yt,yb;
   REAL Poly[16];

   // if (Debug) printf("IsoArea %d %10.2f %10.2f %10.2f %10.2f %d\n",y,vl,vh,Top[0],Bot[0],Dim);
   for (x=0 ; x < Dim-1 ; x++)
   {
      Flag[x] = 0;
      if (Top[x  ] <  vl) Flag[x] |= TOLELO;
      if (Top[x  ] >= vh) Flag[x] |= TOLEHI;
      if (Top[x+1] <  vl) Flag[x] |= TORILO;
      if (Top[x+1] >= vh) Flag[x] |= TORIHI;
      if (Bot[x  ] <  vl) Flag[x] |= BOLELO;
      if (Bot[x  ] >= vh) Flag[x] |= BOLEHI;
      if (Bot[x+1] <  vl) Flag[x] |= BORILO;
      if (Bot[x+1] >= vh) Flag[x] |= BORIHI;
   }

   x = 0;

   while (x < Dim-1)
   {
      xl = VGAX *  x   ;
      xr = VGAX * (x+1);
      yt = VGAY *  y   ;
      yb = VGAY * (y+1);
      f  = Flag[x];

      if (f == 0)
      {
	 x++;
	 while (x < Dim-1 && Flag[x] == 0)
	 {
	    x++;
	    xr = VGAX * x;
	 }
	 Poly[0] = Poly[6] = xl;
	 Poly[1] = Poly[3] = yt;
	 Poly[2] = Poly[4] = xr;
	 Poly[5] = Poly[7] = yb;

	 FillPoly(4,Poly);

      }
      else if (f == (TOLELO | TORILO | BOLELO | BORILO) ||
	       f == (TOLEHI | TORIHI | BOLEHI | BORIHI))   x++;
      else if (Top[x] < vl && Top[x+1] >= vl && Bot[x] >= vl && Bot[x+1] < vl)
      {
	 Poly[1] = Poly[3] = yt;
	 Poly[2] = Poly[4] = xr;
         Poly[0] = IPX(Top[x  ],vl,Top[x+1]);
         Poly[5] = IPY(Top[x+1],vl,Bot[x+1]);
         FillPoly(3,Poly);
         Poly[1] = Poly[3] = yb;
         Poly[2] = Poly[4] = xl;
         Poly[0] = IPX(Bot[x  ],vl,Bot[x+1]);
         Poly[5] = IPY(Top[x  ],vl,Bot[x  ]);
         FillPoly(3,Poly);
         ++x;
      }
      else if (Top[x] >= vl && Top[x+1] < vl && Bot[x] < vl && Bot[x+1] >= vl)
      {
	 Poly[0] = Poly[4] = xl;
	 Poly[1] = Poly[3] = yt;
         Poly[2] = IPX(Top[x  ],vl,Top[x+1]);
         Poly[5] = IPY(Top[x  ],vl,Bot[x  ]);
         FillPoly(3,Poly);
         Poly[0] = Poly[2] = xr;
         Poly[3] = Poly[5] = yb;
         Poly[1] = IPY(Top[x+1],vl,Bot[x+1]);
         Poly[4] = IPX(Bot[x  ],vl,Bot[x+1]);
         FillPoly(3,Poly);
         ++x;
      }
      else
      {
	 p = 0;

	 if (Top[x] < vl)
	 {
	    if (Top[x+1] >= vl)
	    {
	       Poly[p++] = IPX(Top[x],vl,Top[x+1]);
	       Poly[p++] = yt;
	    }

	    if (Top[x+1] >= vh)
	    {
	       Poly[p++] = IPX(Top[x],vh,Top[x+1]);
	       Poly[p++] = yt;
	    }
	 }
	 else if (Top[x] >= vh)
	 {
	    if (Top[x+1] < vh)
	    {
	       Poly[p++] = IPX(Top[x],vh,Top[x+1]);
	       Poly[p++] = yt;
	    }
	    if (Top[x+1] < vl)
	    {
	       Poly[p++] = IPX(Top[x],vl,Top[x+1]);
	       Poly[p++] = yt;
	    }
	 }
	 else
	 {
	    Poly[p++] = xl;
	    Poly[p++] = yt;
	    if (Top[x+1] < vl)
	    {
	       Poly[p++] = IPX(Top[x],vl,Top[x+1]);
	       Poly[p++] = yt;
	    }
	    if (Top[x+1] >= vh)
	    {
	       Poly[p++] = IPX(Top[x],vh,Top[x+1]);
	       Poly[p++] = yt;
	    }
	 }

	 if (Top[x+1] < vl)
	 {
	    if (Bot[x+1] >= vl)
	    {
	       Poly[p++] = xr;
	       Poly[p++] = IPY(Top[x+1],vl,Bot[x+1]);
	    }

	    if (Bot[x+1] >= vh)
	    {
	       Poly[p++] = xr;
	       Poly[p++] = IPY(Top[x+1],vh,Bot[x+1]);
	    }
	 }
	 else if (Top[x+1] >= vh)
	 {
	    if (Bot[x+1] < vh)
	    {
	       Poly[p++] = xr;
	       Poly[p++] = IPY(Top[x+1],vh,Bot[x+1]);
	    }
	    if (Bot[x+1] < vl)
	    {
	       Poly[p++] = xr;
	       Poly[p++] = IPY(Top[x+1],vl,Bot[x+1]);
	    }
	 }
	 else
	 {
	    Poly[p++] = xr;
	    Poly[p++] = yt;
	    if (Bot[x+1] < vl)
	    {
	       Poly[p++] = xr;
	       Poly[p++] = IPY(Top[x+1],vl,Bot[x+1]);
	    }
	    if (Bot[x+1] >= vh)
	    {
	       Poly[p++] = xr;
	       Poly[p++] = IPY(Top[x+1],vh,Bot[x+1]);
	    }
	 }

	 if (Bot[x+1] < vl)
	 {
	    if (Bot[x] >= vl)
	    {
	       Poly[p++] = IPX(Bot[x],vl,Bot[x+1]);
	       Poly[p++] = yb;
	    }
	    if (Bot[x] >= vh)
	    {
	       Poly[p++] = IPX(Bot[x],vh,Bot[x+1]);
	       Poly[p++] = yb;
	    }
	 }
	 else if (Bot[x+1] >= vh)
	 {
	    if (Bot[x] < vh)
	    {
	       Poly[p++] = IPX(Bot[x],vh,Bot[x+1]);
	       Poly[p++] = yb;
	    }
	    if (Bot[x] < vl)
	    {
	       Poly[p++] = IPX(Bot[x],vl,Bot[x+1]);
	       Poly[p++] = yb;
	    }
	 }
	 else
	 {
	    Poly[p++] = xr;
	    Poly[p++] = yb;
	    if (Bot[x] < vl)
	    {
	       Poly[p++] = IPX(Bot[x],vl,Bot[x+1]);
	       Poly[p++] = yb;
	    }
	    if (Bot[x] >= vh)
	    {
	       Poly[p++] = IPX(Bot[x],vh,Bot[x+1]);
	       Poly[p++] = yb;
	    }
	 }

	 if (Bot[x] < vl)
	 {
	    if (Top[x] >= vl)
	    {
	       Poly[p++] = xl;
	       Poly[p++] = IPY(Top[x],vl,Bot[x]);
	    }

	    if (Top[x] >= vh)
	    {
	       Poly[p++] = xl;
	       Poly[p++] = IPY(Top[x],vh,Bot[x]);
	    }
	 }
	 else if (Bot[x] >= vh)
	 {
	    if (Top[x] < vh)
	    {
	       Poly[p++] = xl;
	       Poly[p++] = IPY(Top[x],vh,Bot[x]);
	    }
	    if (Top[x] < vl)
	    {
	       Poly[p++] = xl;
	       Poly[p++] = IPY(Top[x],vl,Bot[x]);
	    }
	 }
	 else
	 {
	    Poly[p++] = xl;
	    Poly[p++] = yb;
	    if (Top[x] < vl)
	    {
	       Poly[p++] = xl;
	       Poly[p++] = IPY(Top[x],vl,Bot[x]);
	    }
	    if (Top[x] >= vh)
	    {
	       Poly[p++] = xl;
	       Poly[p++] = IPY(Top[x],vh,Bot[x]);
	    }
	 }
	 FillPoly(p>>1,Poly);
	 x++;
      }
   }
}


/* ======= */
/* isoarea */
/* ======= */

void IsoAreas(struct ColorStrip Strip[])
{
   INT  i;
   INT  y;
   REAL *Top;
   REAL *Bot;

   i = 0;
   while (Strip[i].Name)
   {
      Top = Field;
      Bot = Field + DimX;
      XSetForeground(display,gc,Strip[i].pixel);
      for (y = 0 ; y < DimY-1 ; y++)
      {
	 IsoArea(y,Strip[i].Lo,Strip[i].Hi,Top,Bot,DimX);
	 Top += DimX;
	 Bot += DimX;
      }
      ++i;
   }
}

/* ======= */
/* IsoLine */
/* ======= */

void IsoLine(INT y, REAL v, REAL Top[], REAL Bot[], INT Dim)
{
   INT  f,x;
   REAL xl,xr,yt,yb;
   REAL x1,yo,x2,y2,x3,y3,x4,y4;

   for (x=0 ; x < Dim-1 ; x++)
   {
      Flag[x] = 0;
      if (Top[x  ] <  v) Flag[x] |= TOLELO;
      else               Flag[x] |= TOLEHI;
      if (Top[x+1] <  v) Flag[x] |= TORILO;
      else               Flag[x] |= TORIHI;
      if (Bot[x  ] <  v) Flag[x] |= BOLELO;
      else               Flag[x] |= BOLEHI;
      if (Bot[x+1] <  v) Flag[x] |= BORILO;
      else               Flag[x] |= BORIHI;
   }

   x = 0;

   while (x < Dim-1)
   {
      xl = VGAX *  x   ;
      xr = VGAX * (x+1);
      yt = VGAY *  y   ;
      yb = VGAY * (y+1);
      f  = Flag[x];

      if (f == 0 || f == (TOLELO | TORILO | BOLELO | BORILO) ||
		    f == (TOLEHI | TORIHI | BOLEHI | BORIHI))   x++;
      else if (f == (TOLELO | BORILO | TORIHI | BOLEHI))
      {
         x1 = IPX(Top[x  ],v,Top[x+1]);
         yo = yt;
         x2 = xr;
         y2 = IPY(Top[x+1],v,Bot[x+1]);
         x3 = IPX(Bot[x  ],v,Bot[x+1]);
         y3 = yb;
         x4 = xl;
         y4 = IPY(Top[x  ],v,Bot[x  ]);
	 XDrawLine(display,pix,gc,x1,OffY+yo,x2,OffY+y2);
	 XDrawLine(display,pix,gc,x3,OffY+y3,x4,OffY+y4);
         ++x;
      }
      else if (f == (TORILO | BOLELO | TOLEHI | BORIHI))
      {
         x1 = xl;
         yo = IPY(Top[x  ],v,Bot[x  ]);
         x2 = IPX(Top[x  ],v,Top[x+1]);
         y2 = yt;
         x3 = xr;
         y3 = IPY(Top[x+1],v,Bot[x+1]);
         x4 = IPX(Bot[x  ],v,Bot[x+1]);
         y4 = yb;
	 XDrawLine(display,pix,gc,x1,OffY+yo,x2,OffY+y2);
	 XDrawLine(display,pix,gc,x3,OffY+y3,x4,OffY+y4);
         ++x;
      }
      else
      {
	 x1 = x2 = x3 = x4 = -1;

	 if ((Top[x  ] < v && Top[x+1] >= v) || (Top[x  ] >= v && Top[x+1] < v))
	 {
	    x1 = IPX(Top[x  ],v,Top[x+1]);
	    yo = yt;
	 }
	 if ((Top[x+1] < v && Bot[x+1] >= v) || (Top[x+1] >= v && Bot[x+1] < v))
	 {
	    x2 = xr;
	    y2 = IPY(Top[x+1],v,Bot[x+1]);
	 }
	 if ((Bot[x+1] < v && Bot[x  ] >= v) || (Bot[x+1] >= v && Bot[x  ] < v))
	 {
	    x3 = IPX(Bot[x  ],v,Bot[x+1]);
	    y3 = yb;
	 }
	 if ((Bot[x  ] < v && Top[x  ] >= v) || (Bot[x  ] >= v && Top[x  ] < v))
	 {
	    x4 = xl;
	    y4 = IPY(Top[x  ],v,Bot[x  ]);
	 }

	 if (x1 >= 0 && x2 >= 0 && x3 >= 0 && x4 >= 0)
	 {
	    XDrawLine(display,pix,gc,x1,OffY+yo,x2,OffY+y2);
	    XDrawLine(display,pix,gc,x3,OffY+y3,x4,OffY+y4);
	 }
	 else if (x1 >= 0)
	 {
	    if      (x2 >= 0) XDrawLine(display,pix,gc,x1,OffY+yo,x2,OffY+y2);
	    else if (x3 >= 0) XDrawLine(display,pix,gc,x1,OffY+yo,x3,OffY+y3);
	    else              XDrawLine(display,pix,gc,x1,OffY+yo,x4,OffY+y4);
	 }
	 else if (x2 >= 0)
	 {
	    if (x3 >= 0) XDrawLine(display,pix,gc,x2,OffY+y2,x3,OffY+y3);
	    else         XDrawLine(display,pix,gc,x2,OffY+y2,x4,OffY+y4);
	 }
	 else XDrawLine(display,pix,gc,x3,OffY+y3,x4,OffY+y4);

	 x++;
      }
   }
}


/* ======== */
/* IsoLines */
/* ======== */

void IsoLines(struct ColorStrip Strip[],int Colored)
{
   INT  i;
   INT  y;
   REAL *Top;
   REAL *Bot;

   XSetForeground(display,gc,BlackPix);
   i = 0;
   while (Strip[i].Name)
   {
      Top = Field;
      Bot = Field + DimX;
      if (Colored) XSetForeground(display,gc,Strip[i].pixel);
      for (y = 0 ; y < DimY-1 ; y++)
      {
	 IsoLine(y,Strip[i].Lo,Top,Bot,DimX);
	 Top += DimX;
	 Bot += DimX;
      }
      ++i;
   }
}

int AziPoint(int fpx, int fpy, int *x2, int *y2)
{
   int xc,yc;
   double flam,fsin,fcos,roa,pil,piy;

   xc  = InXSize >> 1;
   yc  = InYSize >> 1;
   roa = MapLR[win].l * InXSize / 360.0;
   pil = 2.0 * M_PI / InXSize;
   piy = M_PI / InYSize;

   flam = (fpx - xc - roa) * pil;
   if (flam < -M_PI) flam += 2.0 * M_PI;
   if (flam >  M_PI) flam -= 2.0 * M_PI;
   if (flam < -M_PI_2 || flam > M_PI_2) return 1;
   fsin = sin(flam);
   fcos = cos((fpy-yc) * piy);
   *x2  = xc + InYSize/2 * fcos * fsin;
   *y2  = yc + InYSize/2 * sin((fpy-yc) * piy);
   return 0;
}

void AziLine(int x1, int y1, int x2, int y2)
{
   int a1,b1,a2,b2;

   if (AziPoint(x1,y1,&a1,&b1)) return;
   if (AziPoint(x2,y2,&a2,&b2)) return;
   XDrawLine(display,pix,gc,a1,b1,a2,b2);
}


/* ======= */
/* MapLine */
/* ======= */

void MapLine(INT y, REAL v, REAL Top[], REAL Bot[], INT Dim)
{
   INT  f,x;
   REAL xl,xr,yt,yb;
   REAL x1,yo,x2,y2,x3,y3,x4,y4;

   for (x=0 ; x < Dim-1 ; x++)
   {
      Flag[x] = 0;
      if (Top[x  ] <  v) Flag[x] |= TOLELO;
      else               Flag[x] |= TOLEHI;
      if (Top[x+1] <  v) Flag[x] |= TORILO;
      else               Flag[x] |= TORIHI;
      if (Bot[x  ] <  v) Flag[x] |= BOLELO;
      else               Flag[x] |= BOLEHI;
      if (Bot[x+1] <  v) Flag[x] |= BORILO;
      else               Flag[x] |= BORIHI;
   }

   x = 0;

   while (x < Dim-1)
   {
      xl = VGAX *  x   ;
      xr = VGAX * (x+1);
      yt = VGAY *  y   ;
      yb = VGAY * (y+1);
      f  = Flag[x];

      if (f == 0 || f == (TOLELO | TORILO | BOLELO | BORILO) ||
		    f == (TOLEHI | TORIHI | BOLEHI | BORIHI))   x++;
      else if (f == (TOLELO | BORILO | TORIHI | BOLEHI))
      {
         x1 = IPX(Top[x  ],v,Top[x+1]);
         yo = yt;
         x2 = xr;
         y2 = IPY(Top[x+1],v,Bot[x+1]);
         x3 = IPX(Bot[x  ],v,Bot[x+1]);
         y3 = yb;
         x4 = xl;
         y4 = IPY(Top[x  ],v,Bot[x  ]);
	 AziLine(x1,OffY+yo,x2,OffY+y2);
	 AziLine(x3,OffY+y3,x4,OffY+y4);
         ++x;
      }
      else if (f == (TORILO | BOLELO | TOLEHI | BORIHI))
      {
         x1 = xl;
         yo = IPY(Top[x  ],v,Bot[x  ]);
         x2 = IPX(Top[x  ],v,Top[x+1]);
         y2 = yt;
         x3 = xr;
         y3 = IPY(Top[x+1],v,Bot[x+1]);
         x4 = IPX(Bot[x  ],v,Bot[x+1]);
         y4 = yb;
	 AziLine(x1,OffY+yo,x2,OffY+y2);
	 AziLine(x3,OffY+y3,x4,OffY+y4);
         ++x;
      }
      else
      {
	 x1 = x2 = x3 = x4 = -1;

	 if ((Top[x  ] < v && Top[x+1] >= v) || (Top[x  ] >= v && Top[x+1] < v))
	 {
	    x1 = IPX(Top[x  ],v,Top[x+1]);
	    yo = yt;
	 }
	 if ((Top[x+1] < v && Bot[x+1] >= v) || (Top[x+1] >= v && Bot[x+1] < v))
	 {
	    x2 = xr;
	    y2 = IPY(Top[x+1],v,Bot[x+1]);
	 }
	 if ((Bot[x+1] < v && Bot[x  ] >= v) || (Bot[x+1] >= v && Bot[x  ] < v))
	 {
	    x3 = IPX(Bot[x  ],v,Bot[x+1]);
	    y3 = yb;
	 }
	 if ((Bot[x  ] < v && Top[x  ] >= v) || (Bot[x  ] >= v && Top[x  ] < v))
	 {
	    x4 = xl;
	    y4 = IPY(Top[x  ],v,Bot[x  ]);
	 }

	 if (x1 >= 0 && x2 >= 0 && x3 >= 0 && x4 >= 0)
	 {
	    AziLine(x1,OffY+yo,x2,OffY+y2);
	    AziLine(x3,OffY+y3,x4,OffY+y4);
	 }
	 else if (x1 >= 0)
	 {
	    if      (x2 >= 0) AziLine(x1,OffY+yo,x2,OffY+y2);
	    else if (x3 >= 0) AziLine(x1,OffY+yo,x3,OffY+y3);
	    else              AziLine(x1,OffY+yo,x4,OffY+y4);
	 }
	 else if (x2 >= 0)
	 {
	    if (x3 >= 0) AziLine(x2,OffY+y2,x3,OffY+y3);
	    else         AziLine(x2,OffY+y2,x4,OffY+y4);
	 }
	 else AziLine(x3,OffY+y3,x4,OffY+y4);

	 x++;
      }
   }
}


/* ======== */
/* MapLines */
/* ======== */

void MapLines(struct ColorStrip Strip[],int Colored)
{
   INT  i;
   INT  y;
   REAL *Top;
   REAL *Bot;

   XSetForeground(display,gc,BlackPix);
   i = 0;
   while (Strip[i].Name)
   {
      Top = Field;
      Bot = Field + DimX;
      if (Colored) XSetForeground(display,gc,Strip[i].pixel);
      for (y = 0 ; y < DimY-1 ; y++)
      {
	 MapLine(y,Strip[i].Lo,Top,Bot,DimX);
	 Top += DimX;
	 Bot += DimX;
      }
      ++i;
   }
}

int polco = 0;

void poltra(float lam, float phi, float *x, float *y)
{
   REAL xnp,ynp,xsp,ysp,lfa;

   lfa = (2.0 * M_PI) / (InXSize-1) ;
   xnp = (InXSize-1) * 0.25;
   ynp = (InYSize-1) * 0.50;
   xsp = (InXSize-1) * 0.75;
   ysp = ynp;

   /* Northern hemisphere */

   if (phi < ynp)
   {
      *x = xnp + sin(lam * lfa) * phi;
      *y = ynp + cos(lam * lfa) * phi;
   }

   /* Souhern hemisphere */

   else
   {
      *x = xsp - sin(lam * lfa) * (phi-ysp);
      *y = ysp + cos(lam * lfa) * (phi-ysp);
   }
}


/* ========== */
/* TracerPlot */
/* ========== */

void TracerPlot(int w)
{
   int j,ipx,ipy,lpy,lpx,jx,jy,ic,ix,iy,ipr;
   int xnp,ynp,xsp,ysp,xc,yc;
   float *u;
   float *v;
   float fpx,fpy,du,dv,ppx,ppy,rx,ry,spfac;
   float lfa,lam,phi,hpy,fcos,fsin;
   float pil,piy;
   float roa;
   float flam;
   float dpx,dpy;

   if (Uindex < 0 || Vindex < 0) return; // Data missing

   Delpar = (8 * 60) / DeltaTime; // Injection interval
   if (LevelChanged[w]) ClearTracer();
   LevelChanged[w] = 0;
   if (!SpeedScale) SpeedScale = FloatAlloc(DimY,"SpeedScale");
   u = Array[Uindex].Data + Indez[w] * DimXY;
   v = Array[Vindex].Data + Indez[w] * DimXY;

   TracerColor = Yellow.pixel;
   if (nstep % (8*Delpar) == 0) ColInd++;
   if (ColInd >= AUTOCOLORS) ColInd = 0;
   TracerColor = Autostrip[ColInd].pixel;

   dpx = (float)(InXSize-1) / (float)(DimX-1);
   dpy = (float)(InYSize-1) / (float)(DimY-1);
   lpx = InXSize;
   lpy = InYSize;
   hpy = 0.5 * dpy;
   xc  = lpx / 2;
   yc  = lpy / 2;
   ipr = InXSize / 200;
   pil = 2.0 * M_PI / lpx;
   piy = M_PI / lpy;
   roa = MapLR[w].l * lpx / 360.0;

   // Compute scale factors

   if (SizeChanged || nstep < 10)
   {
      for (j=0 ; j < DimY; ++j)
      {
         SpeedScale[j] = InXSize * DeltaTime * 60.0 / 40000000.0 
                       / cos((j-0.5*(DimY-1)) * (M_PI/DimY));
      }
      SpeedScale[DimY-1] = SpeedScale[0] = SpeedScale[1];
   }

   // Insert new particles

   if (nstep % Delpar == 0)
   {
      fpx = dpx * (DimX-1) / 2.0;
      for (ipy = 0 ; ipy < lpy ; ipy += dpy)
      {
         ParInd = (ParInd+1) & (MAXPAR-1);
         pax[ParInd] = fpx;
         pay[ParInd] = ipy + hpy;
         pal[ParInd] = TracerColor;
      }
   }

  // Move particles

  for (j=0 ; j < MAXPAR ; ++j)
  if (pax[j] >= 0.0 && pay[j] >= 0.0)
  {
     jx = rx = pax[j] / dpx;
     jy = ry = pay[j] / dpy;
     ic = jx + DimX * jy;
     ix = ic + 1;
     if (jx >= DimX-1) ix -= DimX;
     iy = ic + DimX;
     if (iy >= DimX * DimY) iy -= DimX;
     rx -= jx;
     ry -= jy;

     du = u[ic]+(u[ix]-u[ic])*rx+(u[iy]-u[ic])*ry;
     dv = v[ic]+(v[ix]-v[ic])*rx+(v[iy]-v[ic])*ry;
     spfac = SpeedScale[jy] + (SpeedScale[jy+1]-SpeedScale[jy]) * ry;
     du *= spfac;
     dv *= spfac;
     fpx = pax[j] + du;
     fpy = pay[j] + dv;
     if (fpx <  0.0) fpx += lpx;
     if (fpx >= lpx) fpx -= lpx;
     if (fpy < 0.0 || fpy > lpy-1)
     {
        pay[j] = -1.0;
        continue;
     }
     pax[j] = fpx;
     pay[j] = fpy;
     ipx    =  -1;
     if (MapPro[w] == POLAR)
     {
       lfa = (2.0 * M_PI) / lpx;
       xnp = (InXSize-1) * 0.25;
       ynp = (InYSize-1) * 0.50;
       xsp = (InXSize-1) * 0.75;
       ysp = ynp;
       lam = fpx;
       phi = fpy;
       if (phi < 0.5 * lpy) /* Northern hemisphere */
       {
          fpx = xnp - sin(lam * lfa) * (phi + hpy);
          fpy = ynp - cos(lam * lfa) * (phi + hpy);
       }
       else                 /* Souhern hemisphere */
       {
          fpx = xsp + sin(lam * lfa) * (lpy - phi + hpy);
          fpy = ysp - cos(lam * lfa) * (lpy - phi + hpy);
       }
       ipx = rintf(fpx);
       ipy = rintf(fpy);
     }
     else if (MapPro[w] == AZIMUTHAL)
     {
        flam = (fpx - xc - roa) * pil;
        if (flam < -M_PI) flam += 2.0 * M_PI;
        if (flam >  M_PI) flam -= 2.0 * M_PI;
        if (flam < -M_PI_2 || flam > M_PI_2) continue;
        fsin = sin(flam);
        fcos = cos((fpy-yc) * piy);
        fpx  = xc + lpy/2 * fcos * fsin;
        fpy  = yc + lpy/2 * sin((fpy-yc) * piy);
        ipx  = rintf(fpx);
        ipy  = rintf(fpy);
     }
     else
     {
        ipx = rintf(fpx);
        ipy = rintf(fpy);
     }
     if (ipx >= 0)
     {
        XSetForeground(display,gc,pal[j]);
        if (ipr < 2) XDrawPoint(display,pix,gc,ipx,ipy);
        else if (ipr == 2) XFillRectangle(display,pix,gc,ipx,ipy,2,2);
        else XFillArc(display,pix,gc,ipx,ipy,ipr,ipr,0,360*64);
     }
  }
}


/* ============= */
/* AmplitudePlot */
/* ============= */

void AmplitudePlot(void)
{
   int i,j,m,n;
   int x,y;
   int r,dx,dy;
   int mmax;
   int imax;
   REAL Amax,Fac;

   if (!Ampli) Ampli = FloatAlloc(DimX,"Ampli");
   if (!Acol ) Acol  = SizeAlloc(DimX,sizeof(XColor),"Acol");

   Amax = 0.0;
   imax = 0;
   for (i=0 ; i < DimX ; ++i)
   {
      Ampli[i] = log (1.0 + sqrt(Field[2*i]*Field[2*i]+Field[2*i+1]*Field[2*i+1]));
      Acol[i].pixel = WhitePix;
      if (Ampli[i] > Amax && i > DimY)
      {
         Amax = Ampli[i];
         imax = i;
      }
   }

   Fac = 0.0;
   if (Amax > 1.e-20) Fac = 1.0 / Amax;
   for (i=0 ; i < DimX ; ++i)
   {
      j = AMPLI_COLS * Ampli[i] * Fac;
      if (j >= AMPLI_COLS) j = AMPLI_COLS-1;
      Acol[i].pixel = AmpliStrip[j].pixel;
   }

   dx = WinXSize / 23;
   r  = (dx+dx) / 3;
   mmax = (InYSize / dx) - 1;
   XSetForeground(display,gc,BlackPix);
   for (m=0,i=0 ; m < mmax ; ++m)
   {
      y = FixFontHeight + 4 + m * dx;
      for (n=m ; n < DimY ; ++n,++i)
      if (n < 22)
      {
         x = dx/2 + n * dx;
         XSetForeground(display,gc,Acol[i].pixel);
         XFillArc(display,pix,gc,x,y,r,r,0,360*64);
      }
   }
}


void RedrawIsoWindow(int w)
{
   if (WinPixMap[w].Pix)
   XCopyArea(display,WinPixMap[w].Pix,Win[w],gc,0,0,WinPixMap[w].DimX,WinPixMap[w].DimY,0,0);
}


void CloseWindow(int i)
{
   INTXS X,Y,xp,yp;
   unsigned int border,depth;
   Window Rootwin,Child;

   if (Win[i])
   {
       XGetGeometry(display,Win[i],&Rootwin,&xp,&yp,&WinAtt[i].w,&WinAtt[i].h,&border,&depth);
       XTranslateCoordinates(display,Win[i],Rootwin,0,0,&X,&Y,&Child);
       WinAtt[i].x =  X - WinLM;
       WinAtt[i].y =  Y - WinTM;
       XDestroyWindow(display,Win[i]);
       Win[i] = 0;
       ShowWindowStatus();
   }
}

void ReopenWindow(int i)
{
   if (Win[i] == 0)
   {
      if (WinAtt[i].x < ScreenOffset) WinAtt[i].x = ScreenOffset + OutXSize * (i % WinCols);
      if (WinAtt[i].y < ScreenOffset) WinAtt[i].y = ScreenOffset + OutYSize * (i / WinCols);
      Win[i] = XCreateSimpleWindow(display,RootWindow(display,screen_num), 
               WinAtt[i].x,WinAtt[i].y,WinAtt[i].w,WinAtt[i].h,
               BorderWidth,BlackPix,WhitePix);
      XSetWMProtocols(display,Win[i],&Delwin,1);
      XSetWMProperties(display,Win[i],WindowName+i,NULL,
		NULL,0,&WinSizeHints,&wm_hints,&class_hints);
      XSelectInput(display,Win[i],ButtonPressMask | KeyPressMask | ExposureMask);
      XMapWindow(display,Win[i]);
      ShowWindowStatus();
   }
}

char *HelpText[] =
{
   "Key assignment for azimuthal projection     ",
   "--------------------------------------------",
   "<- cursor left  : increase westward rotation",
   "-> cursor right : increase eastward rotation",
   "============================================",
   "Key assignment for 3D variables             ",
   "--------------------------------------------",
   "^  cursor up    : display next upper level  ",
   "v  cursor down  : display next lower level  ",
   "============================================",
   "Key assignment for Lon/Lev or Hovmoeller    ",
   "--------------------------------------------",
   "^  cursor up    : switch latitude northwards",
   "v  cursor down  : switch latitude southwards",
   "============================================",
   "Mouse button assignment                     ",
   "--------------------------------------------",
   "Left   button   : decrease level or latitude",
   "Right  button   : increase level or latitude",
   "Middle button   : cycle through projections "
};

int HelpLines = sizeof(HelpText) / sizeof(char *);

void DrawHelpWindow(void)
{
   int x,y,j,l;

   XSetForeground(display,gc,Yellow.pixel);
   XSetBackground(display,gc,DarkBlue.pixel);
   XSetFont(display,gc,BigFont->fid);

   x = BigFontWidth / 2;
   for (j=0 ; j < HelpLines ; ++j)
   {
      y = BigFontHeight / 2 + BigFont->ascent + j * BigFontHeight;
      l = strlen(HelpText[j]);
      XDrawImageString(display,HelpWindow,gc,x,y,HelpText[j],l);
   }
}

void DisplayHelpWindow(void)
{
   int x,y,w,h;
   char *HelpTitle = "Help";

   w = BigFontWidth  * (strlen(HelpText[0])+1);
   h = BigFontHeight * (HelpLines + 1);
   x = (ScreenW - w) / 2;
   y = (ScreenH - h) / 2;
   HelpWindow = XCreateSimpleWindow(display,RootWindow(display,screen_num), 
                x,y,w,h,BorderWidth,Yellow.pixel,DarkBlue.pixel);
   XStringListToTextProperty(&HelpTitle,1,&HelpName);
   XSetWMProtocols(display,HelpWindow,&Delwin,1);
   XSetWMProperties(display,HelpWindow,&HelpName,NULL,NULL,0,NULL,NULL,NULL);
   XSelectInput(display,HelpWindow,ButtonPressMask | ExposureMask);
   XMapWindow(display,HelpWindow);
   XMoveWindow(display,HelpWindow,x,y);
   DrawHelpWindow();
}

int HitGridButton(XEvent* xe)
{
   return (
   (xe->xbutton.button == Button1) &&
   (xe->xbutton.x >= Grid_XL  ) &&
   (xe->xbutton.x <= Grid_XH  ) &&
   (xe->xbutton.y >= Grid_YL  ) &&
   (xe->xbutton.y <= Grid_YH  ));
}

int HitPauseButton(XEvent* xe)
{
   return (
   (xe->xbutton.button == Button1) &&
   (xe->xbutton.x >= Pause_XL-2  ) &&
   (xe->xbutton.x <= Pause_XH+2  ) &&
   (xe->xbutton.y >= Pause_YL-2  ) &&
   (xe->xbutton.y <= Pause_YH+2  ));
}

int HitStartButton(XEvent* xe)
{
   return (
   (xe->xbutton.button == Button1) &&
   (xe->xbutton.x >= Start_XL-2  ) &&
   (xe->xbutton.x <= Start_XH+2  ) &&
   (xe->xbutton.y >= Start_YL-2  ) &&
   (xe->xbutton.y <= Start_YH+2  ));
}

int HitStopButton(XEvent* xe)
{
   return (
   (xe->xbutton.button == Button1) &&
   (xe->xbutton.x >= Stop_XL-2   ) &&
   (xe->xbutton.x <= Stop_XH+2   ) &&
   (xe->xbutton.y >= Stop_YL-2   ) &&
   (xe->xbutton.y <= Stop_YH+2   ));
}

int HitHelpButton(XEvent* xe)
{
   return (
   (xe->xbutton.button == Button1) &&
   (xe->xbutton.x >= Help_XL-2   ) &&
   (xe->xbutton.x <= Help_XH+2   ) &&
   (xe->xbutton.y >= Help_YL-2   ) &&
   (xe->xbutton.y <= Help_YH+2   ));
}

int HitPlusButton(XEvent* xe)
{
   return (
   (xe->xbutton.button == Button1) &&
   (xe->xbutton.x >= Plus_XL  ) &&
   (xe->xbutton.x <= Plus_XH  ) &&
   (xe->xbutton.y >= Plus_YL  ) &&
   (xe->xbutton.y <= Plus_YH  ));
}

int HitFFWDButton(XEvent* xe)
{
   return (
   (xe->xbutton.button == Button1) &&
   (xe->xbutton.x >= FFWD_XL  ) &&
   (xe->xbutton.x <= FFWD_XH  ) &&
   (xe->xbutton.y >= FFWD_YL  ) &&
   (xe->xbutton.y <= FFWD_YH  ));
}

int HitMinusButton(XEvent* xe)
{
   return (
   (xe->xbutton.button == Button1) &&
   (xe->xbutton.x >= Minus_XL  ) &&
   (xe->xbutton.x <= Minus_XH  ) &&
   (xe->xbutton.y >= Minus_YL  ) &&
   (xe->xbutton.y <= Minus_YH  ));
}

int HitFBWDButton(XEvent* xe)
{
   return (
   (xe->xbutton.button == Button1) &&
   (xe->xbutton.x >= FBWD_XL  ) &&
   (xe->xbutton.x <= FBWD_XH  ) &&
   (xe->xbutton.y >= FBWD_YL  ) &&
   (xe->xbutton.y <= FBWD_YH  ));
}

int HitUpperPanel(XEvent* xe)
{
   return (
   (xe->xbutton.button == Button1) &&
   (xe->xbutton.x >= Parc_XL  ) &&
   (xe->xbutton.x <= Parc_XH  ) &&
   (xe->xbutton.y >= Parc_YL-Parc_YD  ) &&
   (xe->xbutton.y <= Parc_YL  ));
}

int HitLowerPanel(XEvent* xe)
{
   return (
   (xe->xbutton.button == Button1) &&
   (xe->xbutton.x >= Parc_XL  ) &&
   (xe->xbutton.x <= Parc_XH  ) &&
   (xe->xbutton.y >= Parc_YH  ) &&
   (xe->xbutton.y <= Parc_YH+Parc_YD));
}

int HitWindowSelect(XEvent* xe)
{
   if (
   (xe->xbutton.button == Button1) &&
   (xe->xbutton.x >= Wbox_XL) &&
   (xe->xbutton.x <= Wbox_XH) &&
   (xe->xbutton.y >= Wbox_YL) &&
   (xe->xbutton.y <= Wbox_YH))
      return (xe->xbutton.y / FixFontHeight); /* Window number */
   else return -1;
}

void OnMouseClick(void)
{
   int w;

   if (HitStopButton(&CowEvent))
   {
      Shutdown = 1;
   }
   else if (HitStartButton(&CowEvent))
   {
      Paused = 0;
      XSetWMName(display,Cow,&WinconName1);
   }
   else if (HitGridButton(&CowEvent))
   {
      Grid = !Grid;
      ShowGridStatus();
      for (w=0 ; w < NUMWIN ; ++w) RedrawFlag[w] = 1;
   }
   else if (HitPauseButton(&CowEvent))
   {
      if (Paused)
      {
         Paused = 0;
         XSetWMName(display,Cow,&WinconName1);
      }
      else
      {
         Paused = 1;
         XSetWMName(display,Cow,&WinconPause);
      }
   }
   else if (HitHelpButton(&CowEvent))
   {
      if (!HelpWindow) DisplayHelpWindow();
   }
   else if ((w = HitWindowSelect(&CowEvent)) >= 0)
   {
      if (Win[w]) CloseWindow(w);
      else        ReopenWindow(w);
   }
   else if (HitPlusButton(&CowEvent))
   {
      Parc[cpi].Val += Parc[cpi].Inc;
      if (Parc[cpi].Val > Parc[cpi].Max) Parc[cpi].Val = Parc[cpi].Max;
   }
   else if (HitFFWDButton(&CowEvent))
   {
      Parc[cpi].Val += 10.0 * Parc[cpi].Inc;
      if (Parc[cpi].Val > Parc[cpi].Max) Parc[cpi].Val = Parc[cpi].Max;
   }
   else if (HitMinusButton(&CowEvent))
   {
       Parc[cpi].Val -= Parc[cpi].Inc;
      if (Parc[cpi].Val < Parc[cpi].Min) Parc[cpi].Val = Parc[cpi].Min;
   }
   else if (HitFBWDButton(&CowEvent))
   {
       Parc[cpi].Val -= 10.0 * Parc[cpi].Inc;
      if (Parc[cpi].Val < Parc[cpi].Min) Parc[cpi].Val = Parc[cpi].Min;
   }
   else if (HitUpperPanel(&CowEvent) && cpi > 0) --cpi;
   else if (HitLowerPanel(&CowEvent) && cpi < Parcs-1) ++cpi;
}


void SwitchIndez(int w, int d)
{
   int jlat;
   REAL delphi;
   char *lpt;
   char levstr[20];

   if (Indez[w] < 0 || Indez[w] >= MaxZ[w])
   {
      if (Debug) printf("### Indez[%d] = %d\n",w,Indez[w]);
      if (Debug) printf("###  MaxZ[%d] = %d\n",w, MaxZ[w]);
      Indez[w] = 0;
   }
   if (Indez[w]+d <  0      ) return;
   if (Indez[w]+d >= MaxZ[w]) return;
   Indez[w] += d;
   lpt = strstr(WindowTitle[w],"Level");
   if (lpt)
   {
      sprintf(levstr," %d",Indez[w]+1);
      strcpy(lpt+5,levstr);
      XStringListToTextProperty(&WindowTitle[w],1,WindowName+w);
      XSetWMProperties(display,Win[w],WindowName+w,NULL,
		NULL,0,&WinSizeHints,&wm_hints,&class_hints);
      LevelChanged[w] = 1;
      return;
   }
   lpt = strstr(WindowTitle[w],"Latitude");
   if (lpt)
   {
      delphi = 180.0 / MaxZ[w];
      jlat = 90.0 -Indez[w] * delphi - delphi * 0.5;
      if (jlat > 0) sprintf(levstr," %dN", jlat);
      else          sprintf(levstr," %dS",-jlat);
      strcpy(lpt+8,levstr);
      XStringListToTextProperty(&WindowTitle[w],1,WindowName+w);
      XSetWMProperties(display,Win[w],WindowName+w,NULL,
		NULL,0,&WinSizeHints,&wm_hints,&class_hints);
      return;
   }
}


void HandleEvents(void)
{
   int w,Key,KeyIndex;
   XEvent WinEvent;

   if (nstep == 1)
   {
      RedrawControlWindow();
      Paused = 1; /* Start in pause mode */
   }


   do
   {
      /* Check for Termination */

      if (XCheckTypedWindowEvent(display,Cow,ClientMessage,&WinEvent))
      {
         /* printf("delwin %d\n",WinEvent.xclient.data.l[0]); */
         Shutdown = 1;
         return;
      }

      /* Check for user request to close windows */

      for (w=0 ; w < NumWin ; ++w)
      if (Win[w])
      {
         if (XCheckTypedWindowEvent(display,Win[w],ClientMessage,&WinEvent))
         CloseWindow(w);
      }

      /* Check for user request to close help window */

      if (HelpWindow && XCheckTypedWindowEvent(display,HelpWindow,ClientMessage,&WinEvent))
      {
         XDestroyWindow(display,HelpWindow);
         HelpWindow = 0;
      }

      /* Check for mouse click */

      if (XCheckTypedWindowEvent(display,Cow,ButtonPress,&CowEvent))
      {
         OnMouseClick();
         RedrawControlWindow();
      }

      /* Check for mouseclicks and expose events */

      for (w=0 ; w < NumWin ; ++w)
      if (Win[w])
      {
         if (XCheckTypedWindowEvent(display,Win[w],ButtonPress,&CowEvent))
         {
           if (CowEvent.xbutton.button == Button1) SwitchIndez(w,-1);
           if (CowEvent.xbutton.button == Button3) SwitchIndez(w, 1);
           if (CowEvent.xbutton.button == Button2)
           {
              if (++MapPro[w] >= MAXMAPS) MapPro[w] = 0;
              if (MapPro[w] == AZIMUTHAL &&
                 (WinAtt[w].Plot != MAPHOR && WinAtt[w].Plot != MAPTRA))
                 MapPro[w] = 0;
           }
           RedrawFlag[w] = 1;
         }
         if (XCheckTypedWindowEvent(display,Win[w],KeyPress,&CowEvent))
         {
            KeyIndex = (CowEvent.xkey.keycode - FirstKey) * SymsPerKey;
            Key = KeyMap[KeyIndex];
            if (Debug) printf("Windows %d got keyindex %d with symbol %x\n",w,KeyIndex,Key);
            if (Key == ROTATE_LEFT  && MapLR[w].f <  5) MapLR[w].f++;
            if (Key == ROTATE_RIGHT && MapLR[w].f > -5) MapLR[w].f--;
            MapLR[w].r = MapLR[w].f * RotInc;
            if (Key == XK_Up  ) SwitchIndez(w,-1);
            if (Key == XK_Down) SwitchIndez(w, 1);
         }
         if (Paused)
         if (XCheckTypedWindowEvent(display,Win[w],Expose,&CowEvent))
                     RedrawIsoWindow(w);
      }

      if (XCheckTypedWindowEvent(display,Cow,Expose,&CowEvent))
         RedrawControlWindow();

      if (HelpWindow && XCheckTypedWindowEvent(display,HelpWindow,Expose,&CowEvent))
         DrawHelpWindow();

   } while (Paused && !Shutdown);
}

void guiclose_(void)
{
   int w;

   if (Debug) printf("guiclose for MRpid %d\n",MRpid);
   if (MRpid >= 0) return; // Don't wait if multiple instances
   XSetWMName(display,Cow,&WinconName3);

   do
   {
      for (w=0 ; w < NumWin ; ++w)
      if (Win[w])
      {
         if (XCheckTypedWindowEvent(display,Win[w],Expose,&CowEvent))
                     RedrawIsoWindow(w);
      }

      if (XCheckTypedWindowEvent(display,Cow,Expose,&CowEvent))
         RedrawControlWindow();
   }
   while (!(XCheckTypedWindowEvent(display,Cow,ButtonPress,&CowEvent) &&
            HitStopButton(&CowEvent)));
     
// XCloseDisplay(display); // segmentation fault on sun compiler!
}

void SaveConfig(void)
{
   int i;
   FILE *scp;
   INTXS X,Y,xp,yp;
   unsigned int w,W,H,border,depth;
   Window Rootwin,Child;
   XWindowAttributes RootAtt;

   scp = fopen(GUI_config,"w");
   if (scp == NULL) return;

   /* Save window properties */

   fprintf(scp,"Hamburg GUI Config File Version 16\n");
   fprintf(scp,"Screen: %dx%d\n\n",ScreenW,ScreenH);
   fprintf(scp,"WinRows = %d\n",WinRows);
   fprintf(scp,"WinCols = %d\n",WinCols);
   for (w=0 ; w < NumWin ; ++w)
   {
      fprintf(scp,"\n[Window %02d]\n",w);
      fprintf(scp,"Array:%s\n",WinAtt[w].array_name);
      fprintf(scp,"Plot:%s\n",IsoNames[WinAtt[w].Plot]);
      if (WinAtt[w].Plot == ISOHOR || WinAtt[w].Plot == MAPHOR ||
          WinAtt[w].Plot == ISOTRA || WinAtt[w].Plot == MAPTRA)
      {
         fprintf(scp,"Projection:%s\n",ProNames[MapPro[w]]);
         fprintf(scp,"Rotation factor:%d\n",MapLR[w].f);
      }
      fprintf(scp,"Palette:%s\n",PalNames[WinAtt[w].Palette]);
      fprintf(scp,"Title:%s\n",WindowTitle[w]+wto);
      if (Win[w])
      {
         XGetGeometry(display,Win[w],&Rootwin,&xp,&yp,&W,&H,&border,&depth);
         XTranslateCoordinates(display,Win[w],Rootwin,0,0,&X,&Y,&Child);
         fprintf(scp,"Geometry: %4d %4d %4d %4d\n",W,H,X-WinLM,Y-WinTM);
         if (Debug)
         {
            printf("Geometry Window %d   [%8x]: %4d / %4d   %4d x %4d x %2d | %d\n",
                  w,(int)Win[w],xp,yp,W,H,depth,border);
            printf("Translated UL corner: %4d / %4d\n",X,Y);
         }
      }
      else
         fprintf(scp,"Inactive: %4d %4d %4d %4d\n",
                 WinAtt[w].w,WinAtt[w].h,WinAtt[w].x,WinAtt[w].y);
   }
   XGetGeometry(display,Cow,&Rootwin,&xp,&yp,&W,&H,&border,&depth);
   XTranslateCoordinates(display,Cow,Rootwin,0,0,&X,&Y,&Child);
   fprintf(scp,"\n[Control Window]\n");
   fprintf(scp,"Geometry: %4d %4d %4d %4d\n",W,H,X-WinLM,Y-WinTM);

   /* Scalar attributes for timeseries and tables */

   fprintf(scp,"\n# Scalar attributes for timeseries and table window\n");

   for (i=0 ; i < Parcs ; ++i)
   {
      fprintf(scp,"\n[Scalar %02d]\n",i);
      fprintf(scp,"Name:%s\n",TSName[i]);
      fprintf(scp,"Sub:%s\n",TSubsc[i]);
      fprintf(scp,"Unit:%s\n",TSUnit[i]);
      fprintf(scp,"Scale:%s\n",TScale[i]);
   }

   /* Scalar attributes for timeseries and tables */

    fprintf(scp,"\n# Parameter attributes for change menu\n");

   for (i=0 ; i < NUMSCALARS ; ++i)
   {
      if (strlen(Parc[i].Name) == 0) break;
      fprintf(scp,"\n[Parameter %02d]\n",i);
      fprintf(scp,"ParName:%s\n",Parc[i].Name);
      fprintf(scp,"ParInc:%10.4f\n",Parc[i].Inc);
      fprintf(scp,"ParMin:%10.4f\n",Parc[i].Min);
      fprintf(scp,"ParMax:%10.4f\n",Parc[i].Max);
   }

   fclose(scp);
}


void lp2ps(REAL *a, REAL *b)
{
   int i,j,k,l,k1,l1;
   REAL xnp,ynp,xsp,ysp,dx,dy,x,y,lfa,ua,ub;

   lfa = (DimX-1) / (2.0 * M_PI);
   xnp = (DimX-1) * 0.25;
   ynp = (DimY-1) * 0.50;
   xsp = (DimX-1) * 0.75;
   ysp = ynp;

   /* Northern hemisphere */

   for (j = 0 ; j <  DimY   ; ++j)
   for (i = 0 ; i <= DimX/2 ; ++i)
   {
      dx = xnp - i;
      dy = ynp - j;
      x = atan2(dx,dy) * lfa;
      y = sqrt(dx * dx + dy * dy);
      if (x < 0.0) x += (DimX-1);
      k = x;
      l = y;
      k1 = k + 1;
      if (k1 >= DimX) k1 = 0;
      if (l <= DimY/2+1)
      {
         ua = (k+1-x) *  a[k+l*DimX] + (x-k) * a[k1+l*DimX];
         ub = (k+1-x) * a[k+(l+1)*DimX] + (x-k) * a[k+1+(l+1)*DimX];
         b[i+j*DimX] = (l+1-y) * ua + (y-l) * ub;
      }
      else b[i+j*DimX] = - 99999.0;
   }

   /* Souhern hemisphere */

   for (j = 0        ; j < DimY ; ++j)
   for (i = DimX/2+1 ; i < DimX ; ++i)
   {
      dx = i - xsp;
      dy = ysp - j;
      x = atan2(dx,dy) * lfa;
      y = (DimY-1) - sqrt(dx * dx + dy * dy);
      if (x < 0.0) x += (DimX-1);
      k = x;
      l = y;
      k1 = k + 1;
      if (k1 >= DimX) k1 = 0;
      if (l >= DimY/2-3)
      {
         ua = (k+1-x) *  a[k+l*DimX] + (x-k) * a[k1+l*DimX];
         ub = (k+1-x) * a[k+(l+1)*DimX] + (x-k) * a[k+1+(l+1)*DimX];
         b[i+j*DimX] = (l+1-y) * ua + (y-l) * ub;
      }
      else b[i+j*DimX] = - 99999.0;
   }
}


/* ======================================================= */
/* lp2az - transform lon/lat array to azimuthal projection */
/* ======================================================= */

void lp2az(REAL *a, REAL *b, float laz)
{
   int x,y,dxc,k,k1,l;
   REAL dx,dy,lam,rho,rad,phi,ua,ub,l00,p00,xpi,ypi;

   rad = DimY >> 1;
   dxc = DimX >> 1;
   xpi = DimX / M_PI * 0.5;
   ypi = DimY / M_PI;
   l00 =  (int)((laz * DimX) / 360 + dxc) % DimX;
   p00 = DimY / 2;

   for (y = 0 ; y < DimY ; ++y)
   {
      dy  = y - rad;
      phi = y;
      for (x = 0 ; x < DimX ; ++x)
      {
         dx = x - dxc;
         rho = sqrt(dx * dx + dy * dy);
         if (rho < rad)
         {
            lam = l00 + xpi * atan2(dx / rad, cos(asin(rho / rad)));
            phi = p00 + ypi *  asin(dy / rad);
            k = lam;
            k1 = k + 1;
            if (k1 >= DimX) k1 = 0;
            l = phi;
            ua = (k+1-lam) * a[k+ l   *DimX] + (lam-k) * a[k1 + l   *DimX];
            ub = (k+1-lam) * a[k+(l+1)*DimX] + (lam-k) * a[k+1+(l+1)*DimX];
            b[x+y*DimX] = (l+1-y) * ua + (y-l) * ub;
         }
         else b[x+y*DimX] = - 99999.0;
      }
   }
}

void ShowGridCS(void)
{
   int jlev,jlat,x,y,len;
   float dx,dy;
   char Text[80];

   /* Grid for zonal mean cross sections */

   XSetForeground(display,gc,WhitePix);
   XSetBackground(display,gc,BlackPix);
   XSetLineAttributes(display,gc,1,LineOnOffDash,CapButt,JoinRound);
   dy = (InYSize - OffY) / (DimY - 1.0);
   for (jlev = 1 ; jlev < DimY-1 ; ++jlev)
   {
      y = OffY + jlev * dy;
      XDrawLine(display,pix,gc,OffX,y,InXSize-1,y);
   }
   dx = (InXSize - OffX) / 6.0; /* Every 30 degrees */
   for (jlat = 1 ; jlat < 6 ; ++jlat)
   {
      x = OffX + jlat * dx;
      XDrawLine(display,pix,gc,x,OffY,x,InYSize-1);
   }
   if (GridLabel)
   {
      XDrawImageString(display,pix,gc,InXSize/2-FixFontWidth,InYSize-FixFont->descent,"EQ",2);
      XDrawImageString(display,pix,gc,OffX,InYSize-FixFont->descent,"N",1);
      XDrawImageString(display,pix,gc,InXSize-FixFontWidth,InYSize-FixFont->descent,"S",1);
   }
   XSetLineAttributes(display,gc,1,LineSolid,CapButt,JoinRound);
}

void ShowGridCol(void)
{
   int jlev,jtim,x,y,len;
   float dx,dy;
   char Text[80];

   /* Grid for column time series */

   XSetForeground(display,gc,WhitePix);
   XSetBackground(display,gc,BlackPix);
   XSetLineAttributes(display,gc,1,LineOnOffDash,CapButt,JoinRound);
   dy = (InYSize - OffY) / (DimY - 1.0);
   for (jlev = 1 ; jlev < DimY-1 ; ++jlev)
   {
      y = OffY + jlev * dy;
      XDrawLine(display,pix,gc,OffX,y,InXSize-1,y);
   }
   dx = (InXSize - OffX) / 4.0; /* 4 slices */
   for (jtim = 1 ; jtim < 4 ; ++jtim)
   {
      x = OffX + jtim * dx;
      XDrawLine(display,pix,gc,x,OffY,x,InYSize-1);
   }
   if (GridLabel)
   {
      XDrawImageString(display,pix,gc,InXSize/2-3*FixFontWidth/2,InYSize-FixFont->descent,"t-2",3);
      XDrawImageString(display,pix,gc,InXSize/4-3*FixFontWidth/2,InYSize-FixFont->descent,"t-3",3);
      XDrawImageString(display,pix,gc,3*InXSize/4-3*FixFontWidth/2,InYSize-FixFont->descent,"t-1",3);
   }
   XSetLineAttributes(display,gc,1,LineSolid,CapButt,JoinRound);
}

void ShowGridLonsi(void)
{
   int jlev,jlat,x,y,len;
   float dx,dy;
   char Text[80];

   /* Grid for Longitude Sigma section */

   XSetForeground(display,gc,WhitePix);
   XSetBackground(display,gc,BlackPix);
   XSetLineAttributes(display,gc,1,LineOnOffDash,CapButt,JoinRound);
   dy = (InYSize - OffY) / (DimY - 1.0);
   for (jlev = 1 ; jlev < DimY-1 ; ++jlev)
   {
      y = OffY + jlev * dy;
      XDrawLine(display,pix,gc,OffX,y,InXSize-1,y);
   }
   dx = (InXSize - OffX) / 6.0; /* Every 60 degrees */
   for (jlat = 1 ; jlat < 6 ; ++jlat)
   {
      x = OffX + jlat * dx;
      XDrawLine(display,pix,gc,x,OffY,x,InYSize-1);
   }
   if (GridLabel)
   {
      XDrawImageString(display,pix,gc,InXSize/2-3*FixFontWidth/2,InYSize-FixFont->descent,"180",3);
      XDrawImageString(display,pix,gc,OffX,InYSize-FixFont->descent,"0",1);
      XDrawImageString(display,pix,gc,InXSize-3*FixFontWidth,InYSize-FixFont->descent,"360",3);
   }
   XSetLineAttributes(display,gc,1,LineSolid,CapButt,JoinRound);
}

void ShowGridCyl(void)
{
   int jlon,jlat,x,y;
   float dx,dy;

   /* Grid for cylinder projection */

   XSetForeground(display,gc,WhitePix);
   XSetBackground(display,gc,BlackPix);
   XSetLineAttributes(display,gc,1,LineOnOffDash,CapButt,JoinRound);
   dy = (InYSize - OffY) / 6.0;
   for (jlat = 1 ; jlat < 6 ; ++jlat)
   {
      y = OffY + jlat * dy;
      XDrawLine(display,pix,gc,OffX,y,InXSize-1,y);
   }
   dx = (InXSize - OffX) / 6.0; /* Every 30 degrees */
   for (jlon = 1 ; jlon < 6 ; ++jlon)
   {
      x = OffX + jlon * dx;
      XDrawLine(display,pix,gc,x,OffY,x,InYSize-1);
   }
   if (GridLabel)
   {
      XDrawImageString(display,pix,gc,OffX,OffY+(InYSize-OffY)/2+FixFontHeight/2-FixFont->descent,"EQ",2);
      XDrawImageString(display,pix,gc,OffX,OffY+FixFontHeight-FixFont->descent,"N",1);
      XDrawImageString(display,pix,gc,OffX,InYSize-FixFont->descent,"-180",4);
      XDrawImageString(display,pix,gc,InXSize-3*FixFontWidth,InYSize-FixFont->descent,"180",3);
      XDrawImageString(display,pix,gc,OffX+(InXSize-OffX)/2-3*FixFontWidth/2,InYSize-FixFont->descent,"0",1);
   }
   XSetLineAttributes(display,gc,1,LineSolid,CapButt,JoinRound);
}

void ShowGridHov(void)
{
   int jlon,jlat,x,y,len;
   float dx,dy;
   char Text[80];

   /* Grid for Hovmoeller */

   XSetForeground(display,gc,WhitePix);
   XSetBackground(display,gc,BlackPix);
   XSetLineAttributes(display,gc,1,LineOnOffDash,CapButt,JoinRound);
   dx = (InXSize - OffX) / 6.0;
   for (jlon = 1 ; jlon < 6 ; ++jlon)
   {
      x = OffX + jlon * dx;
      XDrawLine(display,pix,gc,x,OffY,x,InYSize-1);
   }
   dy = (InYSize - OffY) / 5.0; /* Every 5 days */
   for (jlat = 1 ; jlat < 5 ; ++jlat)
   {
      y = OffY + jlat * dy;
      XDrawLine(display,pix,gc,OffX,y,InXSize-1,y);
   }
   if (GridLabel)
   {
      XDrawImageString(display,pix,gc,OffX,InYSize-FixFont->descent,"0",1);
      XDrawImageString(display,pix,gc,InXSize-3*FixFontWidth,InYSize-FixFont->descent,"360",3);
      XDrawImageString(display,pix,gc,OffX+(InXSize-OffX)/2-3*FixFontWidth/2,InYSize-FixFont->descent,"180",3);
      XDrawImageString(display,pix,gc,OffX,OffY+FixFontHeight-FixFont->descent,"t0",2);
      XDrawImageString(display,pix,gc,OffX,OffY+2*InYSize/5+FixFontHeight/2-FixFont->descent,"t-2",3);
      XDrawImageString(display,pix,gc,OffX,OffY+4*InYSize/5+FixFontHeight/2-FixFont->descent,"t-4",3);
   }
   XSetLineAttributes(display,gc,1,LineSolid,CapButt,JoinRound);
}

void ShowGridHovT(void)
{
   int jlon,jlat,x,y,len;
   float dx,dy;
   char Text[80];

   /* Grid for Hovmoeller */

   XSetForeground(display,gc,WhitePix);
   XSetBackground(display,gc,BlackPix);
   XSetLineAttributes(display,gc,1,LineOnOffDash,CapButt,JoinRound);
   dy = (InYSize - OffY) / 6.0;
   for (jlat = 1 ; jlat < 6 ; ++jlat)
   {
      y = OffY + jlat * dy;
      XDrawLine(display,pix,gc,OffX,y,InXSize-1,y);
   }
   dx = (InXSize - OffX) / 5.0; /* Every 5 days */
   for (jlat = 1 ; jlat < 5 ; ++jlat)
   {
      x = OffX + jlat * dx;
      XDrawLine(display,pix,gc,x,OffY,x,InYSize-1);
   }
   if (GridLabel)
   {
      XDrawImageString(display,pix,gc,OffX,OffY+FixFontHeight-FixFont->descent,"0",1);
      XDrawImageString(display,pix,gc,OffX,OffY+(InYSize-OffY)/2+FixFontHeight/2-FixFont->descent,"180",3);
      XDrawImageString(display,pix,gc,OffX,InYSize-FixFont->descent,"360",3);
      XDrawImageString(display,pix,gc,InXSize-2*FixFontWidth,InYSize-FixFont->descent,"t0",2);
      XDrawImageString(display,pix,gc,3*InXSize/5-3*FixFontWidth/2,InYSize-FixFont->descent,"t-2",3);
      XDrawImageString(display,pix,gc,  InXSize/5-3*FixFontWidth/2,InYSize-FixFont->descent,"t-4",3);
   }
   XSetLineAttributes(display,gc,1,LineSolid,CapButt,JoinRound);
}

void ShowGridPolar(void)
{
   int dx,dy,x,y,xh,yh,ox,oy;
   XPoint pxy[3];

   /* Grid for polar projection */

   dx = (InXSize - OffX) / 2;
   dy = (InYSize - OffY);
   xh = dx / 2;
   yh = dy / 2;
   ox = dx / 3.414;
   oy = dy / 3.414;

   if (Grid)
   {
      XSetForeground(display,gc,WhitePix);
      XSetBackground(display,gc,BlackPix);
      XSetLineAttributes(display,gc,1,LineOnOffDash,CapButt,JoinRound);
   
      /* Northern Hemisphere */
   
      XDrawArc(display,pix,gc,OffX+dx/6,OffY+dy/6,2*dx/3,2*dy/3,0,360*64);
      XDrawArc(display,pix,gc,OffX+dx/3,OffY+dy/3,dx/3,dy/3,0,360*64);
   
      XDrawLine(display,pix,gc,OffX,yh,InXSize,yh);
      XDrawLine(display,pix,gc,OffX+xh,OffY,OffX+xh,InYSize);
   
      /* Southern Hemisphere */
   
      XDrawArc(display,pix,gc,OffX+7*dx/6,OffY+dy/6,2*dx/3,2*dy/3,0,360*64);
      XDrawArc(display,pix,gc,OffX+4*dx/3,OffY+dy/3,dx/3,dy/3,0,360*64);
   
      x = OffX + dx + xh;
      XDrawLine(display,pix,gc,x,OffY,x,InYSize);
   
      if (GridLabel)
      {
         x  = OffX + xh - FixFontWidth/2;
         y  = OffY + yh + FixFontHeight/2 - FixFont->descent;
         XDrawImageString(display,pix,gc,x,y,"N",1);
         XDrawImageString(display,pix,gc,x+dx,y,"S",1);
      }
      XSetLineAttributes(display,gc,1,LineSolid,CapButt,JoinRound);
   }

   /* Octagon mask */

   XSetForeground(display,gc,BlackPix);
   pxy[0].x = pxy[1].x = OffX;
   pxy[0].y = pxy[2].y = OffY;
   pxy[1].y = OffY + oy;
   pxy[2].x = OffX + ox;
   XFillPolygon(display,pix,gc,pxy,3,Convex,CoordModeOrigin);
   pxy[0].x = OffX + dx - ox;
   pxy[1].x = OffX + dx;
   pxy[2].x = OffX + dx + ox;
   XFillPolygon(display,pix,gc,pxy,3,Convex,CoordModeOrigin);
   pxy[0].x = InXSize - ox;
   pxy[1].x = InXSize;
   pxy[2].x = InXSize;
   XFillPolygon(display,pix,gc,pxy,3,Convex,CoordModeOrigin);

   pxy[0].x = pxy[1].x = OffX;
   pxy[0].y = pxy[2].y = InYSize;
   pxy[1].y = InYSize - oy;
   pxy[2].x = OffX + ox;
   XFillPolygon(display,pix,gc,pxy,3,Convex,CoordModeOrigin);
   pxy[0].x = OffX + dx - ox;
   pxy[1].x = OffX + dx;
   pxy[2].x = OffX + dx + ox;
   XFillPolygon(display,pix,gc,pxy,3,Convex,CoordModeOrigin);
   pxy[0].x = InXSize - ox;
   pxy[1].x = InXSize;
   pxy[2].x = InXSize;
   XFillPolygon(display,pix,gc,pxy,3,Convex,CoordModeOrigin);
}


void AutoPalette(int w, struct ColorStrip Strip[], REAL f[], int Dim)
{
   int i,nbands;
   double fmin, fmax, frange, delta, fdelta, xdelta, Lo;

   fmin = fmax = f[0];
   for (i=0; i < Dim ; ++i)
   {
      if (fmin > f[i]) fmin = f[i];
      if (fmax < f[i]) fmax = f[i];
   }
   frange = fmax - fmin;
   Lo     = 0.0;
   xdelta = 0.1;
   if (frange > 1.0e-10)
   {
       delta  = frange / AUTOCOLORS;
       if (delta > AutoDelta[w] / 1.4 && delta < 1.4 * AutoDelta[w])
       {
          xdelta = AutoXdelt[w];
          delta  = AutoDelta[w];
          Lo     = AutoLo[w];
       }
       else
       {
          fdelta = pow(10.0,rint(log10(delta)));
          xdelta = fdelta;
          nbands = frange / xdelta;
          if (nbands < AUTOCOLORS / 2) xdelta = 0.5 * fdelta;
          if (nbands > AUTOCOLORS    ) xdelta = 2.0 * fdelta;
          nbands = frange / xdelta;
          if (nbands > AUTOCOLORS    ) xdelta = 5.0 * fdelta;
          Lo = xdelta * floor(fmin / xdelta);
          AutoDelta[w] = delta;
          AutoXdelt[w] = xdelta;
          AutoLo[w]    = Lo;
/*
          if (Debug)
          {
             printf("Autopalette\n");
             printf(" Range: %14.6e\n",frange);
             printf(" delta: %14.6e\n",delta);
             printf("fdelta: %14.6e\n",fdelta);
             printf("xdelta: %14.6e\n",xdelta);
          }
*/
       }
   }
   Strip[0].Lo = Lo;
   for (i=0 ; i < AUTOCOLORS ; ++i)
   {
       Strip[i+1].Lo = Strip[i].Hi = Strip[i].Lo + xdelta;
       // if (Debug) printf("Auto[%2d] = %10.4f - %10.4f\n",i,Strip[i].Lo,Strip[i].Hi);
   }
   return;
}


// Interface routine to model code in FORTRAN

void guiput_(char *aname, float *array, int *dimx, int *dimy, int *dimz)
{
   int i,nb,nf;

   if (Debug) printf("guiput(%s,%12.4e,%d,%d,%d)\n",
                     aname,*array,*dimx,*dimy,*dimz);
   nf = (*dimx) * (*dimz);
   if (*dimy > 0) nf *= (*dimy);
   else           nf *= 2;
   nb = nf * sizeof(float);
   for (i=0 ; i < NumArrays ; ++i)
      if (!strcmp(aname,Array[i].Name)) break;
   if (i == NumArrays)
   {
      ++NumArrays;
      strcpy(Array[i].Name,aname);
      Array[i].Data = FloatAlloc(nf,aname);
      Array[i].DimX = *dimx;
      Array[i].DimY = *dimy;
      Array[i].DimZ = *dimz;
      if (strcmp(aname,"GU") == 0) Uindex = i;
      if (strcmp(aname,"GV") == 0) Vindex = i;
   }
   if (Array[i].Data) memcpy(Array[i].Data,array,nb);
   Array[i].Flag = 1; // Data have changed
}


void guisend_(char *aname, float *array, int *dimx, int *dimy, int *dimz)
{
   int i,nb,nf;

   nf = (*dimx) * (*dimz);
   if (*dimy > 0) nf *= (*dimy);
   else           nf *= 2;
   nb = nf * sizeof(float);
   for (i=0 ; i < NumArrays ; ++i)
      if (!strcmp(aname,Array[i].Name)) break;
   if (i == NumArrays)
   {
      ++NumArrays;
      strcpy(Array[i].Name,aname);
      Array[i].DimX = *dimx;
      Array[i].DimY = *dimy;
      Array[i].DimZ = *dimz;
   }
   Array[i].Data = array;
   Array[i].Flag = 1; // Data have changed
}


// Do the plot

void iso(int w,int PicType,REAL *field,int dimx,int dimy,int dimz,int pal)
{
   char Text[128];

   int i,j,k,len,lens,xp,yp,status,x;
   int y,dx,r,width,height;
   INTXU border,depth;
   REAL f,o,ra,rb;
   int CapLines;
   float *tspt;
   XEvent WinEvent;
   Window Rootwin,Child;

   // if (Debug) printf("iso(%d,%s,%12.4e,%d,%d,%d,%d)\n",w,IsoNames[PicType],*field,dimx,dimy,dimz,pal);
   if (Win[w] == 0) return;
   win = w;
   if (SkipFreq > 1 && (nstep % SkipFreq) != 0 && 
      (PicType == ISOCS || PicType == ISOHOR || PicType == MAPHOR || PicType == ISOCOL)) return;
   XGetGeometry(display,Win[w],&Rootwin,&xp,&yp,&WinXSize,&WinYSize,
		&border,&depth);
   WinAtt[w].w = WinXSize;
   WinAtt[w].h = WinYSize;
   InXSize = WinXSize - OffX;
   InYSize = WinYSize - OffY;
   if (PicType != ISOTRA && PicType != MAPTRA) InYSize -= 20; // Room for colorbar
   DimX  = dimx;
   DimY  = dimy;
   DimZ  = dimz;
   if (DimY < 0) DimY = -DimY; // Get NTP1 for ISOSH
   DimXY = DimX * DimY;
   VGAX  = (InXSize-1.0) / (DimX-1.0);
   if (DimY > 1) VGAY  = (InYSize-1.0) / (DimY-1.0);
   else          VGAY  = 1.0;
   Field = field;

   if (PicType == ISOTRA || PicType == MAPTRA)
   {
      DimX  = dimx; // NLON
      DimY  = dimy; // NLAT
      DimZ  = dimz; // NLEV
      DimXY = DimX * DimY;
      VGAX  = (InXSize-1.0) / (DimX-1.0);
      VGAY  = (InYSize-1.0) / (DimY-1.0);
      if (MaxZ[w] < DimZ)
      {
         MaxZ[w]   = DimZ;
         Indez[w]  = DimZ / 4;
         SwitchIndez(w,0); // Initialize title
      }
   }

   if (PicType == ISOHOR || PicType == MAPHOR)
   {
      DimX  = dimx; // NLON
      DimY  = dimy; // NLAT
      DimZ  = dimz; // NLEV
      DimXY = DimX * DimY;
      VGAX  = (InXSize-1.0) / (DimX-1.0);
      VGAY  = (InYSize-1.0) / (DimY-1.0);
      if (MaxZ[w] < DimZ)
      {
         MaxZ[w]   = DimZ;
         Indez[w]  = DimZ / 4;
         SwitchIndez(w,0); // Initialize title
         TSdata[w] = FloatAlloc(DimXY,"ISOLON");
      }
      if (MapPro[w] == POLAR)
      {
         Field = TSdata[w];
         lp2ps(field + Indez[w] * DimXY,Field);
      }
/*
      else if (MapPro[w] == AZIMUTHAL)
      {
         Field = TSdata[w];
         lp2az(field + Indez[w] * DimXY,Field,MapLR[w].l);
      }
*/
      else Field = field + Indez[w] * DimXY;
   }

   if (PicType == ISOLON)
   {
      DimX  = dimx; // NLON
      DimY  = dimz; // NLEV
      DimZ  = dimy; // NLAT
      DimXY = DimX * DimY;
      VGAX  = (InXSize-1.0) / (DimX-1.0);
      VGAY  = (InYSize-1.0) / (DimY-1.0);
      if (!TSdata[w])
      {
         TSdata[w] = FloatAlloc(DimXY,"ISOLON");
         MaxZ[w]   = DimZ;
         Indez[w]  = DimZ / 4;
         SwitchIndez(w,0); // Initialize title
      }
      Field = TSdata[w];
      for (j=0,i=0,k=DimX*Indez[w] ; j < DimY ; ++j,i+=DimX,k+=DimX*DimZ)
      {
         memcpy(Field+i,field+k,DimX * sizeof(REAL)); // copy latitude
      }
   }

   if (PicType == ISOHOV)
   {
      DimX  = dimx; // NLON
      DimY  = DimT; // time
      DimZ  = dimy; // NLAT
      DimXY = DimX * DimY;
      VGAX  = (InXSize-1.0) / (DimX-1.0);
      VGAY  = (InYSize-1.0) / (DimY-1.0);
      if (!TSdata[w])
      {
         TSdata[w] = FloatAlloc(DimXY,"ISOHOV");
         MaxZ[w]   = DimZ;
         Indez[w]  = DimZ / 4;
         SwitchIndez(w,0); // Initialize title
      }
      Field = TSdata[w];
      memmove(Field+DimX,Field,(DimXY-DimX) * sizeof(REAL)); // scroll array
      memcpy(Field,field+Indez[w]*DimX,DimX * sizeof(REAL)); // add line
   }

   // Advance write pointer HovInx until end of DimX
   // then scroll array and reset pointer

   if (PicType == ISOTS)
   {
      DimX  = DimT; // time
      DimY  = dimx; // # of variables
      DimZ  = dimz;
      DimXY = DimX * DimY;
      VGAX  = (InXSize-1.0) / (DimX-1.0);
      VGAY  = (InYSize-1.0) / (DimY-1.0);
      if (!TSdata[w])
      {
         TSdata[w] = FloatAlloc(DimXY+DimX,"ISOTS");
      }
      tspt  = TSdata[w];
      ++HovInx[w];
      if (HovInx[w] >= DimX)
      {
          HovInx[w] = 0;
          memmove(tspt,tspt+DimX,DimXY * sizeof(REAL)); // scroll array
      }
      for (i=HovInx[w]+DimX-1,j=0 ; j < DimY ; i+=DimX,++j)
         tspt[i] = field[j];                                  // add column
      Field = tspt+HovInx[w];
   }

   if (PicType == ISOCOL)
   {
      DimX  = DimT; // time
      DimY  = dimx; // level
      DimZ  = dimy; // clickable index
      DimXY = DimX * DimY;
      VGAX  = (InXSize-1.0) / (DimX-1.0);
      VGAY  = (InYSize-1.0) / (DimY-1.0);
      if (!TSdata[w])
      {
         TSdata[w] = FloatAlloc(DimXY+DimX,"ISOCOL");
         MaxZ[w]   = DimZ;
         Indez[w]  = DimZ / 2;
         SwitchIndez(w,0); // Initialize title
      }
      tspt  = TSdata[w];
      ++HovInx[w];
      if (HovInx[w] >= DimX)
      {
          HovInx[w] = 0;
          memmove(tspt,tspt+DimX,DimXY * sizeof(REAL)); // scroll array
      }
      for (i=HovInx[w]+DimX-1,j=0 ; j < DimY ; i+=DimX,++j)
         tspt[i] = field[j+Indez[w]*DimY];                                  // add column
      Field = tspt+HovInx[w];
   }

   if (DimX > FlagSize)
   {
      if (Flag) free(Flag);
      Flag = IntAlloc(DimX,"Flag");
      FlagSize = DimX;
   }

   if (PicType == ISOSH) pal = 1; // ISOSH has its own palette
   if (pal < 1 || pal >= NUMPAL)
   {
      Cstrip = Autostrip;
      AutoPalette(w,Cstrip,Field,DimXY);
      Lines  = AUTOCOLORS;
      // printf("Autopalette %d for Win %d\n",pal,w);
   }
   else
   {
      Cstrip = Pallet[pal];
      Lines  = LineCo[pal];
   }

   SizeChanged = (WinPixMap[w].DimX != WinXSize || WinPixMap[w].DimY != WinYSize);
   if (SizeChanged)
   {
      if (WinPixMap[w].Pix) XFreePixmap(display,WinPixMap[w].Pix);
      WinPixMap[w].Pix = XCreatePixmap(display,Win[w],WinXSize,WinYSize,ScreenD);
      if (Debug)
         printf("CreatePixmap  %10x %6d bytes\n",
               (unsigned int)WinPixMap[w].Pix,WinXSize*WinYSize);
      WinPixMap[w].DimX = WinXSize;
      WinPixMap[w].DimY = WinYSize;
   }
   pix = WinPixMap[w].Pix; /* Set current pixmap */

   /* Draw colour bar */
   
   if ((SizeChanged || Cstrip == Autostrip) &&
       (PicType == ISOCS  || PicType == ISOHOR ||
       	PicType == ISOHOV || PicType == ISOLON ||
        PicType == ISOCOL || PicType == MAPHOR ))
   {
      XSetForeground(display,gc,BlackPix);
      XFillRectangle(display,pix,gc,0,InYSize,WinXSize,WinYSize);

      CapLines = WinXSize / (20 + 3 * FixFontWidth);
      if (CapLines > Lines) CapLines = Lines;
      for (i=0 ; i < CapLines ; ++i)
      {
         XSetForeground(display,gc,Cstrip[i].pixel);
         XFillRectangle(display,pix,gc,OffX+5+i*(WinXSize-20)/(CapLines-1),
			OffY+InYSize+5,10,10);
         XSetForeground(display,gc,BlackPix);
         XDrawRectangle(display,pix,gc,OffX+5+i*(WinXSize-20)/(CapLines-1),
			OffY+InYSize+5,10,10);
      }
      XSetForeground(display,gc,WhitePix);
      XSetBackground(display,gc,BlackPix);
      rb = Cstrip[CapLines-2].Hi - Cstrip[0].Hi;
      for (i=0 ; i < CapLines-1 ; ++i)
      {
	 Text[0] = 0;
	 ra = Cstrip[i].Hi;
	 j = ra;
	 if (j > -1000 && j < 10000) sprintf(Text,"%d",j);
         if (ra < 100.0 && ra > -9.99)
         {
            sprintf(Text,"%4.1f",ra);
            if (!strcmp(Text+2,".0")) Text[2] = 0;
         }
	 if (ra < 10.0 && ra >= 0.0)
         {
            sprintf(Text,"%4.2f",ra);
            if (!strcmp(Text+1,".00")) Text[1] = 0;
         }
         len    = strlen(Text);
	 if (len)
         {
            width  = XTextWidth(FixFont,Text,len);
            height = FixFont->ascent + FixFont->descent;
            xp     = OffX + 10 + (i+0.5) * (WinXSize-20)/(CapLines-1) - width/2;
            yp     = WinYSize - height + 10;
            XDrawImageString(display,pix,gc,xp,yp,Text,len);
	 }
      }
   }

   /* Draw mode legend */
   
   if (SizeChanged && PicType == ISOSH)
   {
      dx = WinXSize / 23;
      XSetForeground(display,gc,BlackPix);
      XFillRectangle(display,pix,gc,0,0,WinXSize,WinYSize);

      XSetForeground(display,gc,LightGreen.pixel);
      XSetBackground(display,gc,BlackPix);
      for (i=0 ; i < 21 ; i+=2)
      {
         sprintf(Text,"%d",i);
         len    = strlen(Text);
         width  = XTextWidth(FixFont,Text,len);
         height = FixFont->ascent + FixFont->descent;
         xp     = dx + i * dx - width/2 - FixFontWidth/2 + 1;
         yp     = FixFontHeight;
         XDrawImageString(display,pix,gc,xp,yp,Text,len);
      }
      XSetForeground(display,gc,LightBlue.pixel);
      for (i=0 ; i < 21 ; i+=2)
      {
         sprintf(Text,"%d",i);
         len    = strlen(Text);
         width  = XTextWidth(FixFont,Text,len);
         height = FixFont->ascent + FixFont->descent;
         xp     = WinXSize - width - 2;
         yp     = 2 * FixFontHeight + i * dx;
         if (yp+FixFont->descent > InYSize) break;
         XDrawImageString(display,pix,gc,xp,yp,Text,len);
      }
      strcpy(Text,"n/m");
      len    = strlen(Text);
      width  = XTextWidth(FixFont,Text,len);
      xp     = WinXSize - width - 2;
      yp     = FixFontHeight;
      XSetForeground(display,gc,WhitePix);
      XDrawImageString(display,pix,gc,xp,yp,Text,len);
      strcpy(Text,"High ");
      len    = strlen(Text);
      width  = XTextWidth(FixFont,Text,len);
      xp     = WinXSize/2 - AMPLI_COLS * 10 - width;
      yp     = WinYSize - FixFontHeight + 10;
      r      = FixFontHeight-2;
      XSetForeground(display,gc,WhitePix);
      if (xp > 0) XDrawImageString(display,pix,gc,xp,yp,Text,len);
      strcpy(Text," Low");
      len    = strlen(Text);
      width  = XTextWidth(FixFont,Text,len);
      xp     = WinXSize/2 + AMPLI_COLS * 10;
      if (xp + width < WinXSize) XDrawImageString(display,pix,gc,xp,yp,Text,len);
      yp     = InYSize;
      XDrawLine(display,pix,gc,OffX,yp,WinXSize,yp);
      yp     = WinYSize - r -  FixFont->descent;

      for (i=0 ; i < AMPLI_COLS ; ++i)
      {
         xp = WinXSize/2 - AMPLI_COLS * 10 + i * 20;
         XSetForeground(display,gc,AmpliStrip[AMPLI_COLS-i-1].pixel);
         XFillArc(display,pix,gc,xp,yp,r,r,0,360*64);
      }
   }

   if (SizeChanged && PicType == ISOTS) /* Timeseries Caption */
   {
      XSetForeground(display,gc,BlackPix);
      XFillRectangle(display,pix,gc,0,InYSize,WinXSize,WinYSize);
      XSetBackground(display,gc,BlackPix);

      xp = FixFontWidth;
      yp = WinYSize - FixFontHeight + 10;
      for (j=0 ; j < DimY ; ++j)
      {
         if (TSName[j][0])
         {
            strcpy(Text,TSName[j]);
            len    = strlen(Text);
            width  = XTextWidth(FixFont,Text,len);
            XSetForeground(display,gc,TSColor[j]);
            XDrawImageString(display,pix,gc,xp,yp,Text,len);
   	    xp += width;
         }
         if (TSubsc[j][0])
         {
            XSetFont(display, gc, SubFont->fid);
            strcpy(Text,TSubsc[j]);
            len    = strlen(Text);
            width  = XTextWidth(SubFont,Text,len);
            XSetForeground(display,gc,TSColor[j]);
            XDrawImageString(display,pix,gc,xp,yp,Text,len);
   	    xp += width;
            XSetFont(display, gc, FixFont->fid);
         }
         xp += FixFontWidth;
      }
   }

   if (SizeChanged && PicType == ISOTAB) /* Table Caption */
   {
      XSetForeground(display,gc,BlackPix);
      XFillRectangle(display,pix,gc,0,0,WinXSize,WinYSize);
      XSetBackground(display,gc,BlackPix);
   }

   if (PicType == ISOTS) 
   {
      if (TSxp[w] == NULL) TSxp[w] = SizeAlloc(DimX * DimY , sizeof(XPoint),"TSxp");
      if (Dmin[w] == NULL) Dmin[w] = FloatAlloc(DimY ,"Dmin");
      if (Dmax[w] == NULL) Dmax[w] = FloatAlloc(DimY ,"Dmax");
      XSetForeground(display,gc,BlackPix);
      XFillRectangle(display,pix,gc,0,0,WinXSize,InYSize);

      for (j=0 ; j < DimY ; ++j) Dmin[w][j] = Dmax[w][j] = field[j];
      if (nstep > 2)
      {
         for (j=0 ; j < DimY ; ++j)
         {
            for (i=1 ; i < DimX ; ++i)
            {
               Dmin[w][j] = MIN(Dmin[w][j],Field[i+j*DimX]);
               Dmax[w][j] = MAX(Dmax[w][j],Field[i+j*DimX]);
            }
         }
      }

      if (nstep > 2)
      for (j=0 ; j < DimY ; ++j)
      {
         XSetForeground(display,gc,TSColor[j]);
         o = Dmin[w][j];
         if ((Dmax[w][j] - Dmin[w][j]) > 1.0e-20) f = (InYSize-2) / (Dmax[w][j] - Dmin[w][j]);
         else f = 1.0;
   
         for (i=1 ; i < DimX ; ++i)
         {
            TSxp[w][i].x = VGAX * i;
   	    TSxp[w][i].y = InYSize - 1 - f * (Field[i+j*DimX] - o);
         }
         XDrawLines(display,pix,gc,TSxp[w]+1,DimX-1,CoordModeOrigin);
      }
   }

   if (PicType == ISOTAB) 
   {
      XSetForeground(display,gc,BlackPix);
      XSetBackground(display,gc,BlackPix);
      XFillRectangle(display,pix,gc,0,0,WinXSize,InYSize);
      XSetFont(display, gc, BigFont->fid);
      yp = BigFontHeight;
      for (j=0 ; j < DimX ; ++j)
      if (TSName[j][0])
      {
         XSetForeground(display,gc,TSColor[j]);
         xp = BigFontWidth;
         strcpy(Text,TSName[j]);
         len = strlen(Text);
         XDrawImageString(display,pix,gc,xp,yp,Text,len);
   	 xp += XTextWidth(BigFont,Text,len);
         if (TSubsc[j][0])
         {
            XSetFont(display, gc, FixFont->fid);
            XDrawImageString(display,pix,gc,xp,yp,TSubsc[j],strlen(TSubsc[j]));
            XSetFont(display, gc, BigFont->fid);
         }
         xp = 8 * BigFontWidth;
         if (TSUnit[j][0])
         {
   	    sprintf(Text,"= %7.3f [%s]",field[j],TSUnit[j]);
            if (TScale[j][0]) strcat(Text," 10");
            len = strlen(Text);
            XDrawImageString(display,pix,gc,xp,yp,Text,len);
            if (TScale[j][0] > ' ')
            {
               lens = strlen(TScale[j]);
               XSetFont(display, gc, FixFont->fid);
               XDrawImageString(display,pix,gc,xp+len*BigFontWidth,yp-BigFontHeight+FixFontHeight,TScale[j],lens);
               XSetFont(display, gc, BigFont->fid);
            }
         }
   	 yp += BigFontHeight;
      }
      XSetFont(display, gc, FixFont->fid);
   }

   if (PicType == ISOTRA)
   {
      XSetForeground(display,gc,BlackPix);
      XFillRectangle(display,pix,gc,0,0,WinXSize,InYSize);
      if (SizeChanged) ClearTracer();
      TracerPlot(w);
   }
   if (PicType == MAPTRA)
   {
      if (SizeChanged) ClearTracer();
      if (MapHR.X)
      {
         if (MapPro[w] == AZIMUTHAL)
         {
            if (RedrawFlag[w] || MapLR[w].w != WinXSize || MapLR[w].h != InYSize ||
               (nstep % rmui == 0 && MapLR[w].f != 0))
            {
               RedrawFlag[w] = 0;
               MapLR[w].w    = WinXSize;
               MapLR[w].h    = InYSize;
               MapLR[w].l   += MapLR[w].r;
               if (MapLR[w].l < -180.0) MapLR[w].l += 360.0;
               if (MapLR[w].l >  180.0) MapLR[w].l -= 360.0;
               AzimuthalImage(&MapHR,&MapLR[w]);
            }
         }
         else
         {
            if (RedrawFlag[w] || MapLR[w].w != WinXSize || MapLR[w].h != InYSize)
            {
               RedrawFlag[w] = 0;
               MapLR[w].w  = WinXSize;
               MapLR[w].h  = InYSize;
               if (MapPro[w] == 1) PolarImage(&MapHR,&MapLR[w]);
               else                ScaleImage(&MapHR,&MapLR[w]);
            }
         }
         XPutImage(display,pix,gc,MapLR[w].X,0,0,0,0,MapLR[w].w,MapLR[w].h);
      }
      else
      {
         XSetForeground(display,gc,BlackPix);
         XFillRectangle(display,pix,gc,0,0,WinXSize,InYSize);
      }
      TracerPlot(w);
   }
   if (PicType == ISOCS || PicType == ISOHOR || PicType == ISOLON || PicType == ISOCOL)
   {
      XSetForeground(display,gc,BlackPix);
      XFillRectangle(display,pix,gc,0,0,WinXSize,InYSize);
      IsoAreas(Cstrip);
      IsoLines(Cstrip,0);
   }
   if (PicType == MAPHOR)
   {
      if (MapHR.X)
      {
         if (MapPro[w] == AZIMUTHAL)
         {
            if (RedrawFlag[w] || MapLR[w].w != WinXSize || MapLR[w].h != InYSize ||
               (nstep % rmui == 0 && MapLR[w].f != 0))
            {
               RedrawFlag[w] = 0;
               MapLR[w].w    = WinXSize;
               MapLR[w].h    = InYSize;
               MapLR[w].l   += MapLR[w].r;
               if (MapLR[w].l < -180.0) MapLR[w].l += 360.0;
               if (MapLR[w].l >  180.0) MapLR[w].l -= 360.0;
               AzimuthalImage(&MapHR,&MapLR[w]);
            }
            XPutImage(display,pix,gc,MapLR[w].X,0,0,0,0,MapLR[w].w,MapLR[w].h);
            MapLines(Cstrip,1);
         }
         else
         {
            if (RedrawFlag[w] || MapLR[w].w != WinXSize || MapLR[w].h != InYSize)
            {
               RedrawFlag[w] = 0;
               MapLR[w].w  = WinXSize;
               MapLR[w].h  = InYSize;
               if (MapPro[w] == 1) PolarImage(&MapHR,&MapLR[w]);
               else                ScaleImage(&MapHR,&MapLR[w]);
            }
            XPutImage(display,pix,gc,MapLR[w].X,0,0,0,0,MapLR[w].w,MapLR[w].h);
            IsoLines(Cstrip,1);
         }
      }
      else
      {
         XSetForeground(display,gc,BlackPix);
         XFillRectangle(display,pix,gc,0,0,WinXSize,InYSize);
         IsoLines(Cstrip,1);
      }
   }
   if (PicType == ISOHOV)
   {
      IsoAreas(Cstrip);
   }
   
   if (PicType == ISOSH)
   {
      AmplitudePlot();     
   }

   if (Grid && PicType == ISOCS) ShowGridCS();
   if (Grid && PicType == ISOLON) ShowGridLonsi();
   if (Grid && PicType == ISOHOV && MapPro[w] == 0) ShowGridHov();
   if (Grid && PicType == ISOHOV && MapPro[w] == 1) ShowGridHovT();
   if (PicType == ISOHOR &&  MapPro[w] == 1) ShowGridPolar();
   if (Grid && PicType == ISOHOR && MapPro[w] == 0) ShowGridCyl();
   if (Grid && PicType == ISOCOL) ShowGridCol();

   XCopyArea(display,pix,Win[w],gc,0,0,WinXSize,WinYSize,0,0);
}


/* ==================================================================== */
/* iguistep - this function is called from the model for every timestep */
/* ==================================================================== */

int iguistep_(float pparcs[],int kdatim[])
{
   int i,j,w;               // Loop indices
   struct timeval TimeVal;  // Retrieve time info

   if (Debug) printf("iguistep(%12.2e,%d-%d-%d)\n",
                     pparcs[0],kdatim[0],kdatim[1],kdatim[2]);
   nstep++;
   memcpy(ndatim,kdatim,sizeof(ndatim));
   DeltaTime = ndatim[4] - LastMinute;
   if (DeltaTime < 0) DeltaTime += 60;
   if (DeltaTime ==0) DeltaTime  = 60;
   LastMinute = ndatim[4];

   // Compute frames per second every 20 calls

   gettimeofday(&TimeVal,NULL);
   LastSecond = ThisSecond;
   ThisSecond = TimeVal.tv_sec;

   if ((SecEvent = ThisSecond) > LastSecond)
   {
      fps      = nstep - LastStep;
      LastStep = nstep;
      rmui     = fps / rmuf;
      if (rmui < 1) rmui = 1;
   }

   SkipFreq = 1 + fps / 30; // Reduce plot rate on fast cpu's
   if (SkipFreq < 0 || SkipFreq > 10) SkipFreq = 0;

   for (w=0 ; w < NumWin ; ++w)
   {
      for (j=0 ; j < NumArrays ; ++j)
      {
         if (!strcmp(WinAtt[w].array_name,Array[j].Name))
         {
            if (Array[j].Flag || RedrawFlag[w])
            iso(w,WinAtt[w].Plot,Array[j].Data,Array[j].DimX,Array[j].DimY,
			         Array[j].DimZ,WinAtt[w].Palette);
         }
      }
   }

   for (i=0 ; i < Parcs ; ++i) Parc[i].Val = pparcs[i];
   ShowStep();
   ShowParcs();
   HandleEvents();
   if (XCheckTypedWindowEvent(display,Cow,Expose,&CowEvent))
      RedrawControlWindow();
   for (i=0 ; i < Parcs ; ++i) pparcs[i] = Parc[i].Val;

   XSync(display,1);
   if (Shutdown && MRpid < 1) SaveConfig();
   if (Debug) printf("iguistep returns %d\n",Shutdown);
   return Shutdown;
}

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


void INITGUI(int *model, int *debug, int *lats, int *mrpid, int *mrnum)
{
   initgui_(model,debug,lats,mrpid,mrnum);
}

void GUICLOSE(void)
{
   guiclose_();
}

int IGUISTEP(float pparcs[], int kdatim[])
{
   return iguistep_(pparcs,kdatim);
}

void GUIPUT(char *aname, float *array, int *dimx, int *dimy, int *dimz)
{
   guiput_(aname, array, dimx, dimy, dimz);
}

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


void initgui(int *model, int *debug, int *lats, int *mrpid, int *mrnum)
{
   initgui_(model,debug,lats,mrpid,mrnum);
}

void guiclose(void)
{
   guiclose_();
}

int iguistep(float pparcs[], int kdatim[])
{
   return iguistep_(pparcs,kdatim);
}

void guiput(char *aname, float *array, int *dimx, int *dimy, int *dimz)
{
   guiput_(aname, array, dimx, dimy, dimz);
}

int iguinan_(double *p)
{
   return isnan(*p);
}

int nresources(double *ut, double *st, long *mem, long *par, long *paf,
              long *swa, long *dr, long *dw)
{
   return nresources_(ut,st,mem,par,paf,swa,dr,dw);
}
