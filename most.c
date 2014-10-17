/*

MOST - Model Starter - Version 17 - Edilbert Kirk
---------------------------------------------------------------------------
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
---------------------------------------------------------------------------
*/

#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>
#include <X11/keysym.h>
#include <X11/XKBlib.h>

#define INT int

#define LINEMAX 256
typedef char String[LINEMAX];
String Buffer;

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#define FULLARC (360 * 64)


#define WINDOW_WIDTH 820

/* Models */

#define PUMA   0
#define SAM    1
#define PLASIM 2
#define MODELS 3

/* Resolutions */

#define RES_T15 0
#define RES_T21 1
#define RES_T31 2
#define RES_T42 3
#define RES_T63 4
#define RESOLUTIONS 5

int PlasimSteps[RESOLUTIONS] =
{
   60,  // T15
   45,  // T21
   30,  // T31
   20,  // T42
   10   // T63
};

int PumaSteps[RESOLUTIONS] =
{
   60,  // T15
   60,  // T21
   45,  // T31
   30,  // T42
   15   // T63
};

int ResLat[RESOLUTIONS] =
{
   24,
   32,
   48,
   64,
   96
};

int Resolution = -1;

char *ShortModelName[MODELS] =
{
   "puma",
   "sam",
   "plasim"
};

char *FullModelName[MODELS] =
{
   "PUMA",
   "SAM",
   "Planet Simulator"
};

int Model = PUMA;

#define ALL    -1
#define EARTH   0
#define MARS    1
#define EXO     2
#define PLANETS 3

char *PlanetName[PLANETS] =
{
   "Earth",
   "Mars",
   "Exo"
};

int Planet = EARTH;

double *OroEarth;  // Earth Orography T42 (128 x 64)
double *OroMars;   // Mars  Orography T42 (128 x 64)
double *OroPrep;   // Preprocessed Orography

Pixmap OpmEarth;
Pixmap OpmMars;
Pixmap OpmPrep;

int OroSize = 129 * 64 * sizeof(double);;

double RevGra = 1.0 / 9.81;

/* Object types */

#define SEL_TEXT   1
#define SEL_CHECK  2
#define SEL_INT    3
#define SEL_REAL   4
#define SEL_TEVA   5
#define SEL_PLANET 6

struct ItemStruct
{
   char   list[80];      // Name of namelist
   char   name[80];      // Name of variable
   char   text[80];      // Text    value
   float  pvec[PLANETS]; // Planetary parameter
   double rval;          // Real    value
   int    ival;          // Integer value
   int    model;         // Model
   int    flag;          // Flag
};

struct SelStruct
{
   struct SelStruct  *Next; // Link to next Sel
   struct SelStruct  *Prev; // Link to previous Sel
   struct ItemStruct *Item; // Link to item
   int type;                // Type (TEXT,CHECK,INT,REAL)
   int x;                   // x coordinate of box
   int y;                   // y coordinate of box
   int h;                   // height of box
   int w;                   // width of box
   int xo;                  // x offset
   int xt;                  // x coordinate for text
   int yt;                  // y coordinate for text
   int lt;                  // length of text
   int iv;                  // integer value of box content
   int div;                 // default integer value
   int no;                  // 1: not selectable
   int teco;                // text colour
   int edco;                // edit column (cursor)
   int hide;                // 1: don't show
   float fv;                // floating value of box content
   float dfv;               // default floating value
   float fpl[PLANETS];      // planet dependent floats
   int   *piv;              // pointer to linked integer variable
   float *pfv;              // pointer to linked float   variable
   char text[80];           // text to display
   char teva[16];           // box content
};

struct SelStruct SelStart;
struct SelStruct SelModels[MODELS];
struct SelStruct *CursorSel;
struct SelStruct *ComEnd;
struct SelStruct *SelAno;
struct SelStruct *SelTgr;
struct SelStruct *SelOro;
struct SelStruct *SelMod;
struct SelStruct *SelRes; 
struct SelStruct *SelOce;
struct SelStruct *SelLsg;
struct SelStruct *SelIce;
struct SelStruct *SelAnn;
struct SelStruct *SelPlanet[PLANETS];
struct SelStruct *SelSYear;
struct SelStruct *SelCPU;
struct SelStruct *SelMulti;
struct SelStruct *SelLat2;

#define DIMLOGO 3

struct SymbolStruct 
{
   int     x;     // x pos
   int     y;     // y pos
   int     w;     // width
   int     h;     // height
   int     b;     // background color
   int     f;     // foreground color
   char  **i;     // image
   char   *d;     // image data
   XImage *X;     // X-Image
   char   t[2][80]; // caption
} Logo[DIMLOGO];
  
int Logos;

struct ColorStrip
{
   double  Lo;
   double  Hi;
   char *Name;
   unsigned long pixel;
};

unsigned long WinBG = 0; // Black background
unsigned long TextC = 0; // Text color
unsigned long HeadC = 0; // Heading color

#define DIMFRAME 3

int Frames;

struct FrameStruct 
{
   void (*Action)();
   Pixmap pixmap;   // pixmap
   int    x;        // x pos
   int    y;        // y pos
   int    w;        // width
   int    h;        // height
   int    b;        // background color
   int    f;        // foreground color
   int    xs;       // x pos  selection
   int    ys;       // y pos  selection
   int    ws;       // width  selection
   int    hs;       // height selection
   char   t[3][80]; // text
} Frame[DIMFRAME];
  
int FrameNo;
int Buttons;
int Button1Down = 1;

#define DIMBUTTON 5

struct ButtonStruct 
{
   void (*Action)();
   int    x;        // x pos
   int    y;        // y pos
   int    w;        // width
   int    h;        // height
   int    b;        // background color
   int    f;        // foreground color
   char   t[3][80]; // text
} Button[DIMBUTTON];
  
int Buttons;


/* Coordinates */

int opbox_x;
int opbox_y;
int opbox_w;
int opbox_h;

int nlbox_x;
int nlbox_y;
int nlbox_w;
int nlbox_h;

int nlpos_x;
int nledi_x;

int sfbox_x;
int sfbox_y;
int sfbox_w;
int sfbox_h;
int sfbox_b;

int bubox_x;
int bubox_y;
int bubox_w;
int bubox_h;

int CursorCol;
int LastSel;
int EdiKeyCode;

#define NL_MAX_ITEMS 100

int NL_items = 0;
struct ItemStruct NL_list[NL_MAX_ITEMS];

int BigEndian;
int Oce;
int Ice;
int Lsg;
int SimStart;
int SimYears;
int nreadsr;
int noutput;
int ndebug;
int nprec;
int ngui;
int noro = 1;
int OroAno;
int OroClear;
int OroAqua;
int TgrAno;
int nac;
int Preprocessed;
int SAMindex;
int ScreenHeight;
int Expert = 1;
int LsgEnabled;
int ModeRadiusSq;
int ForceRebuild;
int dxsh;
int dxs2;
int Yoden;

/* Special parameter */

int Latitudes       = 32; // Number of latitudes in atmosphere
int Latitude2       = 32; // Number of latitudes for 2nd. instance
int Levels          = 10; // Number of levels in atmosphere
int Cores           =  1; // Number of cores for parallel version
int Multirun        =  1; // Number of coupled runs
int MultirunEnabled =  0; // Multirun module present?
int Truncation      = 21; // Spectral truncation computed from Latitudes

unsigned Seed;

char cfg_file[256]      = "most_last_used.cfg";
char exec_name[256]     = "most_puma.x";
char exec_nam2[256]     = "most_puma.x";  // 2nd. instance
char exec_ppp[256]      = "most_ppp.x";
char oro_name[256]      = "N064_surf_0129.sra"; // Orography in T42 resolution
char namelist_name[256] = "puma_namelist";
char diag_name[256]     = "puma_diag";
char outp_name[256]     = "puma_output";
char build_name[256]    = "most_puma_build";
char build_ppp[256]     = "most_ppp_build";
char run_name[256]      = "most_puma_run";
char run_ppp[256]       = "most_ppp_run";
char res_name[256]      = "resmod.f90";
char hostname[256]      = "localhost";
char mpirun[256]        = "mpirun";       // "openmpirun" for Open MPI

int ScreenN;
time_t CurrentDate;

Display *display;
unsigned int ScreenW,ScreenH;
int SmallScreen = 1;

static char *progname = "Most";

int DimX;
int DimY;
int DimTr = 21;
int DimSH = 11 * 23;
int DimSE = 11 * 23 + 2;
int ScreenD;
int OffX;
int OffY = 0;
int WinXSize;
int WinYSize;
int WhitePix;
int BlackPix;

int *Ampli;
int *ModeX;
int *ModeY;
int *ModeM;
int *ModeN;

Colormap colormap;

XColor xcolor1,xcolor2;
XColor Red,Green,Blue,Grey,LightRed,DarkRed,LightBlue,DarkBlue;
XColor LightGreen,DarkGreen,Yellow,Cyan,Dummy;

Window Cow; // Control bar
XTextProperty WinconName1;

XEvent WinEvent;

XSizeHints CowSizeHints;
int count;
XEvent report;
GC gc;
char ModFontName[80] = "-misc-fixed-bold-r-normal--15-*-*-*-*-*-*-*";
char FixFontName[80] = "-misc-fixed-bold-r-normal--15-*-*-*-*-*-*-*";
char BigFontName[80] = "-misc-fixed-medium-r-normal--20-*-*-*-*-*-*-*";
//char ModFontName[80] = "9x15bold";
//char FixFontName[80] = "9x15bold";
//char BigFontName[80] = "10x20";
XFontStruct *ModFont;
XFontStruct *FixFont;
XFontStruct *BigFont;
XFontStruct *SamFont;
int ModFontHeight;
int FixFontHeight;
int BigFontHeight;
int ModFontWidth;
int FixFontWidth;
int BigFontWidth;
int ModFontAscent;
int FixFontAscent;
int BigFontAscent;
int EdiFirstKey;
int EdiLastKey;
int EdiSymsPerKey;
int NumLockMask;
int ModeSwitchMask;

int Debug = 0;

KeySym *EdiKeymap;
XModifierKeymap *Mok;

char *display_name = NULL;
XWMHints wm_hints;
XClassHint class_hints;
Atom Delwin;

char *mona[12] =
{
   "Jan","Feb","Mar",
   "Apr","May","Jun",
   "Jul","Aug","Sep",
   "Oct","Nov","Dec"
};

char datch[32];

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

struct MapImageStruct MapHRE; // Hires (2560x1280) Earth image
struct MapImageStruct MapHRM; // Hires (         ) Mars  image
struct MapImageStruct MapLRE; // Lores azimuthal   Earth image
struct MapImageStruct MapLRM; // Lores azimuthal   Mars  image
struct MapImageStruct MapLRK; // Lores azimuthal   Kepler-16
struct MapImageStruct MapLRL; // Lores azimuthal   Kepler-186b

/* Isoarea data */

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

int nlon = 128;
int nlat =  64;
int nlev =   5;

int *Flag;

double VGAX;
double VGAY;

Pixmap pix;

/* Orography units are [m2/s2] */

struct ColorStrip OroStrip[] =
{
   {-99999.0,   80.0,"DarkBlue"},
   {    80.0,  800.0,"green"},
   {   800.0, 1600.0,"brown"},
   {  1600.0, 2400.0,"orange"},
   {  2400.0,99999.0,"white"},
   {     0.0,     0.0,NULL}
};

struct ColorStrip OroMarsStrip[] =
{
   {-99999.0,-16000.0,"yellow"},
   {-16000.0,-12000.0,"RoyalBlue4"},
   {-12000.0, -8000.0,"RoyalBlue3"},
   { -8000.0, -4000.0,"RoyalBlue2"},
   { -4000.0,     0.0,"RoyalBlue1"},
   {     0.0,  4000.0,"dark green"},
   {  4000.0,  8000.0,"light green"},
   {  8000.0, 12000.0,"brown"},
   { 12000.0, 16000.0,"orange"},
   { 16000.0, 99999.0,"white"},
   {     0.0,     0.0,NULL}
};

struct ColorStrip GibbStrip[] =
{
   {-99999.0,    0.0,"DarkBlue"},
   {     0.0,99999.0,"green"},
   {     0.0,     0.0,NULL}
};


struct ColorStrip TStrip[] =
{
   {-99.0,-15.0,"RoyalBlue4"},
   {-15.0,-20.0,"RoyalBlue3"},
   {-10.0, -5.0,"RoyalBlue2"},
   { -5.0, -0.1,"RoyalBlue1"},
   { -0.1,  0.1,"yellow"},
   {  0.1,  5.0,"IndianRed1"},
   {  5.0, 10.0,"IndianRed2"},
   { 10.0, 15.0,"IndianRed3"},
   { 15.0, 20.0,"IndianRed4"},
   { 20.0,999.0,"red"},
   {  0.0,  0.0,NULL}
};

void DoNothing(void)
{
   printf("Do nothing\n");
}


void ChangeModel(int NewMo)
{
   int i;
   struct SelStruct *Sel;

   if (NewMo < 0) // locate actice model
   {
      NewMo = PUMA; // default if none is specified
      for (i=0 , Sel = SelMod ; i < MODELS; ++i , Sel = Sel->Next)
      {
         if (Sel->iv == 1) NewMo = i;
      }
   }

   ComEnd->Next = &SelModels[NewMo];
   ComEnd->Next->Prev = ComEnd;

   for (i=0 , Sel = SelMod ; i < MODELS; ++i , Sel = Sel->Next)
   {
      if (i == NewMo) Sel->iv = 1;
      else            Sel->iv = 0;
   }
   if (NewMo == PUMA || NewMo == SAM)
   {
       SelAno->hide = 0;
       SelAnn->no   = 0;
       SelOce->no   = 1;
       SelIce->no   = 1;
       SelOro->no   = 0;
       SelSYear->no = 1;
       if (SelLsg) SelLsg->no   = 1;
       for (i=0 ; i < PLANETS ; ++i) SelPlanet[i]->no = 1;
   }
   if (NewMo == PLASIM)
   {
       SelAno->hide = 1;
       SelAnn->no   = 1;
       SelOro->no   = 1;
       SelSYear->no = 0;
       if (Expert)
       {
          SelOce->no   = 0;
          SelIce->no   = 0;
          if (SelLsg) SelLsg->no = 0;
          for (i=0 ; i < PLANETS ; ++i) SelPlanet[i]->no = 0;
       }
   }
   Model = NewMo;
}


void WriteSettings(void)
{
   int   mxtl;
   char  tb[80];
   FILE *fp;
   struct SelStruct *Sel;

   fp = fopen("most_last_used.cfg","w");
   if (!fp) return;

   fprintf(fp,"[MoSt 17 - configuration file] created %s\n",ctime(&CurrentDate));

   mxtl = 0;
   for (Sel = &SelStart ; Sel ; Sel = Sel->Next)
   {
      if (Sel->type != SEL_TEXT ) mxtl = MAX(mxtl,strlen(Sel->text));
   }

   for (Sel = &SelStart ; Sel ; Sel = Sel->Next)
   {
      if (Sel->type == SEL_TEXT ) fprintf(fp,"\n[%s]\n",Sel->text);
      else
      {
         memset(tb,' ',sizeof(tb));
         strncpy(tb,Sel->text,strlen(Sel->text));
         tb[mxtl] = 0;
         fprintf(fp,"%s = ",tb);
      }
           if (Sel->type == SEL_CHECK) fprintf(fp,"%6d\n",Sel->iv);
      else if (Sel->type == SEL_INT  ) fprintf(fp,"%6d\n",Sel->iv);
      else if (Sel->type == SEL_TEVA ) fprintf(fp,"%s\n",Sel->teva);
      else if (Sel->type == SEL_REAL || Sel->type == SEL_PLANET)
      {
         if (Sel->type == SEL_PLANET) Sel->fv = Sel->fpl[Planet];
         if (fabs(Sel->fv) >= 1.0 && fabs(Sel->fv) < 9000.0)  fprintf(fp,"%10.4f\n",Sel->fv);
         else                                                 fprintf(fp,"%10.4e\n",Sel->fv);
      }
   }
   fclose(fp);
}


void FormatReal(float fv, char *text)
{
   if (fv == 0.0                      ||
      (fv >   0.001 && fv < 99999.0)  ||
      (fv > -9999.0 && fv <  -0.001))
   {
      sprintf(text,"%11.4f",fv);
           if (!strcmp(text+8,"000")) text[ 8] = 0;
      else if (!strcmp(text+9,"00" )) text[ 9] = 0;
      else if (text[10] == '0')       text[10] = 0;
   }
   else
   {
      sprintf(text,"%11.4e",fv);
   }
}


void UpdateSelections(struct SelStruct *Sel)
{
   char text[80];
   int selx;

   selx = 0;

   for (; Sel ; Sel = Sel->Next)
   {
      if (Sel == ComEnd->Next) selx = 1;
      if (selx) Sel->x = nledi_x;
      if (Sel->type == SEL_INT )
      {
         sprintf(text,"%10d",Sel->iv);
         strcpy(Sel->teva,text+10-Sel->edco);
         if (Sel->piv) *Sel->piv = Sel->iv;
      }
      if (Sel->type == SEL_REAL || Sel->type == SEL_PLANET)
      {
         if (Sel->type == SEL_PLANET) Sel->fv = Sel->fpl[Planet];
         FormatReal(Sel->fv,text);
         strcpy(Sel->teva,text);
      }
   }
   if (SelLat2) SelLat2->no = (SelMulti->iv != 2) ;
}


void ChangePlanet(int NewPlanet)
{
   struct SelStruct *Sel;

   SelPlanet[Planet]->iv = 0;
   Planet = NewPlanet;
   SelPlanet[Planet]->iv = 1;
   strcpy(SelModels[PLASIM].text,PlanetName[Planet]);
   SelModels[PLASIM].lt = strlen(PlanetName[Planet]);

   for (Sel = &SelStart ; Sel ; Sel = Sel->Next)
   {
      if (Sel->type == SEL_PLANET)
      {
         Sel->fv = Sel->fpl[Planet];
      }
   }
   UpdateSelections(ComEnd->Next);
}


int ReadSettings(char *fn)
{
   int pl;
   char  tb[80];
   char *eq;
   FILE *fp;
   struct SelStruct *Sel;

   fp = fopen(fn,"r");
   if (!fp)
   {
      printf("\nUsing default configuration\n");
      return 0;
   }

   fgets(tb,sizeof(tb),fp);
   if (strncmp(tb,"[MoSt",5))
   {
      printf("\nFileheader is: %s\n",tb);
      printf("Expected     : [MoSt ...]\n");
      printf("Using default configuration\n");
      fclose(fp);
      return 0;
   }

   while (!feof(fp))
   {
      eq = strchr(tb,'=');
      if (eq)
      {
         for (Sel = &SelStart ; Sel ; Sel = Sel->Next)
         {
            if (!strncmp(tb,Sel->text,strlen(Sel->text)))
            {
               if (!isalnum(tb[strlen(Sel->text)]))
               {
                  if (Sel->type == SEL_REAL || Sel->type == SEL_PLANET)
                  {
                     Sel->fv = atof(eq+1);
                     if (Sel->pfv) *Sel->pfv = Sel->fv;
                     if (Sel->type == SEL_PLANET) Sel->fpl[Planet] = Sel->fv;
                  }
                  else
                  {
                     Sel->iv = atoi(eq+1);
                     if (Sel->piv) *Sel->piv = Sel->iv;
                  }
               }
            }
         }
      }
      else if (!strncmp(tb,"[Modules]",9)) ChangeModel(-1);
      else
      {
         for (pl=0 ; pl < PLANETS ; ++pl)
         {
            if (!strncmp(tb+1,PlanetName[pl],strlen(PlanetName[pl])) && pl != Planet)
               ChangePlanet(pl);
         }
      }
      fgets(tb,sizeof(tb),fp);
   }
   fclose(fp);
   
   return 1;
}


void Exit(void)
{
   int i;

   WriteSettings();
   // XNextEvent(display,&WinEvent);
   // XCloseDisplay(display);
   exit(0);
}


void AbortMessage(char *s)
{
   printf("\n\n%s\n",s);
   exit(1);
}

void Abort(void)
{
   Exit();
}


void InitNextSelection(struct SelStruct *Sel, int dy, char *t)
{
   Sel->type = Sel->Prev->type;
   Sel->y    = Sel->Prev->y + dy;
   Sel->h    = Sel->Prev->h;
   Sel->w    = Sel->Prev->w;
   Sel->xo   = Sel->Prev->xo;
   Sel->yt   = Sel->Prev->yt + dy;
   strcpy(Sel->text,t);
   Sel->lt   = strlen(Sel->text);
   Sel->teco = Sel->Prev->teco;
   Sel->edco = Sel->Prev->edco;
}


void CursorOn(void)
{
   int x,y,w,h;

   if (CursorSel)
   {
      x = CursorSel->x + CursorCol * FixFontWidth;
      y = CursorSel->y;
      w = CursorSel->w;
      h = CursorSel->h;

      XSetForeground(display,gc,Red.pixel);
      XDrawRectangle(display,Cow,gc,x,y,FixFontWidth,h);
      XDrawRectangle(display,Cow,gc,x-1,y-1,FixFontWidth+2,h+2);
   }
}
   

void RemoveBlanks(char *p)
{
   int l;

   l = strlen(p);
   while (l && p[l-1] == ' ') p[--l] = 0;
}

      
struct SelStruct *NewSel(struct SelStruct *OldSel)
{
   struct SelStruct *Sel;

   Sel = calloc(1,sizeof(struct SelStruct));
   OldSel->Next = Sel;
   Sel->Prev = OldSel;
   return Sel;
}


void NL_i(int m, char *list, char *name, int i)
{
   if (NL_items > NL_MAX_ITEMS-2)
   {
      printf("\n*** Error ***\n");
      printf("Number of namelist items exceeds %d \n",NL_MAX_ITEMS);
      printf("Increase NL_MAX_ITEMS in most,c\n");
      exit(1);
   }
   NL_list[NL_items].model  = m;
   strcpy(NL_list[NL_items].list,list);
   strcpy(NL_list[NL_items].name,name);
   NL_list[NL_items].flag = SEL_INT;
   NL_list[NL_items].ival = i;
   NL_list[NL_items].rval = 0.0;
   NL_list[NL_items].text[0] = 0;
   NL_items++;
}


void NL_r(int m, char *list, char *name, double r)
{
   if (NL_items > NL_MAX_ITEMS-2)
   {
      printf("\n*** Error ***\n");
      printf("Number of namelist items exceeds %d \n",NL_MAX_ITEMS);
      printf("Increase NL_MAX_ITEMS in most,c\n");
      exit(1);
   }
   NL_list[NL_items].model = m;
   strcpy(NL_list[NL_items].list,list);
   strcpy(NL_list[NL_items].name,name);
   NL_list[NL_items].flag = SEL_REAL;
   NL_list[NL_items].ival = 0;
   NL_list[NL_items].rval = r;
   NL_list[NL_items].text[0] = 0;
   NL_items++;
}


void NL_p(char *name, float r[])
{
   if (NL_items > NL_MAX_ITEMS-2)
   {
      printf("\n*** Error ***\n");
      printf("Number of namelist items exceeds %d \n",NL_MAX_ITEMS);
      printf("Increase NL_MAX_ITEMS in most,c\n");
      exit(1);
   }
   NL_list[NL_items].model = PLASIM;
   strcpy(NL_list[NL_items].list,"planet");
   strcpy(NL_list[NL_items].name,name);
   NL_list[NL_items].flag = SEL_PLANET;
   NL_list[NL_items].ival = 0;
   NL_list[NL_items].rval = 0.0;
   NL_list[NL_items].text[0] = 0;
   memcpy(NL_list[NL_items].pvec,r,sizeof(float) * PLANETS);
   NL_items++;
}


void NL_t(int m, char *list, char *name, char *t)
{
   if (NL_items > NL_MAX_ITEMS-2)
   {
      printf("\n*** Error ***\n");
      printf("Number of namelist items exceeds %d \n",NL_MAX_ITEMS);
      printf("Increase NL_MAX_ITEMS in most,c\n");
      exit(1);
   }
   NL_list[NL_items].model = m;
   strcpy(NL_list[NL_items].list,list);
   strcpy(NL_list[NL_items].name,name);
   NL_list[NL_items].flag = SEL_TEVA;
   NL_list[NL_items].ival = 0;
   NL_list[NL_items].rval = 0.0;
   strcpy(NL_list[NL_items].text,t);
   NL_items++;
}


//                             Earth       Mars          Exo
float eccen_vec[PLANETS] = {   0.016715,   0.09341233,   0.0};
float mvelp_vec[PLANETS] = { 102.7     , 336.04084   ,   0.0};
float obliq_vec[PLANETS] = {  23.44    ,  25.19      ,   0.0};
float gsol0_vec[PLANETS] = {1367.0     , 595.0       ,1000.0};

void InitNamelist(void)
{
   // Planet Simulator

   NL_p("ECCEN" , eccen_vec);
   NL_p("MVELP" , mvelp_vec);
   NL_p("OBLIQ" , obliq_vec);
   NL_p("GSOL0" , gsol0_vec);

   NL_i(PLASIM,"planet" ,"NFIXORB" ,  0);
   NL_r(PLASIM,"radmod" ,"CO2"     ,360.0);
   NL_i(PLASIM,"plasim" ,"KICK"    ,  1);
   NL_i(PLASIM,"plasim" ,"MPSTEP"  ,  0);
   NL_i(PLASIM,"plasim" ,"NAQUA"   ,  0);
   NL_i(PLASIM,"plasim" ,"NDIAG"   ,  0);
   NL_i(PLASIM,"plasim" ,"NGUIDBG" ,  0);
   NL_i(PLASIM,"plasim" ,"NQSPEC"  ,  1);
   NL_i(PLASIM,"plasim" ,"NVEG"    ,  0);
   NL_i(PLASIM,"plasim" ,"NWPD"    ,  1);
   NL_i(PLASIM,"plasim" ,"NPRINT"  ,  0);
   NL_i(PLASIM,"plasim" ,"NSYNC"   ,  1);
   NL_i(PLASIM,"rainmod","NCLOUDS" ,  1);
   NL_i(PLASIM,"rainmod","NSTORAIN",  0);
   NL_r(PLASIM,"plasim" ,"SYNCSTR", 0.0);

   // SAM

   NL_i(SAM,"sam","KICK"   ,  1);
   NL_i(SAM,"sam","NAFTER" , 24);
   NL_i(SAM,"sam","NDEL"   ,  8);
   NL_i(SAM,"sam","NDIAG"  ,240);
   NL_i(SAM,"sam","NEXP"   ,  4);
   NL_i(SAM,"sam","NGUIDBG",  0);
   NL_i(SAM,"sam","NTSPD"  ,  0);
   NL_r(SAM,"sam","DISP"   ,0.0);
   NL_r(SAM,"sam","ROTSPD" ,1.0);

   // PUMA

   NL_i(PUMA,"puma","KICK"   ,  1);
   NL_i(PUMA,"puma","MPSTEP" ,  0);
   NL_i(PUMA,"puma","NDEL"   ,  6);
   NL_i(PUMA,"puma","NDHEAT" ,  0);
   NL_i(PUMA,"puma","NDIAG"  ,  0);
   NL_i(PUMA,"puma","NEWSR"  ,  0);
   NL_i(PUMA,"puma","NGUIDBG",  0);
   NL_i(PUMA,"puma","NHELSUA",  0);
   NL_i(PUMA,"puma","NSYNC"  ,  0);
   NL_i(PUMA,"puma","NWPD"   ,  1);
   NL_r(PUMA,"puma","DTEP"   , 60);
   NL_r(PUMA,"puma","DTNS"   ,-70.0);
   NL_r(PUMA,"puma","DTROP"  , 12000.0);
   NL_r(PUMA,"puma","DTTRP"  ,  2.0);
   NL_r(PUMA,"puma","SYNCSTR",  0.0);
   NL_r(PUMA,"puma","ROTSPD" ,  1.0);
   NL_r(PUMA,"puma","TGR"    , 288.0);
};


void ReadNamelist(void)
{
}


void NamelistSelector(int model)
{
   int i,yn,n,ml;
   struct SelStruct *Sel;
   int ledi_x,lbox_w,lbox_h;

   yn = opbox_y;
   n  = 0;
   ml = 0;
   Sel = &SelModels[model];

   if (model == PLASIM)
   {
      Sel->type = SEL_TEXT;
      Sel->teco = HeadC;
      Sel->y    = opbox_y;
      Sel->xt   = nlpos_x;
      Sel->yt   = Sel->y + BigFontAscent + 1;
      strcpy(Sel->text,PlanetName[Planet]);
      Sel->lt   = strlen(Sel->text);
      yn += BigFontHeight + 1;
      Sel = NewSel(Sel);
   }

   for (i=0 ; i < NL_items ; ++i)
   if (NL_list[i].model == model)
   {
      Sel->Item = NL_list + i;
      strcpy(Sel->text,NL_list[i].name);
      Sel->lt   = strlen(Sel->text);
      Sel->type = NL_list[i].flag;
      Sel->teco = BlackPix;
      Sel->h    = FixFontHeight + 1;
      Sel->y    = yn;
      Sel->xt   = nlpos_x;
      Sel->yt   = Sel->y + FixFont->ascent + 1;
      if (Sel->lt > ml) ml = Sel->lt;
      if (Sel->lt > ml) ml = Sel->lt;

      if (NL_list[i].flag == SEL_INT)
      {
         Sel->w    =  6 * FixFontWidth + 2;
         Sel->edco =  6;
         Sel->div  = Sel->iv = NL_list[i].ival;
      }
      else if (NL_list[i].flag == SEL_REAL)
      {
         Sel->edco = 11;
         Sel->w    = 11 * FixFontWidth + 2;
         Sel->dfv  = Sel->fv = NL_list[i].rval;
      }
      else if (NL_list[i].flag == SEL_PLANET)
      {
         Sel->edco = 11;
         Sel->w    = 11 * FixFontWidth + 2;
         Sel->dfv  = Sel->fv = NL_list[i].pvec[Planet];
         memcpy(Sel->fpl,NL_list[i].pvec,sizeof(float) * PLANETS);
      }
      else
      {
         Sel->edco = 11;
         Sel->w    = 11 * FixFontWidth + 2;
         strncpy(Sel->teva,NL_list[i].text,15);
      }

      yn += FixFontHeight + 1;
      Sel = NewSel(Sel);
      ++n;
   }
   Sel->Prev->Next = NULL; // Mark end of chain

   /* Compute size of namelist box */

   ledi_x = nlpos_x + (ml + 1) * FixFontWidth;
   lbox_w = (ml + 14) * FixFontWidth;
   lbox_h = BigFontHeight + (n+2) * FixFontHeight;

   if (nledi_x < ledi_x) nledi_x = ledi_x;
   if (nlbox_w < lbox_w) nlbox_w = lbox_w;
   if (nlbox_h < lbox_h) nlbox_h = lbox_h;

   UpdateSelections(ComEnd->Next);
}


void LoadFonts(void)
{
   if (ScreenH < 700)
   {
      strcpy(ModFontName,"5x8");
      strcpy(FixFontName,"6x10");
      strcpy(BigFontName,"7x13");
   }
   else if (ScreenH < 800)
   {
      strcpy(ModFontName,"6x10");
      strcpy(FixFontName,"7x13");
      strcpy(BigFontName,"9x15bold");
   }
   if ((ModFont = XLoadQueryFont(display,ModFontName)) == NULL)
   {
      printf("%s: Cannot open %s font\n",progname,ModFontName);
      exit(-1);
   }
   if ((FixFont = XLoadQueryFont(display,FixFontName)) == NULL)
   {
      printf("%s: Cannot open %s font\n",progname,FixFontName);
      exit(-1);
   }
   if ((BigFont = XLoadQueryFont(display,BigFontName)) == NULL)
   {
      printf("%s: Cannot open %s font\n",progname,BigFontName);
      exit(-1);
   }
   SamFont = XLoadQueryFont(display,"rk24");
   ModFontWidth  = XTextWidth(ModFont,"X",1);
   FixFontWidth  = XTextWidth(FixFont,"X",1);
   BigFontWidth  = XTextWidth(BigFont,"X",1);
   FixFontAscent = FixFont->ascent;
   ModFontAscent = ModFont->ascent;
   BigFontAscent = BigFont->ascent;
   ModFontHeight = ModFont->ascent + ModFont->descent;
   FixFontHeight = FixFont->ascent + FixFont->descent;
   BigFontHeight = BigFont->ascent + BigFont->descent;
   if (Debug)
   {
      printf("ModFont %2d x %2d %s\n",ModFontWidth,ModFontHeight,ModFontName);
      printf("FixFont %2d x %2d %s\n",FixFontWidth,FixFontHeight,FixFontName);
      printf("BigFont %2d x %2d %s\n",BigFontWidth,BigFontHeight,BigFontName);
   }
}

void InitSelections(void)
{
   int i,l,n,dyn,dys,mw,x,dx;
   char text[80];
   struct SelStruct *Sel;
   struct tm *dati;
   FILE *fp;

   opbox_x  = FixFontWidth / 2;

   dyn = FixFontHeight + 1;
   dys = (3 * dyn) / 2;

   // Initialize Sel at anchor SelStart

   Sel = &SelStart;

   // Model

   Sel->type = SEL_TEXT;
   Sel->teco = HeadC;
   Sel->y    = opbox_y;
   Sel->yt   = Sel->y + BigFontAscent + 1;
   strcpy(Sel->text,"Model");
   Sel->lt   = strlen(Sel->text);

   // PUMA

   Sel = NewSel(Sel);
   InitNextSelection(Sel,BigFontHeight,FullModelName[PUMA]);
   Sel->type = SEL_CHECK;
   Sel->teco = BlackPix;
   Sel->h    = FixFontHeight + 1;
   Sel->w    = FixFontHeight + 1;
   Sel->yt   = Sel->y + FixFontAscent + 1;
   Sel->div  = Sel->iv   =  1;
   SelMod = Sel;

   // SAM

   Sel = NewSel(Sel);
   InitNextSelection(Sel,dyn,"SAM");

   // Planet Simulator

   Sel = NewSel(Sel);
   InitNextSelection(Sel,dyn,FullModelName[PLASIM]);

   // Earth

   Sel = NewSel(Sel);
   InitNextSelection(Sel,dyn,"Earth");
   SelPlanet[EARTH] = Sel;
   Sel->div  = Sel->iv   =  1;
   Sel->no = 1;
   Sel->xo = (5 * FixFontWidth) / 2;

   // Mars

   Sel = NewSel(Sel);
   InitNextSelection(Sel,dyn,"Mars");
   SelPlanet[MARS] = Sel;
   Sel->no  = 1;

   // Exo

   Sel = NewSel(Sel);
   InitNextSelection(Sel,dyn,"Exo");
   SelPlanet[EXO] = Sel;
   Sel->no  = 1;

   // Modules

   Sel = NewSel(Sel);
   InitNextSelection(Sel,dys,"Modules");
   Sel->type = SEL_TEXT;
   Sel->teco = HeadC;
   Sel->h    = 0;
   Sel->w    = 0;
   Sel->yt   = Sel->y + BigFontAscent + 1;
   Sel->xo   = 0;
   
   // Mixed Layer Ocean

   Sel = NewSel(Sel);
   InitNextSelection(Sel,BigFontHeight,"ML  Ocean");
   Sel->type = SEL_CHECK;
   Sel->teco = BlackPix;
   Sel->h    = FixFontHeight + 1;
   Sel->w    = FixFontHeight + 1;
   Sel->yt   = Sel->y + FixFontAscent + 1;
   Sel->div  = Sel->iv   =  0;
   Sel->no   = 1;
   SelOce    = Sel;
   Sel->piv  = &Oce;

   // LSG Ocean

   if (LsgEnabled)
   {
      Sel = NewSel(Sel);
      InitNextSelection(Sel,FixFontHeight,"LSG Ocean");
      Sel->yt   = Sel->y + FixFontAscent + 1;
      Sel->div  = Sel->iv   =  0;
      Sel->no   = 1;
      SelLsg    = Sel;
      Sel->piv  = &Lsg;
   }

   // Ice

   Sel = NewSel(Sel);
   InitNextSelection(Sel,FixFontHeight,"Sea Ice");
   Sel->type = SEL_CHECK;
   Sel->teco = BlackPix;
   Sel->h    = FixFontHeight + 1;
   Sel->w    = FixFontHeight + 1;
   Sel->yt   = Sel->y + FixFontAscent + 1;
   Sel->div  = Sel->iv   =  0;
   Sel->no   = 1;
   SelIce    = Sel;
   Sel->piv  = &Ice;

   // Hardware

   Sel = NewSel(Sel);
   InitNextSelection(Sel,dys,"Parallelism");
   Sel->type = SEL_TEXT;
   Sel->teco = HeadC;
   Sel->h    = 0;
   Sel->w    = 0;
   Sel->yt   = Sel->y + BigFontAscent + 1;
   

   // Number of CPUs

   Sel = NewSel(Sel);
   InitNextSelection(Sel,dys,"Cores");
   Sel->h    = FixFontHeight + 1;
   Sel->w    = FixFontHeight + 1;
   Sel->yt   = Sel->y + FixFont->ascent + 1;
   Sel->w    = 4 * FixFontWidth + 2;
   Sel->type = SEL_INT;
   Sel->edco = 4;
   Sel->div  = Sel->iv = Cores;
   Sel->piv  = &Cores;
   SelCPU    = Sel;
   fp = fopen("most_compiler_mpi","r");
   if (fp) fclose(fp);
   else Sel->no = 1;

   // Number of synchronous runs

   if (MultirunEnabled && Expert)
   {
      Sel = NewSel(Sel);
      InitNextSelection(Sel,dys,"Instances");
      Sel->h    = FixFontHeight + 1;
      Sel->w    = FixFontHeight + 1;
      Sel->yt   = Sel->y + FixFont->ascent + 1;
      Sel->w    = 4 * FixFontWidth + 2;
      Sel->type = SEL_INT;
      Sel->edco = 4;
      Sel->div  = Sel->iv = Multirun;
      Sel->piv  = &Multirun;
      SelMulti  = Sel;
      fp = fopen("most_compiler_mpi","r");
      if (fp) fclose(fp);
      else Sel->no = 1;
   }

   // Resolution

   Sel = NewSel(Sel);
   InitNextSelection(Sel,dys,"Resolution");
   Sel->type = SEL_TEXT;
   Sel->h    = 0;
   Sel->w    = 0;
   Sel->yt   = Sel->y + BigFont->ascent + 1;

   // Horizontal resolution

   if (Expert)
   {
      Sel = NewSel(Sel);
      InitNextSelection(Sel,dys,"Latitudes #1");
      Sel->h    = FixFontHeight + 1;
      Sel->yt   = Sel->y + FixFont->ascent + 1;
      Sel->w    = 4 * FixFontWidth + 2;
      Sel->edco = 4;
      Sel->type = SEL_INT;
      Sel->div  = Sel->iv = Latitudes;
      Sel->piv  = &Latitudes;

      Sel = NewSel(Sel);
      InitNextSelection(Sel,dys,"Latitudes #2");
      Sel->div  = Sel->iv = Latitudes;
      Sel->piv  = &Latitude2;
      Sel->no   = 1;
      SelLat2   = Sel;
   }
   else
   {
      Sel = NewSel(Sel);
      InitNextSelection(Sel,dys,"T21   [64x32]");
      Sel->type = SEL_CHECK;
      Sel->h    = FixFontHeight + 1;
      Sel->w    = FixFontHeight + 1;
      Sel->yt   = Sel->y + FixFont->ascent + 1;
      Sel->div  = Sel->iv   =  1;
      SelRes    = Sel;

      Sel = NewSel(Sel);
      InitNextSelection(Sel,dyn,"T31   [96x48]");
      Sel->type = SEL_CHECK;
      Sel->h    = FixFontHeight + 1;
      Sel->w    = FixFontHeight + 1;
      Sel->div  = Sel->iv   =  0;

      Sel = NewSel(Sel);
      InitNextSelection(Sel,dyn,"T42  [128x64]");
      Sel->type = SEL_CHECK;
      Sel->h    = FixFontHeight + 1;
      Sel->w    = FixFontHeight + 1;
      Sel->div  = Sel->iv   =  0;
   }

   // Vertical resolution

   if (Expert)
   {
      Sel = NewSel(Sel);
      InitNextSelection(Sel,dys,"Levels");
      Sel->type = SEL_INT;
      Sel->div  = Sel->iv = Levels;
      Sel->piv  = &Levels;
   }

   // Options

   Sel = NewSel(Sel);
   InitNextSelection(Sel,dys,"Options");
   Sel->type = SEL_TEXT;
   Sel->h    = 0;
   Sel->w    = 0;
   Sel->yt   = Sel->y + BigFont->ascent + 1;

   // Global debug switch

   Sel = NewSel(Sel);
   InitNextSelection(Sel,dys,"Debug mode");
   Sel->type = SEL_CHECK;
   Sel->h    = FixFontHeight + 1;
   Sel->w    = FixFontHeight + 1;
   Sel->yt   = Sel->y + FixFont->ascent + 1;
   Sel->div  = Sel->iv   =  0;
   Sel->piv  = &ndebug;

   // Precision switch

   Sel = NewSel(Sel);
   InitNextSelection(Sel,dyn,"Double Precision");
   Sel->piv  = &nprec;

   // Global output switch

   Sel = NewSel(Sel);
   InitNextSelection(Sel,dys,"Write Output");
   Sel->type = SEL_CHECK;
   Sel->h    = FixFontHeight + 1;
   Sel->w    = FixFontHeight + 1;
   Sel->yt   = Sel->y + FixFont->ascent + 1;
   Sel->div  = Sel->iv   =  0;
   Sel->piv  = &noutput;

   // GUI

   Sel = NewSel(Sel);
   InitNextSelection(Sel,dyn,"Run with GUI");
   Sel->type = SEL_CHECK;
   Sel->h    = FixFontHeight + 1;
   Sel->w    = FixFontHeight + 1;
   Sel->div  = Sel->iv   =  1;
   Sel->piv  = &ngui;

   // Orography

   Sel = NewSel(Sel);
   InitNextSelection(Sel,dyn,"Orography");
   Sel->type = SEL_CHECK;
   Sel->div  = Sel->iv   =  1;
   Sel->piv  = &noro;
   SelOro    = Sel;

   // Annual cycle

   Sel = NewSel(Sel);
   InitNextSelection(Sel,dyn,"Annual cycle");
   Sel->type = SEL_CHECK;
   Sel->div  = Sel->iv   =  1;
   SelAnn = Sel;

   // Experiment

   Sel = NewSel(Sel);
   InitNextSelection(Sel,dys,"Simulation");
   Sel->type = SEL_TEXT;
   Sel->h    = 0;
   Sel->w    = 0;
   Sel->yt   = Sel->y + BigFont->ascent + 1;

   // Simulation Start

   dati = localtime(&CurrentDate);
   Sel = NewSel(Sel);
   InitNextSelection(Sel,dys,"Start year");
   Sel->h    = FixFontHeight + 1;
   Sel->yt   = Sel->y + FixFont->ascent + 1;
   Sel->w    = 4 * FixFontWidth + 2;
   Sel->edco = 4;
   Sel->type = SEL_INT;
   Sel->div  = Sel->iv = 1;
   Sel->piv  = &SimStart;
   SelSYear  = Sel;

   // Years

   Sel = NewSel(Sel);
   InitNextSelection(Sel,dyn,"Years to run");
   Sel->h    = FixFontHeight + 1;
   Sel->yt   = Sel->y + FixFont->ascent + 1;
   Sel->w    = 4 * FixFontWidth + 2;
   Sel->type = SEL_INT;
   Sel->div  = Sel->iv   = 10;
   Sel->piv  = &SimYears;

   opbox_h = Sel->y + dys - opbox_y;

   // Orography edit

   Sel = NewSel(Sel);
   InitNextSelection(Sel,dys,"Change [gpm]");
   Sel->h    = FixFontHeight + 1;
   Sel->yt   = Sel->y + FixFont->ascent + 1;
   Sel->w    = 6 * FixFontWidth + 2;
   Sel->type = SEL_INT;
   Sel->div  = Sel->iv   = 0;
   Sel->piv  = &OroAno;
   SelAno    = Sel;

   // Ground temperature edit

/*
   Sel = NewSel(Sel);
   InitNextSelection(Sel,dys,"Change [K]");
   Sel->h    = FixFontHeight + 1;
   Sel->yt   = Sel->y + FixFont->ascent + 1;
   Sel->w    = 6 * FixFontWidth + 2;
   Sel->type = SEL_INT;
   Sel->div  = Sel->iv   = 0;
   Sel->piv  = &TgrAno;
   SelTgr    = Sel;
*/

   // Mark end of common selection chain

   ComEnd = Sel;

   mw = 0;
   for (Sel = &SelStart ; Sel ; Sel = Sel->Next)
      if (Sel->w > mw) mw = Sel->w;

   for (Sel = &SelStart ; Sel ; Sel = Sel->Next)
   {
      Sel->xt = Sel->xo + FixFontWidth * 2 + mw;
      if (!Sel->x) Sel->x  = Sel->xt - Sel->w - FixFontWidth;
      if (Sel->type == SEL_INT )
      {
         sprintf(text,"%10d",Sel->iv);
         strcpy(Sel->teva,text+10-Sel->edco);
      }
   }

   for (x = 0 , Sel = &SelStart ; Sel ; Sel = Sel->Next)
   {
      if (Sel->type == SEL_TEXT) l = Sel->lt * BigFontWidth;
      else                       l = Sel->lt * FixFontWidth;
      if (x < l) x = l;
   }
   opbox_w = x + mw + 2 * FixFontWidth;

   nlbox_x = opbox_x + opbox_w + FixFontWidth;
   nlbox_y = opbox_y;
   nlpos_x = nlbox_x + FixFontWidth;
}

void CheckMark(int x, int y, int d)
{
   XSetForeground(display,gc,LightGreen.pixel);
   XFillArc(display,Cow,gc,x+1,y+1,d-5,d-5,0,FULLARC);
}


void DisplayTeVa(int x, int y, int w, char *t)
{
   int l;

   l = strlen(t);
   if (l > w / FixFontWidth) l = w / FixFontWidth;
   XSetForeground(display,gc,WhitePix);
   XSetBackground(display,gc,BlackPix);
   if (l) XDrawImageString(display,Cow,gc,x,y,t,l);
}


void ShowSelection(struct SelStruct *Sel)
{
   int x,y,h,w;

   if (Sel->hide) return;
   h = Sel->h;
   w = Sel->w;
   x = Sel->x;
   y = Sel->y;

   if (!Sel->no)
   {
      if (Sel->type == SEL_CHECK)
      {
         XSetForeground(display,gc,TextC);
         XDrawArc(display,Cow,gc,x,y,w-2,h-2,0,FULLARC);
         if (Sel->iv) CheckMark(x,y,w+1);
      }
      else DisplayTeVa(x+1,Sel->yt,w,Sel->teva);
   }

   XSetForeground(display,gc,TextC);
   XSetBackground(display,gc,WinBG);
   if (Sel->type == SEL_TEXT)
   {
      XSetFont(display, gc, BigFont->fid);
      XSetForeground(display,gc,Sel->teco);
      XDrawImageString(display,Cow,gc,Sel->xt,Sel->yt,Sel->text,Sel->lt);
      XSetFont(display, gc, FixFont->fid);
   }
   else
      XDrawImageString(display,Cow,gc,Sel->xt,Sel->yt,Sel->text,Sel->lt);
}

void InitLogo(void)
{
   int n;

   // KlimaCampus

   n = 0;
   Logo[n].x = FixFontWidth;
   Logo[n].y = 8;

   // PUMA

   ++n;
   strcpy(Logo[n].t[0],FullModelName[PUMA]);

   // Plasim

   ++n;
   strcpy(Logo[n].t[0],FullModelName[PLASIM]);

/*
   Logo[n].x = 350;
   Logo[n].y = 8;
   Logo[n].w = DIMX_LOGO_PLASIM;
   Logo[n].h = DIMY_LOGO_PLASIM;
   Logo[n].b = Blue.pixel;
   Logo[n].f = WhitePix;
   strcpy(Logo[n].t[0],"Planet");
   strcpy(Logo[n].t[1],"Simulator");
   Logo[n].i    = pixelplasim;
*/

   Logos = n;
}

void GenerateNames(void)
{
   Truncation = (2 * Latitudes - 1) / 3;
   sprintf(namelist_name,"%s_namelist",ShortModelName[Model]);
   if (Model == PUMA)
   {
      if (Cores < 2) strcpy(exec_name,"most_puma.x");
      else           strcpy(exec_name,"most_puma_mpi.x");
   }
   else if (Model == SAM)
   {
      if (Cores < 2) strcpy(exec_name,"most_sam.x");
      else           strcpy(exec_name,"most_sam_mpi.x");
   }
   else
      sprintf(exec_name,"most_%s_t%d_l%d_p%d.x",
      ShortModelName[Model],Truncation,Levels,Cores);

   sprintf(diag_name,"%s_diag",ShortModelName[Model]);
   sprintf(outp_name,"%s_output",ShortModelName[Model]);
   sprintf(build_name,"most_%s_build",ShortModelName[Model]);
   sprintf(build_ppp,"most_%s_build","ppp");
   sprintf(run_name,"most_%s_run",ShortModelName[Model]);
   if (Latitudes < 1000) sprintf(oro_name,"N%3.3d_surf_0129.sra",Latitudes);
   else                  sprintf(oro_name,"N%d_surf_0129.sra",Latitudes);
}


int WriteRunScript(int model)
{
   int i;
   int porm;
   FILE *fp;
   char command[256];
   char run[256];

   strcpy(exec_nam2,exec_name); // Duplicate exec name

   if (model == PUMA) // Add Latitudes and Levels as arguments
   {
      sprintf(exec_name+strlen(exec_name)," %d %d",Latitudes,Levels);
      sprintf(exec_nam2+strlen(exec_nam2)," %d %d",Latitude2,Levels);
   }

   if (model == SAM) // Add Latitudes as arguments
   {
      sprintf(exec_name+strlen(exec_name)," %d",Latitudes);
      sprintf(exec_nam2+strlen(exec_nam2)," %d",Latitude2);
   }

   // porm is for Parallel OR Multiprocessing

   porm = MAX(Cores,Multirun);

   strcpy(run,ShortModelName[model]);
   strcat(run,"/run/");
   strcat(run,run_name);

   fp = fopen(run,"w");
   if (fp == NULL)
   {
      printf("Could not open file <%s> for writing\n",run);
      return 0; /* Failure */
   }
   fputs("#!/bin/bash\n",fp);
   fprintf(fp,"# run-script generated by Most %s",ctime(&CurrentDate));
   fputs("EXP=MOST    # Name your experiment here\n",fp);
   fprintf(fp,"[ $# == 1 ] && cd $1\n");
   fprintf(fp,"rm -f %s_restart\n",ShortModelName[model]);
   fputs("rm -f Abort_Message\n",fp);
   fputs("YEAR=0\n",fp);
   fprintf(fp,"YEARS=%d\n",SimYears);
   if (Multirun > 1) fprintf(fp,"INSTANCES=%d\n",Multirun);
   if (ngui) fputs("# Remove '#' from 'while' and 'end' lines for restart loop\n",fp);
   if (ngui) fputs("# ",fp); /* deactivate loop for GUI case */
   fputs("while [ $YEAR -lt $YEARS ]\n",fp);
   if (ngui) fputs("# ",fp); /* deactivate loop for GUI case */
   fputs("do\n",fp);
   fputs("   YEAR=`expr $YEAR + 1`\n",fp);
   fputs("   DATANAME=`printf '%s.%03d' $EXP $YEAR`\n",fp);
   fputs("   DIAGNAME=`printf '%s_DIAG.%03d' $EXP $YEAR`\n",fp);
   fputs("   RESTNAME=`printf '%s_REST.%03d' $EXP $YEAR`\n",fp);
   if (porm < 2)
   {
      fprintf(fp,"   ./%s\n",exec_name);
   }
   else 
   {
      if (Multirun > 1) 
      {    
         fprintf(fp,"   %s",mpirun);
         fprintf(fp," -np 1 %s : -np 1 %s\n",exec_name,exec_nam2);
      }    
      else 
      {    
         fprintf(fp,"   %s -np %d %s\n",mpirun,porm,exec_name);
      }    
   }
   fputs("   [ -e Abort_Message ] && exit 1\n",fp);
   if (Multirun > 1)
   {
      fputs("   INST=0\n",fp);
      fputs("   while [ $INST -lt $INSTANCES ]\n   do\n",fp);
      fputs("      EXT=`printf '%02d' $INST`\n",fp);
      fprintf(fp,"      [ -e %s_$EXT ] && mv %s_$EXT ${EXT}_$DATANAME\n",
              outp_name,outp_name);
      fprintf(fp,"      mv %s_$EXT   ${EXT}_$DIAGNAME\n",diag_name);
      fprintf(fp,"      cp %s_status_$EXT %s_restart_$EXT\n",
              ShortModelName[model],ShortModelName[model]);
      fprintf(fp,"      mv %s_status_$EXT ${EXT}_$RESTNAME\n",ShortModelName[model]);
      fputs("      INST=`expr $INST + 1`\n",fp);
      fputs("   done\n",fp);
   }
   else
   {
      fprintf(fp,"   [ -e %s ] && mv %s $DATANAME\n",outp_name,outp_name);
      fprintf(fp,"   [ -e %s ] && mv %s $DIAGNAME\n",diag_name,diag_name);
      fprintf(fp,"   cp %s_status %s_restart\n",ShortModelName[model],ShortModelName[model]);
      fprintf(fp,"   mv %s_status $RESTNAME\n",ShortModelName[model]);
   }
   if (ngui) fputs("# ",fp); /* deactivate loop for GUI case */
   fputs("done\n",fp);
   fclose(fp);
   sprintf(command,"chmod a+x %s",run);
   system(command);
   return 1; /* Success */
}


void BlueMessage(char *mess)
{
   int l,x,y,w,h;

   XSetFont(display, gc, FixFont->fid);
   l = strlen(mess);
   w = FixFontWidth * (l + 4);
   h = FixFontHeight * 3;
   if (w > WinXSize) w = WinXSize;
   XSetForeground(display,gc,Blue.pixel);
   x = (WinXSize - w) / 2;
   y = (WinYSize - h) / 2;
   XFillRectangle(display,Cow,gc,x,y,w,h);
   XSetForeground(display,gc,WhitePix);
   XSetBackground(display,gc,Blue.pixel);
   XDrawRectangle(display,Cow,gc,x-1,y-1,w+1,h+1);
   x += 2 * FixFontWidth;
   y += 2 * FixFontHeight;
   XDrawImageString(display,Cow,gc,x,y,mess,l);
   XSync(display,0);
}


void RedMessage(char *mess)
{
   int l,x,y,w,h;

   XSetFont(display, gc, FixFont->fid);
   l = strlen(mess);
   w = FixFontWidth * (l + 4);
   h = FixFontHeight * 3;
   if (w > WinXSize) w = WinXSize;
   XSetForeground(display,gc,Red.pixel);
   x = (WinXSize - w) / 2;
   y = (WinYSize - h) / 2;
   XFillRectangle(display,Cow,gc,x,y,w,h);
   XSetForeground(display,gc,WhitePix);
   XSetBackground(display,gc,Red.pixel);
   XDrawRectangle(display,Cow,gc,x-1,y-1,w+1,h+1);
   x += 2 * FixFontWidth;
   y += 2 * FixFontHeight;
   XDrawImageString(display,Cow,gc,x,y,mess,l);
   XSync(display,0);
}


int WriteRunPPP(void)
{
   int i;
   FILE *fp;
   char command[256];
   char run[256];

   strcpy(run,ShortModelName[Model]);
   strcat(run,"/run/");
   strcat(run,run_ppp);

   fp = fopen(run,"w");
   if (fp == NULL)
   {
      printf("Could not open file <%s> for writing\n",run);
      return 0; /* Failure */
   }
   fputs("#!/bin/bash\n",fp);
   fprintf(fp,"# run-script generated by Most %s",ctime(&CurrentDate));
   fputs("[ $# == 1 ] && cd $1\n",fp);
   if (OroClear) fprintf(fp,"rm -f %s\n",oro_name);
   else if (!OroAqua)
      fprintf(fp,"[ ! -e %s -a -e ../dat/%s ] && cp -p ../dat/%s .\n",
              oro_name,oro_name,oro_name);
   fprintf(fp,"   ./%s >!ppp.out\n",exec_ppp);
   fclose(fp);
   sprintf(command,"chmod a+x %s",run);
   system(command);
   strcpy(command,run);
   BlueMessage("Running PUMA Preprocessor");
   sprintf(command,"%s/run/%s %s/run",ShortModelName[Model],run_ppp,ShortModelName[Model]);
   system(command);
   OroAno = 0; /* Anomaly applied -> reset it */
   SelAno->iv = 0;
   return 1; /* Success */
}


void WriteResmod(char *res)
{
   FILE *fp;
   int OldLat,OldLev,OldPro;
   char Line[80];

   // Read existing file if there and check for changes

   fp = fopen(res,"r");
   if (fp)
   {
      fgets(Line,80,fp); // header line 1
      fgets(Line,80,fp); // header line 2
      fgets(Line,80,fp);
      OldLat = atoi(Line+27);
      fgets(Line,80,fp);
      OldLev = atoi(Line+27);
      fgets(Line,80,fp);
      OldPro = atoi(Line+27);
      fclose(fp);
      if (OldLat == Latitudes && OldLev == Levels && OldPro == Cores) return;
   }
   fp = fopen(res,"w");
   if (fp == NULL)
   {
      printf("Could not open file <%s> for writing\n",res);
      exit(1); /* Failure */
   }
   fprintf(fp,"      module resmod ! generated by MoSt %s\n",ctime(&CurrentDate));
   fprintf(fp,"      parameter(NLAT_ATM = %d)\n",Latitudes);
   fprintf(fp,"      parameter(NLEV_ATM = %d)\n",Levels);
   fprintf(fp,"      parameter(NPRO_ATM = %d)\n",Cores);
   fprintf(fp,"      end module resmod\n");
   fclose(fp);
}


int Build(int model)
{
   int i,l,x,y,ierr;
   int porm;
   FILE *fp;
   char script_backup[256];
   char command[256];
   char message[256];
   char bld[256];
   char res[256];
   char *shomo;

   porm = MAX(Cores,Multirun);

   shomo = ShortModelName[model];
   strcpy(bld,shomo);
   strcat(bld,"/bld/");
   strcpy(res,bld);
   strcat(bld,build_name);
   strcat(res,res_name);
   strcpy(script_backup,bld);
   strcat(script_backup,".bak");
   rename(bld,script_backup);

   if (model == PLASIM) WriteResmod(res);

   fp = fopen(bld,"w");
   if (fp == NULL)
   {
      printf("Could not open file <%s> for writing\n",bld);
      return 0; /* Failure */
   }
   fputs("#!/bin/bash\n",fp);
   fprintf(fp,"# compile-script generated by Most %s",ctime(&CurrentDate));

   fputs("[ $# == 1 ] && cd $1\n",fp);
   if (ForceRebuild) fputs("rm -f *.x *.o *.mod\n",fp);
   fputs("cp -p ../src/* .\n",fp);
   if (ngui == 0) // Use stub routines and switch off X11 lib
   {
      putenv("GUIMOD=guimod_stub");
      putenv("PUMAX=pumax_stub");
      putenv("GUILIB=");    
   }
   if (Planet == MARS) putenv("PLAMOD=p_mars");
   if (Planet == EXO ) putenv("PLAMOD=p_exo");
   if (Lsg)
   {
      fputs("cp -p ../../lsg/src/lsgmod.f90 .\n",fp);
      putenv("OCEANCOUP=cpl");
   }
   else
   {
      putenv("OCEANCOUP=cpl_stub");
   }
   if (Latitudes < 4) putenv("FFTMOD=fft991mod");
   if (porm > 1)
   {
      fputs("[ ! -e MPI ] && rm -f *.o *.mod *.x\n",fp);
      fputs("touch MPI\n",fp);
      fputs("cat ../../most_compiler_mpi",fp);
      if (Multirun > 1) putenv("MPIMOD=mpimod_multi");
   }
   else
   {
      fputs("[ -e MPI ] && rm -f *.o *.mod *.x MPI\n",fp);
      fputs("cat ../../most_compiler",fp);
   }
   if (ndebug) fputs(" ../../most_debug_options",fp);
   if (nprec)
   {
      fputs(" ../../most_precision_options",fp);
   }
   fprintf(fp," make_%s > makefile\n",shomo);
   fputs("make -e\n",fp);

   fprintf(fp,"[ $? == 0 ] && cp %s.x ../bin/%s\n",shomo,exec_name);

   fclose(fp);
   sprintf(command,"chmod a+x %s",bld);
   system(command);

   sprintf(message,"Building %s - wait a minute!",FullModelName[Model]);
   BlueMessage(message);
   strcat(bld," ");
   strcat(bld,shomo);
   strcat(bld,"/bld");
   if (system(bld))
   {
      RedMessage("Error in build process");
      sleep(5);
      return 0; // error
   }
   sprintf(command,"cp %s/bin/%s %s/run\n",shomo,exec_name,shomo);
   system(command);

   if (model == PUMA || model == PLASIM )
   {
       if (model == PLASIM && Planet == MARS)
       {
          sprintf(command,"cp images/mars.bmp %s/run/map.bmp\n",shomo);
       system(command);
       }
       else if (model == PLASIM && Planet == EARTH)
       {
          sprintf(command,"cp images/earth.bmp %s/run/map.bmp\n",shomo);
          system(command);
       }
   }
   if (model == PLASIM)
   {
      if (Planet == MARS)
         sprintf(command,"cp plasim/dat/T%d_mars/* plasim/run/\n",Truncation);
      else if (Planet == EXO)
         sprintf(command,"cp plasim/dat/T%d_exo/* plasim/run/\n",Truncation);
      else
         sprintf(command,"cp plasim/dat/T%d/* plasim/run/\n",Truncation);
      system(command);
      system(command);
      if (Lsg)
      {    
         sprintf(command,"cp lsg/dat/* %s/run\n",shomo);
         system(command);
         sprintf(command,"cp %s/dat/GUI_LSG.cfg %s/run/GUI.cfg\n",shomo,shomo);
         system(command);
      }
      else
      {
         sprintf(command,"cp %s/dat/GUI.cfg %s/run\n",shomo,shomo);
         system(command);
      }
   }
   else
   {
      sprintf(command,"cp %s/dat/GUI.cfg %s/run\n",shomo,shomo);
      system(command);
   }
   if (Multirun > 1)
   {
      sprintf(command,"cp %s/dat/GUI_0?.cfg %s/run\n",shomo,shomo);
      system(command);
   }
   return 1; /* Success */
}


int BuildPPP(void)
{
   int i,l,x,y;
   FILE *fp;
   char script_backup[256];
   char command[256];
   char message[256];
   char bld[256];
   char res[256];

   strcpy(bld,ShortModelName[Model]);
   strcat(bld,"/bld/");
   strcpy(res,bld);
   strcat(bld,build_ppp);
   strcat(res,res_name);
   strcpy(script_backup,bld);
   strcat(script_backup,".bak");
   rename(bld,script_backup);

   fp = fopen(bld,"w");
   if (fp == NULL)
   {
      printf("Could not open file <%s> for writing\n",bld);
      return 0; /* Failure */
   }
   fputs("#!/bin/bash\n",fp);
   fprintf(fp,"# compile-script generated by Most %s",ctime(&CurrentDate));
   fputs("[ $# == 1 ] && cd $1\n",fp);
   fputs("rm -f *.o *.mod\n",fp);
   fputs("cp ../src/make_ppp .\n",fp);
   fputs("cp ../src/*.f90 .\n",fp);
   fputs("cat ../../most_compiler",fp);
   if (ndebug) fputs(" ../../most_debug_options",fp);
   if (nprec)  fputs(" ../../most_precision_options",fp);
   fprintf(fp," make_ppp > makefile_ppp\n");
   fputs("make -e -f makefile_ppp\n",fp);
   fprintf(fp,"cp ppp.x ../bin/%s\n",exec_ppp);

   fclose(fp);
   sprintf(command,"chmod a+x %s",bld);
   system(command);

   sprintf(command,"%s/bin/%s",ShortModelName[Model],exec_ppp);
   if (!(fp = fopen(command,"r")))
   {
      sprintf(message,"Building preprocessor ppp");
      BlueMessage(message);
      strcat(bld," ");
      strcat(bld,ShortModelName[Model]);
      strcat(bld,"/bld");
      system(bld);
   }
   else fclose(fp); // Executable exists already

   sprintf(command,"cp %s/bin/%s %s/run\n",ShortModelName[Model],exec_ppp,ShortModelName[Model]);
   system(command);
   return 1; /* Success */
}


void FinishLine(void)
{
   char text[80];
   int  OldCores;

   OldCores = Cores;

   if (CursorSel == NULL) return;
   if (CursorSel->type == SEL_INT)
   {
      CursorSel->iv = atoi(CursorSel->teva);
      sprintf(text,"%10d",CursorSel->iv);
      strcpy(CursorSel->teva,text+10-CursorSel->edco);
      if (CursorSel->piv) *CursorSel->piv = CursorSel->iv;
   }
   if (CursorSel->type == SEL_REAL)
   {
      CursorSel->fv = atof(CursorSel->teva);
      FormatReal(CursorSel->fv,text);
      strcpy(CursorSel->teva,text);
   }
   if (CursorSel == SelMulti)  // Enable or disable Lat2
   {
      SelLat2->no = (SelMulti->iv != 2) ;
   }
   if (OldCores == 1 && Cores  > 1) ForceRebuild = 1;
   if (OldCores  > 1 && Cores == 1) ForceRebuild = 1;
   if (Model == PLASIM && OldCores != Cores) ForceRebuild = 1;
}


int CheckPumaNamelist(void)
{
   int i,safe_ntspd;
   int *ntspd;
   int *nyoden = NULL;
   double s;
   struct SelStruct *Sel;
   FILE *fp;

   FinishLine();

   for (Sel = &SelStart ; Sel ; Sel = Sel->Next)
   {
      if (Sel->piv) *Sel->piv = Sel->iv;
      if (!strcmp(Sel->text,"NTSPD"       )) ntspd    = &Sel->iv;
      if (!strcmp(Sel->text,"NYODEN"      )) nyoden   = &Sel->iv;
      if (!strcmp(Sel->text,"Orography"   )) nreadsr  = Sel->iv;
      if (!strcmp(Sel->text,"Annual cycle")) nac      = Sel->iv;
   }

   // Check for resolution defines as T value

   if (Resolution > 0) Latitudes = ResLat[Resolution];

   // Check # of latitudes for correct values (FFT requirements)

/*
        if (Latitudes >= 2048) Latitudes = 2048; // T1365
   else if (Latitudes >= 1024) Latitudes = 1024; // T682
   else if (Latitudes >=  512) Latitudes =  512; // T341
   else if (Latitudes >=  256) Latitudes =  256; // T170
   else if (Latitudes >=  192) Latitudes =  192; // T127
   else if (Latitudes >=  160) Latitudes =  160; // T106
   else if (Latitudes >=  128) Latitudes =  128; // T85
   else if (Latitudes >=   96) Latitudes =   96; // T63
   else if (Latitudes >=   64) Latitudes =   64; // T42
   else if (Latitudes >=   48) Latitudes =   48; // T31
   else if (Latitudes >=   32) Latitudes =   32; // T21
   else if (Latitudes >=   24) Latitudes =   24; // T15
   else                        Latitudes =   16; // T10
*/

   if (Debug) printf("Latitudes = %d\n",Latitudes);

   // Force # of levels to 20 for predefined Yoden profiles

   if (nyoden && *nyoden > 0)
   {
      Levels = 20;
      Yoden  = *nyoden;
   }

   // Make sure # of CPU's is a power of 2 and in the proper range

   if (Cores > Latitudes/2) Cores = Latitudes/2;
   if (Cores <           1) Cores =           1;
   if ((Latitudes/2) % Cores != 0) Cores = 1;


   // Multirun is currently implemented for two runs only

   if (Multirun >= 2)
   {
      Multirun   = SelMulti->iv = 2;
      Cores = SelCPU->iv   = 1;
   }
   return 0; /* Success */
}


int CheckPlasimNamelist(void)
{
   int i,safe_mpstep;
   int *mpstep;
   double s;
   struct SelStruct *Sel;
   FILE *fp;

   FinishLine();

   for (Sel = &SelStart ; Sel ; Sel = Sel->Next)
   {
      if (Sel->piv) *Sel->piv = Sel->iv;
      if (!strcmp(Sel->text,"MPSTEP")) mpstep = &Sel->iv;
      if (!strcmp(Sel->text,"Orography"   )) nreadsr  = Sel->iv;
      if (!strcmp(Sel->text,"Annual cycle")) nac      = Sel->iv;
   }

   // Check for resolution defines as T value

   if (Resolution > 0) Latitudes = ResLat[Resolution];

   // LSG works currently only with T21 PlaSim
   
   if (Lsg)
   {
       Resolution = RES_T21;
       Latitudes = ResLat[RES_T21];
   }

   // Check # of latitudes for correct values (FFT requirements)

        if (Latitudes >=  64) Latitudes =  64; // T42
   else if (Latitudes >=  48) Latitudes =  48; // T31
   else if (Latitudes >=  32) Latitudes =  32; // T21
   else if (Latitudes >=   4) Latitudes =   4; // T2
   else if (Latitudes >=   2) Latitudes =   2; // T1

   // Set mpstep to a safe value

   if (Resolution > 0) *mpstep = PlasimSteps[Resolution];
   else
   {
      s = log(Latitudes / 32.0) / log(2.0);
      safe_mpstep = 45.0 / pow(3.0,s);
      if (*mpstep > safe_mpstep) *mpstep = safe_mpstep;
   }

   // Make sure # of CPU's is a power of 2 and in the proper range

   if (Cores > Latitudes) Cores = Latitudes;
   if (Cores <         1) Cores =         1;

   if (Latitudes % Cores != 0) Cores = 1;

   return 0; /* Success */
}


double *ReadGrid(const char *fn, int klat, double Sc)
{
   int i,j,k,jlev,iteml,gridsize;
   FILE *ifp;
   int head[8];
   double *gp;
   char Line[256];
   char Elem[16];

   nlat = klat;
   nlon = nlat + nlat;
   VGAX = (double)sfbox_w / nlon;
   VGAY = (double)sfbox_h / (nlat-1);
   gridsize = (nlon+1) * nlat * sizeof(double);

   gp = malloc(gridsize);
   ifp = fopen(fn,"r");
   if (!ifp) return NULL;
   fgets(Line,sizeof(Line),ifp);
   for (i = 0 ; i < 8 ; ++i)
   {
      strncpy(Elem,Line+i*10,10);
      Elem[10] = 0;
      head[i] = atoi(Elem);
   }
   if (head[4] == nlon && head[5] == nlat)
   {
      for (j = 0 ; j < nlat ; ++j)
      {
         for (i = 0 ; i < nlon ; i+=8)
         {
            fgets(Line,sizeof(Line),ifp);
            iteml = strlen(Line) / 8;
            for (k=0 ; k < 8 ; ++k)
            {
               strncpy(Elem,Line+k*iteml,iteml);
               Elem[iteml] = 0;
               gp[k+i+j*(nlon+1)] = atof(Elem) * Sc;
            }
         }
         gp[nlon+j*(nlon+1)] = gp[j*(nlon+1)];
      }
   }
   else printf("head[4] = %d   head[5] = %d\n",head[4],head[5]);
   fclose(ifp);
   return gp;
}


void WriteOroNamelist(FILE *fp)
{
   int lon1,lon2,lat1,lat2;
   double w,h;

   w = Frame[0].w;
   h = Frame[0].h;

   // Round to next gridpoint

   lon1 = 1 + (Frame[0].xs * Latitudes*2) / w;
   lon2 = 1 +((Frame[0].xs  + Frame[0].ws) * Latitudes*2) / w;
   lat1 = 1 + (Frame[0].ys * Latitudes) / h;
   lat2 = 1 +((Frame[0].ys  + Frame[0].hs) * Latitudes) / h;

   fprintf(fp," %-8s=%d\n","NSRV"   ,   0);
   fprintf(fp," %-8s=%d\n","NORO"   ,noro);
   fprintf(fp," %-8s=%d\n","OROANO" ,OroAno);
   fprintf(fp," %-8s=%d\n","LON1ORO",lon1);
   fprintf(fp," %-8s=%d\n","LON2ORO",lon2);
   fprintf(fp," %-8s=%d\n","LAT1ORO",lat1);
   fprintf(fp," %-8s=%d\n","LAT2ORO",lat2);
   fprintf(fp," %-8s=%d\n","TGRANO" ,TgrAno);

   lon1 = 1 + (Frame[1].xs * Latitudes*2) / Frame[1].w;
   lon2 = 1 +((Frame[1].xs  + Frame[1].ws) * Latitudes*2) / Frame[1].w;
   lat1 = 1 + (Frame[1].ys * Latitudes) / Frame[1].h;
   lat2 = 1 +((Frame[1].ys  + Frame[1].hs) * Latitudes) / Frame[1].h;

   fprintf(fp," %-8s=%d\n","LON1TGR",lon1);
   fprintf(fp," %-8s=%d\n","LON2TGR",lon2);
   fprintf(fp," %-8s=%d\n","LAT1TGR",lat1);
   fprintf(fp," %-8s=%d\n","LAT2TGR",lat2);
   return; /* Success */
}


void WriteResolutionNamelist(void)
{
   int i;
   char nln[256];
   FILE *fp;

   for (i = 0 ; i < Multirun ; ++i)
   {
      sprintf(nln,"%s/run/resolution_namelist",ShortModelName[Model]);
      if (Multirun > 1) sprintf(nln+strlen(nln),"_%2.2d",i);
      fp = fopen(nln,"w");
      if (fp == NULL)
      {
         printf("Could not open file <%s> for writing\n",nln);
         return; /* Failure */
      }

      fprintf(fp," &RES\n");
      fprintf(fp," NLAT = %d\n" ,Latitudes);
      if (Model != SAM) fprintf(fp," NLEV = %d\n" ,Levels);
      fprintf(fp," /END\n");
      fclose(fp);
   }
}


/* ================= */
/* WritePumaNamelist */
/* ================= */

int WritePumaNamelist(void)  /* also used for model SAM */
{
   int i,j,k,sum,val,imr;
   FILE *fp;
   char backup_name[256];
   char nln[256];
   struct SelStruct *Sel;

   FinishLine();

   // Write file <puma_namelist> or <sam_namelist>

   for (imr = 0 ; imr < Multirun ; ++imr)
   {
      strcpy(nln,ShortModelName[Model]);
      strcat(nln,"/run/");
      strcat(nln,namelist_name);
      if (Multirun > 1) sprintf(nln+strlen(nln),"_%2.2d",imr);
      fp = fopen(nln,"w");
      if (fp == NULL)
      {
         printf("Could not open file <%s> for writing\n",nln);
         return 0; /* Failure */
      }

      fprintf(fp," &%s_nl\n",ShortModelName[Model]);

      if (Preprocessed || noro) nreadsr = 1;
      fprintf(fp," %-8s=%6d\n","NOUTPUT" ,noutput);
      fprintf(fp," %-8s=%6d\n","NGUI"    ,ngui);
      if (nac)
      {
         fprintf(fp," %-8s=%11.4f\n","TAC"  ,360.0);
         fprintf(fp," %-8s=%11.4f\n","PAC"  ,  0.0);
      }
      else
      {
         fprintf(fp," %-8s=%11.4f\n","TAC"  ,  0.0);
         fprintf(fp," %-8s=%11.4f\n","PAC"  ,  0.0);
      }
      if (ngui) fprintf(fp," %-8s=%6d\n","NYEARS",SimYears);
      else      fprintf(fp," %-8s=%6d\n","NYEARS",1);
      fprintf(fp," %-8s=%6d\n","NMONTHS",0);

      for (Sel = ComEnd->Next ; Sel ; Sel = Sel->Next)
      {
         if (Sel->type == SEL_INT)
            fprintf(fp," %-8s=%6d\n",Sel->text,Sel->iv);
         if (Sel->type == SEL_REAL)
            fprintf(fp," %-8s=%s\n",Sel->text,Sel->teva);
      }

      // Check for mode selections

      for (i=0,sum=0 ; i < DimSH ; ++i) sum += Ampli[i];
      if (sum != DimSH) // Some modes are off
      {
         fprintf(fp," NSPECSEL = ");
         i = 0;
         while (i < DimSH)
         {
            val = Ampli[i];
            j = i;
            while (j < DimSH && Ampli[j] == val) ++j;
            k = j - i;
            if (k == 1)    fprintf(fp,"%d",val);
            else           fprintf(fp,"%d*%d",k,val);
            if (j < DimSH) fprintf(fp,",");
            i = j;
         }
         fprintf(fp,"\n");
      }
      fprintf(fp," /END\n");
      fclose(fp);
   }
   return 1; /* Success */
}


/* ================ */
/* WritePPPNamelist */
/* ================ */

int WritePPPNamelist(void)
{
   int i,j,k,sum,val;
   FILE *fp;
   char backup_name[256];
   char nln[256];
   struct SelStruct *Sel;

   FinishLine();

   WriteResolutionNamelist();

   /* 2. Write file <ppp_namelist> */

   strcpy(nln,ShortModelName[Model]);
   strcat(nln,"/run/");
   strcat(nln,"ppp_namelist");
   strcpy(backup_name,nln);
   strcat(backup_name,".bak");
   rename(nln,backup_name);
   fp = fopen(nln,"w");
   if (fp == NULL)
   {
      printf("Could not open file <%s> for writing\n",nln);
      return 0; /* Failure */
   }

   fprintf(fp," &ppp_nl\n");

   if (Preprocessed || noro) nreadsr = 1;
   fprintf(fp," %-8s=%d\n","NOUTPUT" ,noutput);
   fprintf(fp," %-8s=%d\n","NGUI"    ,ngui);
   if (nac)
   {
      fprintf(fp," %-8s=%11.4f\n","TAC"  ,360.0);
      fprintf(fp," %-8s=%11.4f\n","PAC"  ,  0.0);
   }
   else
   {
      fprintf(fp," %-8s=%11.4f\n","TAC"  ,  0.0);
      fprintf(fp," %-8s=%11.4f\n","PAC"  ,  0.0);
   }

   for (Sel = ComEnd->Next ; Sel ; Sel = Sel->Next)
   {
      if (Sel->type == SEL_INT)
         fprintf(fp," %-8s=%6d\n",Sel->text,Sel->iv);
      if (Sel->type == SEL_REAL)
         fprintf(fp," %-8s=%s\n",Sel->text,Sel->teva);
   }

   WriteOroNamelist(fp);
   fprintf(fp," /END\n");
   fclose(fp);
   return 1; /* Success */
}


/* ======== */
/* FillPoly */
/* ======== */

void FillPoly(int n, double Poly[])
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

void IsoArea(INT y, double vl, double vh, double Top[], double Bot[], INT Dim)
{
   INT  f,x,p;
   double xl,xr,yt,yb;
   double Poly[16];

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


/* ======== */
/* IsoAreas */
/* ======== */

void IsoAreas(double *Field, int DimX, int DimY, struct ColorStrip Strip[])
{
   INT  i;
   INT  y;
   double *Top;
   double *Bot;

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


/* ========= */
/* CalcFrame */
/* ========= */

int CalcFrame(int lats)
{
//   int i,l,x,y;
//   double *Grid;
//
//   nlat = lats;
//   nlon = lats * 2;
//   VGAX = (sfbox_w - 1.0) / nlon;
//   VGAY = (sfbox_h - 1.0) / (nlat-1);
//
//   if (!Flag) Flag = malloc((nlon+1) * sizeof(int));
//   if (!Frame[0].Grid) Frame[0].Grid = malloc(nlon * nlat * nlev * sizeof(double));
//   if (!Frame[0].pixmap) Frame[0].pixmap = XCreatePixmap(display,Cow,sfbox_w,sfbox_h,ScreenD);
//   pix = Frame[0].pixmap;
//   Grid = Frame[0].Grid;
//   IsoAreas(Grid,nlon+1,nlat,Frame[0].Strip);
   return 0;
}


Pixmap IsoPixmap(double *Grid, int nlon, int nlat, struct ColorStrip Strip[])
{
   pix = XCreatePixmap(display,Cow,sfbox_w,sfbox_h,ScreenD);
   //XFillRectangle(display,Cow,gc,0,0,Frame[0].w,Frame[0].h);
   Flag = realloc(Flag,(nlon+1) * sizeof(int));
   IsoAreas(Grid,nlon+1,nlat,Strip);
   //XSetForeground(display,gc,Blue.pixel);
   //XDrawRectangle(display,Cow,gc,Frame[0].x-1,Frame[0].y-1,Frame[0].w,Frame[0].h);
   return pix;
}


int PPPCompiled(void)
{
   FILE *fp;
   char fn[128];

   strcpy(fn,"puma/run/");
   strcat(fn,exec_ppp);

   if ((fp = fopen(fn,"r")))
   {
      fclose(fp);
      return 1;
   }
   return 0;
}


void ShowOrography(void)
{
   int i,m,n,k,l,x,y,h,w,lon1,lon2,lat1,lat2,tx,ty,dx;
   int r,len,width,height,xp,yp;
   char Text[80];
   double flon,flat;

   XSetFont(display, gc, FixFont->fid);
   XSetForeground(display,gc,TextC);
   XSetBackground(display,gc,WinBG);

   if (Model == PLASIM)
   {
      if (Planet == MARS)
      {
         XCopyArea(display,OpmMars,Cow,gc,0,0,
                   Frame[0].w,Frame[0].h,Frame[0].x,Frame[0].y);
      }
      else if (Planet == EARTH)
      {
         XCopyArea(display,OpmEarth,Cow,gc,0,0,
                   Frame[0].w,Frame[0].h,Frame[0].x,Frame[0].y);
      }
      return;
   }

   k = 0; // Frame number 0 reserved for topography
   if (noro == 0) // No orography
   {
      XSetForeground(display,gc,Blue.pixel);
      XFillRectangle(display,Cow,gc,Frame[k].x-1,Frame[k].y-1,Frame[k].w,Frame[k].h);
      XSetForeground(display,gc,WhitePix);
      XSetBackground(display,gc,Blue.pixel);
      l = strlen("Aqua Planet Setup");
      tx = Frame[k].x + (Frame[k].w - l * FixFontWidth) / 2;
      ty = Frame[k].y + Frame[k].h / 2 - FixFontHeight;
      XDrawImageString(display,Cow,gc,tx,ty,"Aqua Planet Setup",l);
      l = strlen("To change mark option 'Orography'");
      tx = Frame[k].x + (Frame[k].w - l * FixFontWidth) / 2;
      ty += FixFontHeight;
      XDrawImageString(display,Cow,gc,tx,ty,"To change mark option 'Orography'",l);
      XSetBackground(display,gc,WinBG);
   }
   else
   {
      XCopyArea(display,OpmPrep,Cow,gc,0,0,Frame[0].w,Frame[0].h,Frame[0].x,Frame[0].y);
      XDrawRectangle(display,Cow,gc,Frame[k].x-1,Frame[k].y-1,Frame[k].w,Frame[k].h);
      l = strlen(Frame[k].t[0]);
      tx = Frame[k].x + (Frame[k].w - l * FixFontWidth) / 2;
      ty = Frame[k].y - FixFont->descent - 2;
      XDrawImageString(display,Cow,gc,tx,ty,Frame[k].t[0],l);
   }
   if (noro && (Button1Down || (Frame[k].ws && Frame[k].hs)))
   {
      x = Frame[k].x+Frame[k].xs;
      y = Frame[k].y+Frame[k].ys;
      w = Frame[k].ws;
      h = Frame[k].hs;
      if (w < 0) {x += w; w = -w;}
      if (h < 0) {y += h; h = -h;}
      XSetForeground(display,gc,WhitePix);
      XDrawArc(display,Cow,gc,x,y,w,h,0,FULLARC);
      XSetForeground(display,gc,BlackPix);
      XDrawArc(display,Cow,gc,x-1,y-1,w+2,h+2,0,FULLARC);
      flon = 360.0 / Frame[k].w;
      flat = 180.0 / Frame[k].h;
      lon1 = flon * (x - Frame[k].x);
      lon2 = lon1 + flon * w;
      lat1 = 90.0 - flat * (y - Frame[k].y);
      lat2 = lat1 - flat * h;
      sprintf(Frame[k].t[1],"Lon: (%3d:%3d) Lat: (%3d:%3d)",lon1,lon2,lat2,lat1);
      l = strlen(Frame[k].t[1]);
      tx = Frame[k].x;
      ty = Frame[k].y + Frame[k].h + FixFont->ascent;
      XDrawImageString(display,Cow,gc,tx,ty,Frame[k].t[1],l);
      XSetForeground(display,gc,WhitePix);
      sprintf(Frame[k].t[1],"Mouse marks area             ");
      l = strlen(Frame[k].t[1]);
      tx = Frame[k].x;
      ty = Frame[k].y + Frame[k].h + FixFont->ascent;
      XDrawImageString(display,Cow,gc,tx,ty,Frame[k].t[1],l);
   }
}


void ShowMars(void)
{
   if (Debug) printf("ShowMars\n");
   XPutImage(display,Cow,gc,MapLRM.X,0,0,Frame[1].x,Frame[1].y,MapLRM.w,MapLRM.h);
}


void ShowExo(void)
{
   XPutImage(display,Cow,gc,MapLRL.X,0,0,Frame[0].x,Frame[0].y,MapLRL.w,MapLRL.h);
   XPutImage(display,Cow,gc,MapLRK.X,0,0,Frame[1].x,Frame[1].y,MapLRK.w,MapLRK.h);
}


void ShowEarth(void)
{
   if (Debug) printf("ShowEarth\n");
   XPutImage(display,Cow,gc,MapLRE.X,0,0,Frame[1].x,Frame[1].y,MapLRE.w,MapLRE.h);
}


void ShowModeSelector(void)
{
   int k,dx,r,m,i,n,len,width,height,xp,yp,l,tx,ty;
   char Text[80];

   if (Debug) printf("ShowModeSelector\n");
   k  = 1; // Frame number 1 reserved for mode display
   pix = Frame[k].pixmap;
   strcpy(Frame[1].t[0],"Spherical Harmonics mode selector");
   dx = dxsh;
   r  = (dx+dx) / 3;
   XSetForeground(display,gc,BlackPix);
   XFillRectangle(display,pix,gc,0,0,Frame[k].w,Frame[k].h);
   for (m=0,i=0 ; m <= DimTr ; ++m)
   {
      for (n=m ; n <= DimTr ; ++n,++i)
      {
         if (Ampli[i]) XSetForeground(display,gc,Green.pixel);
         else          XSetForeground(display,gc,  Red.pixel);
         XFillArc(display,pix,gc,ModeX[i]-dxs2,ModeY[i]-dxs2,r,r,0,360*64);
      }
   }

   /* Draw mode legend */
   
   XSetForeground(display,gc,Yellow.pixel);
   XSetBackground(display,gc,BlackPix);
   XSetFont(display, gc, ModFont->fid);
   for (i=0 ; i < 21 ; i+=2)
   {
      sprintf(Text,"%d",i);
      len    = strlen(Text);
      width  = XTextWidth(ModFont,Text,len);
      height = ModFont->ascent + ModFont->descent;
      xp     = dx + i * dx - width/2 - ModFontWidth/2 + 1;
      yp     = ModFontHeight;
      XDrawImageString(display,pix,gc,xp,yp,Text,len);
   }
   strcpy(Text,"PUMA T21 only!");
   len = strlen(Text);
   xp  = dx - ModFontWidth/2 + 1;
   yp  = 2 * ModFontHeight + 10 * dx;
   XDrawImageString(display,pix,gc,xp,yp,Text,len);

   XSetForeground(display,gc,WhitePix);
   strcpy(Text,"MB 1: Toggle mode");
   len = strlen(Text);
   yp += 2*dx;
   XDrawImageString(display,pix,gc,xp,yp,Text,len);

   strcpy(Text,"MB 2: Toggle column");
   len = strlen(Text);
   yp += dx;
   XDrawImageString(display,pix,gc,xp,yp,Text,len);

   strcpy(Text,"MB 3: Toggle line");
   len = strlen(Text);
   yp += dx;
   XDrawImageString(display,pix,gc,xp,yp,Text,len);

   XSetForeground(display,gc,Yellow.pixel);
   strcpy(Text,"n : Total Wavenumber");
   len = strlen(Text);
   xp  = dx - ModFontWidth/2 + 1;
   yp  = 2 * ModFontHeight + 18 * dx;
   XDrawImageString(display,pix,gc,xp,yp,Text,len);

   XSetForeground(display,gc,LightBlue.pixel);
   for (i=0 ; i < 21 ; i+=2)
   {
      sprintf(Text,"%d",i);
      len    = strlen(Text);
      width  = XTextWidth(ModFont,Text,len);
      height = ModFont->ascent + ModFont->descent;
      xp     = Frame[k].w - width - 2;
      yp     = 2 * ModFontHeight + i * dx;
      XDrawImageString(display,pix,gc,xp,yp,Text,len);
   }
   strcpy(Text,"m : Zonal Wavenumber");
   len = strlen(Text);
   xp  = dx - ModFontWidth/2 + 1;
   yp  = 2 * ModFontHeight + 19 * dx;
   XDrawImageString(display,pix,gc,xp,yp,Text,len);

   XSetForeground(display,gc,Green.pixel);
   strcpy(Text,"Mode is on");
   len = strlen(Text);
   xp  = dx - ModFontWidth/2 + 1;
   yp  = ModFontHeight + 16 * dx;
   XFillArc(display,pix,gc,xp,yp+ModFont->descent+2,r,r,0,360*64);
   xp += dx;
   yp += dx;
   XDrawImageString(display,pix,gc,xp,yp,Text,len);

   XSetForeground(display,gc,Red.pixel);
   strcpy(Text,"Mode is off");
   len = strlen(Text);
   xp  = dx - ModFontWidth/2 + 1;
   yp  = ModFontHeight + 17 * dx;
   XFillArc(display,pix,gc,xp,yp+ModFont->descent+2,r,r,0,360*64);
   xp += dx;
   yp += dx;
   XDrawImageString(display,pix,gc,xp,yp,Text,len);

   XSetForeground(display,gc,WhitePix);
   strcpy(Text,"Switch all modes on");
   len = strlen(Text);
   xp  = dx - ModFontWidth/2 + 1;
   yp  = ModFontHeight + 20 * dx;
   XFillArc(display,pix,gc,xp,yp+ModFont->descent+2,r,r,0,360*64);
   xp += dx;
   yp += dx;
   XSetForeground(display,gc,Green.pixel);
   XDrawImageString(display,pix,gc,xp,yp,Text,len);

   XSetForeground(display,gc,WhitePix);
   strcpy(Text,"Switch all modes off");
   len = strlen(Text);
   xp  = dx - ModFontWidth/2 + 1;
   yp  = ModFontHeight + 21 * dx;
   XFillArc(display,pix,gc,xp,yp+ModFont->descent+2,r,r,0,360*64);
   xp += dx;
   yp += dx;
   XSetForeground(display,gc,Red.pixel);
   XDrawImageString(display,pix,gc,xp,yp,Text,len);

   XCopyArea(display,Frame[k].pixmap,Cow,gc,0,0,Frame[k].w,Frame[k].h,Frame[k].x,Frame[k].y);

   XSetForeground(display,gc,Cyan.pixel);
   l = strlen(Frame[k].t[0]);
   tx = Frame[k].x + (Frame[k].w - l * FixFontWidth) / 2;
   ty = Frame[k].y - FixFont->descent - 2;
   XDrawImageString(display,Cow,gc,tx,ty,Frame[k].t[0],l);
   XSetForeground(display,gc,LightBlue.pixel);
   l = strlen(Frame[k].t[1]);
   if (l)
   {
      tx = Frame[k].x + (3 * ModFontWidth) / 2;
      ty = Frame[k].y + Frame[k].h + FixFont->ascent;
      XDrawImageString(display,Cow,gc,tx,ty,Frame[k].t[1],l);
   }
}


void ShowFrame1(void)
{
   if (Model == PLASIM)
   {
           if (Planet == MARS) ShowMars();
      else if (Planet == EXO ) ShowExo();
      else                     ShowEarth();
   }
   else ShowModeSelector();
}


void PreProcess(void)
{
   FinishLine();
   if (Model == PUMA)
   {
      if (CheckPumaNamelist()) return;
      GenerateNames();
      WritePPPNamelist();
      WritePumaNamelist();
      if (!PPPCompiled()) BuildPPP();
      WriteRunPPP();
      Preprocessed = 1;
      if (Latitudes < 1000)
         sprintf(oro_name,"puma/run/N%3.3d_surf_0129.sra",Latitudes);
      else
         sprintf(oro_name,"puma/run/N%d_surf_0129.sra",Latitudes);
      OroPrep = ReadGrid(oro_name,Latitudes,RevGra);
      OpmPrep = IsoPixmap(OroPrep,2*Latitudes,Latitudes,OroStrip);
      free(OroPrep);
      OroClear = 0;
   }
   else
   {
      RedMessage("Preprocessing is for PUMA only");
      sleep(5);
   }
}


void ClearOro(void)
{
   OroClear = 1;
   OroAqua  = 1;
   OroAno   = 0;
   PreProcess();   
}

void WriteNamelistFile(char *nl,  int instance)
{
   char *shomo = ShortModelName[Model];
   struct SelStruct *Sel;
   char fn[256];
   FILE *fp;

   if (Multirun > 1)
      sprintf(fn,"%s/run/%s_namelist_%2.2d",shomo,nl,instance);
   else
      sprintf(fn,"%s/run/%s_namelist",shomo,nl);

   fp = fopen(fn,"w");
   if (fp == NULL)
   {
      printf("Could not open file <%s> for writing\n",fn);
      return; /* Failure */
   }

   // Write namelist file

   fprintf(fp," &%s_nl\n",nl);

   // Write special parameters

   if (!strcmp(nl,"icemod"))
   {
      fprintf(fp," %-12s=%6d\n","NICE",Ice);
   }
   if (!strcmp(nl,"oceanmod"))
   {
      fprintf(fp," %-12s=%6d\n","NOCEAN",Oce);
      fprintf(fp," %-12s=%6d\n","NLSG"  ,Lsg);
   }
   if (!strcmp(nl,"plasim"))
   {
      fprintf(fp," %-12s=%6d\n","NOUTPUT",noutput);
      fprintf(fp," %-12s=%6d\n","NGUI"  ,ngui);
      fprintf(fp," %-12s=%6d\n","N_START_YEAR",SimStart);
      if (Lsg)        fprintf(fp," %-12s=%6d\n","N_DAYS_PER_YEAR",360);
      else if (Planet != MARS) fprintf(fp," %-12s=%6d\n","N_DAYS_PER_YEAR",365);
      if (ngui) fprintf(fp," %-12s=%6d\n","N_RUN_YEARS",SimYears);
      else      fprintf(fp," %-12s=%6d\n","N_RUN_YEARS",1);
      fprintf(fp," %-12s=%6d\n","N_RUN_MONTHS",0);
      fprintf(fp," %-12s=%6d\n","N_RUN_DAYS",0);
   }
   for (Sel = ComEnd->Next ; Sel ; Sel = Sel->Next)
   {
      if (Sel->Item && !strcasecmp(Sel->Item->list,nl))
      {    
         if (Sel->type == SEL_TEVA)
            fprintf(fp," %-12s=%c%s%c\n",Sel->text,'\"',Sel->teva,'\"');
         else
            fprintf(fp," %-12s=%s\n",Sel->text,Sel->teva);
      }    
   }
   fprintf(fp," /END\n");
   fclose(fp);
}


/* =================== */
/* WritePlasimNamelist */
/* =================== */

void WritePlasimNamelist(void)
{
   int imr;

   FinishLine();

   for (imr = 0 ; imr < Multirun ; ++imr)
   {
      WriteNamelistFile("fluxmod" ,imr);
      WriteNamelistFile("icemod"  ,imr);
      WriteNamelistFile("landmod" ,imr);
      WriteNamelistFile("vegmod"  ,imr);
      WriteNamelistFile("miscmod" ,imr);
      WriteNamelistFile("oceanmod",imr);
      WriteNamelistFile("planet"  ,imr);
      WriteNamelistFile("plasim"  ,imr);
      WriteNamelistFile("radmod"  ,imr);
      WriteNamelistFile("rainmod" ,imr);
      WriteNamelistFile("seamod"  ,imr);
      WriteNamelistFile("surfmod" ,imr);
   }
}


void BuildScripts(void)
{
   if (Model == PUMA)
   {
      if (CheckPumaNamelist()) return;
      GenerateNames();
      WritePumaNamelist();
      if (!Build(PUMA)) Exit();
      WriteRunScript(PUMA);
   }
   if (Model == SAM)
   {
      if (CheckPumaNamelist()) return;
      GenerateNames();
      WritePumaNamelist();
      if (!Build(SAM)) Exit();
      WriteRunScript(SAM);
   }
   if (Model == PLASIM)
   {
      CheckPlasimNamelist();
      GenerateNames();
      WritePlasimNamelist();
      if (!Build(PLASIM)) Exit();
      WriteRunScript(PLASIM);
   }
}


void SaveExit(void)
{
   BuildScripts();
   if ((noro || Yoden) && !Preprocessed && Model == PUMA) PreProcess();
   Exit();
}


void SaveRun(void)
{
   char command[256];
   BuildScripts();
   if ((noro || Yoden) && !Preprocessed && Model == PUMA) PreProcess();
   sprintf(command,"%s/run/%s %s/run &",ShortModelName[Model],run_name,ShortModelName[Model]);
   system(command);
   printf("\n=== Success: Launched process %s ===\n\n",run_name);
   Exit();
}


void MarkRectangle(void)
{
   if (Model != PUMA) return;
   if (WinEvent.type == ButtonPress) 
   {
      Frame[FrameNo].xs = WinEvent.xbutton.x - Frame[FrameNo].x;
      Frame[FrameNo].ys = WinEvent.xbutton.y - Frame[FrameNo].y;
      Frame[FrameNo].ws = 0;
      Frame[FrameNo].hs = 0;
   }
   if (WinEvent.type == ButtonRelease) 
   {
      Frame[FrameNo].ws = WinEvent.xbutton.x - Frame[FrameNo].x - Frame[FrameNo].xs;
      Frame[FrameNo].hs = WinEvent.xbutton.y - Frame[FrameNo].y - Frame[FrameNo].ys;
   }
   if (WinEvent.type == MotionNotify && Button1Down) 
   {
      Frame[FrameNo].ws = WinEvent.xbutton.x - Frame[FrameNo].x - Frame[FrameNo].xs;
      Frame[FrameNo].hs = WinEvent.xbutton.y - Frame[FrameNo].y - Frame[FrameNo].ys;
   }
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


/* =============================================== */
/* RectImage - Display map in azimuthal projection */
/* =============================================== */

void RectImage(struct MapImageStruct *s, struct MapImageStruct *d)
{
   int x  ;           // x      pixel coordinate in destination image
   int y  ;           // y      pixel coordinate in destination image

   unsigned int dih;  // destination image height
   unsigned int diw;  // destination image width
   unsigned int dpw;  // destination image padded width

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

   for (y = 0 ; y < dih ; ++y)
   {
      for (x = 0 ; x < diw ; ++x)
      {
         XPutPixel(dX,x,y,XGetPixel(sX,x,y));
      }
   }
}


void ToggleMode(void)
{
   int i,j,m,n,mx,my,dx,dy;

   if (Model == PLASIM) // Rotate image
   {
      if (WinEvent.xbutton.button == Button1) // Rotate right
      {
         if (Planet == MARS)
         {
            MapLRM.l += 10;
            if (MapLRM.l >  180) MapLRM.l -= 360;
            AzimuthalImage(&MapHRM,&MapLRM);
         }
         else
         {
            MapLRE.l += 10;
            if (MapLRE.l >  180) MapLRE.l -= 360;
            AzimuthalImage(&MapHRE,&MapLRE);
         }
      }
      else if (WinEvent.xbutton.button == Button3) // Rotate left
      {
         if (Planet == MARS)
         {
            MapLRM.l -= 10;
            if (MapLRM.l < -180) MapLRM.l += 360;
            AzimuthalImage(&MapHRM,&MapLRM);
         }
         else
         {
            MapLRE.l -= 10;
            if (MapLRE.l < -180) MapLRE.l += 360;
            AzimuthalImage(&MapHRE,&MapLRE);
         }
      }
      return;
   }
   if (WinEvent.type == ButtonPress) 
   {
      mx = WinEvent.xbutton.x - Frame[FrameNo].x;
      my = WinEvent.xbutton.y - Frame[FrameNo].y;
      for (i=0 ; i < DimSE ; ++i)
      {
         dx = mx - ModeX[i];
         dy = my - ModeY[i];
         if ((dx*dx + dy*dy) < ModeRadiusSq) break;
      }
      if (WinEvent.xbutton.button == Button1)
      {
               if (i == DimSH  ) // All modes on
            for (j=0 ; j < DimSH ; ++j) Ampli[j] = 1;
         else if (i == DimSH+1) // All modes off
            for (j=0 ; j < DimSH ; ++j) Ampli[j] = 0;
         else Ampli[i] = 1 - Ampli[i];
      }
      else if (WinEvent.xbutton.button == Button2 && i < DimSH)
      {
         n = ModeN[i];
         for (j=0 ; j < DimSH ; ++j)
            if (ModeN[j] == n) Ampli[j] = 1 - Ampli[j];
      }
      else if (WinEvent.xbutton.button == Button3 && i < DimSH)
      {
         m = ModeM[i];
         for (j=0 ; j < DimSH ; ++j)
            if (ModeM[j] == m) Ampli[j] = 1 - Ampli[j];
      }
   }
   if (WinEvent.type == MotionNotify) 
   {
      mx = WinEvent.xbutton.x - Frame[FrameNo].x;
      my = WinEvent.xbutton.y - Frame[FrameNo].y;
      for (i=0 ; i < DimSH ; ++i)
      {
         dx = mx - ModeX[i];
         dy = my - ModeY[i];
         if ((dx*dx + dy*dy) < ModeRadiusSq)
         {
            sprintf(Frame[1].t[1],"Mode (%2d:%2d)",ModeN[i],ModeM[i]);
            return;
         }
      }
      Frame[1].t[1][0] = 0;
   }
   return;
}


void InitFrames(void)
{
   int i,j,n;
   int moselw;
   FILE *ftp;

   /* Size and position of orography window */

   sfbox_w = 360;
   sfbox_h = 180;
   sfbox_x = WINDOW_WIDTH - sfbox_w - 3 * FixFontWidth;
   sfbox_y = nlbox_y + 4;

   /* Parameter for default orographies */

   nlat = 64;
   nlon = nlat * 2;
   VGAX = (sfbox_w - 1.0) / nlon;
   VGAY = (sfbox_h - 1.0) / (nlat-1);

   /* Read T42 orography for Earth and Mars */

   OroEarth = ReadGrid("plasim/dat/T42/N064_surf_0129.sra",nlat,RevGra);
   OroMars  = ReadGrid("plasim/dat/T42_mars/N064_surf_0129.sra",nlat,1.0/3.74);
   OpmEarth = IsoPixmap(OroEarth,nlon,nlat,OroStrip);
   OpmPrep  = IsoPixmap(OroEarth,nlon,nlat,OroStrip);
   OpmMars  = IsoPixmap(OroMars ,nlon,nlat,OroMarsStrip);
   free(OroEarth);
   free(OroMars);

   /* Check for orography in run directory */

   ftp = fopen(oro_name,"r");
   if (ftp) fclose(ftp);
   else strcpy(oro_name,"puma/dat/N064_surf_0129.sra");

   // Orography

   n = 0;

   Frame[n].x = sfbox_x;
   Frame[n].y = sfbox_y + FixFontHeight;
   Frame[n].w = sfbox_w;
   Frame[n].h = sfbox_h;
   Frame[n].b = DarkGreen.pixel;
   Frame[n].f = WhitePix;
   strcpy(Frame[n].t[0],"Orography");
   strcpy(Frame[n].t[1],"64 x 32");
   strcpy(Frame[n].t[2],"");
   Frame[n].Action = MarkRectangle;
   Frame[n].xs =  110;
   Frame[n].ys =   95;
   Frame[n].ws =   60;
   Frame[n].hs =   40;
   SelAno->x    = Frame[n].x;
   SelAno->y    = Frame[n].y + Frame[n].h + FixFontHeight;
   SelAno->edco = 6;
   SelAno->w    = SelAno->edco * FixFontWidth + 2;
   SelAno->xt   = SelAno->x + (SelAno->edco + 1) * FixFontWidth;
   SelAno->yt   = SelAno->y + FixFont->ascent + 1;
   strcpy(SelAno->teva,"     0");
   Button[4].y = SelAno->y - FixFontHeight/2; // Clear orography button 
   Button[4].x = Frame[n].x + Frame[n].w - 6 * FixFontWidth;

   // T21 mode selector


   ++n;
   Frame[n].x = sfbox_x;
   Frame[n].y = Frame[n-1].y + Frame[n-1].h + 4 * FixFontHeight;
   moselw = sfbox_w;
   if (SmallScreen)
   {
       moselw = WinYSize - Frame[n].y - FixFontHeight - FixFontHeight/2;
       Frame[n].x += (sfbox_x - moselw) / 2;
   }
   Frame[n].w = moselw;
   Frame[n].h = moselw;
   Frame[n].b = BlackPix;
   Frame[n].f = WhitePix;
   strcpy(Frame[n].t[0],"Spherical Harmonics mode selector");
   strcpy(Frame[n].t[1],"");
   strcpy(Frame[n].t[2],"");
   Frame[n].Action = ToggleMode;

   sfbox_b = Frame[n].y + Frame[n].h - Frame[0].y + 2 * FixFontHeight + FixFontHeight/2;

   Frames = n + 1;
}


void InitButtons(void)
{
   int n;

   bubox_w = 64;
   bubox_h = 64;
   bubox_x = WINDOW_WIDTH - 4 * 64 - 6 * FixFontWidth;
   bubox_y = 8;

   // Preprocess

   n = 0;

   Button[n].x = bubox_x;
   Button[n].y = bubox_y;
   Button[n].w = bubox_w;
   Button[n].h = bubox_h;
   Button[n].b = Yellow.pixel;
   Button[n].f = BlackPix;
   strcpy(Button[n].t[0],"Pre-");
   strcpy(Button[n].t[1],"pro-");
   strcpy(Button[n].t[2],"cess");
   Button[n].Action = PreProcess;

   // Save & Exit

   ++n;
   Button[n].x = Button[n-1].x + Button[n-1].w + FixFontWidth;
   Button[n].y = Button[n-1].y;
   Button[n].w = bubox_w;
   Button[n].h = bubox_h;
   Button[n].b = DarkGreen.pixel;
   Button[n].f = WhitePix;
   strcpy(Button[n].t[0],"Save");
   strcpy(Button[n].t[1],"&");
   strcpy(Button[n].t[2],"Exit");
   Button[n].Action = SaveExit;

   // Save & Run

   ++n;
   Button[n].x = Button[n-1].x + Button[n-1].w + FixFontWidth;
   Button[n].y = Button[n-1].y;
   Button[n].w = bubox_w;
   Button[n].h = bubox_h;
   Button[n].b = DarkBlue.pixel;
   Button[n].f = WhitePix;
   strcpy(Button[n].t[0],"Save");
   strcpy(Button[n].t[1],"&");
   strcpy(Button[n].t[2],"Run");
   Button[n].Action = SaveRun;

   // Abort

   ++n;
   Button[n].x = Button[n-1].x + Button[n-1].w + FixFontWidth;
   Button[n].y = Button[n-1].y;
   Button[n].w = bubox_w;
   Button[n].h = bubox_h;
   Button[n].b = Red.pixel;
   Button[n].f = WhitePix;
   strcpy(Button[n].t[1],"Abort");
   Button[n].Action = Abort;

   // Clear orography

   ++n;
   // Button[n].x = Button[n-1].x;
   // Button[n].y = Button[n-1].y; // set in InitFrames
   Button[n].w = 6 * FixFontWidth;
   Button[n].h = 2 * FixFontHeight;
   Button[n].b = DarkBlue.pixel;
   Button[n].f = Yellow.pixel;
   strcpy(Button[n].t[1],"Clear");
   Button[n].Action = ClearOro;
}

/* ================== */
/* CreateModeSelector */
/* ================== */

void CreateModeSelector(int k)
{
   int i,j,m,n;
   int x,y,len,width,height;
   int r,dx,dy,xp,yp;
   char Text[20];

   pix = Frame[k].pixmap = XCreatePixmap(display,Cow,Frame[k].w,Frame[k].h,ScreenD);

   if (!Ampli) Ampli = malloc(DimSE * sizeof(int));
   if (!ModeX) ModeX = malloc(DimSE * sizeof(int));
   if (!ModeY) ModeY = malloc(DimSE * sizeof(int));
   if (!ModeM) ModeM = malloc(DimSE * sizeof(int));
   if (!ModeN) ModeN = malloc(DimSE * sizeof(int));

   dxsh = Frame[k].w / (DimTr + 3);
   dxs2 = dxsh / 2;
   dx   = dxsh;
   r  = (dx+dx) / 3;
   ModeRadiusSq = r*r;

   for (m=0,i=0 ; m <= DimTr ; ++m)
   {
      y = ModFontHeight + 4 + m * dx;
      for (n=m ; n <= DimTr ; ++n,++i)
      {
         x = dx/2 + n * dx;
         Ampli[i] = 1;
         ModeX[i] = x + dx/2;
         ModeY[i] = y + dx/2;
         ModeM[i] = m;
         ModeN[i] = n;
      }
   }
   i  = DimSH; // Button for all modes on
   x  = dx - ModFontWidth/2 + 1;
   y  = ModFontHeight + ModFont->descent + 2 + 20 * dx;
   Ampli[i] = 0;
   ModeX[i] = x + dx/2;
   ModeY[i] = y + dx/2;
   i++;       // Button for all modes off
   y += dx;
   Ampli[i] = 0;
   ModeX[i] = x + dx/2;
   ModeY[i] = y + dx/2;

   /* Draw mode legend */
   
   dx = dxsh;
   XSetForeground(display,gc,Red.pixel);
   XSetBackground(display,gc,BlackPix);
   for (i=0 ; i < 21 ; i+=2)
   {
      sprintf(Text,"%d",i);
      len    = strlen(Text);
      width  = XTextWidth(ModFont,Text,len);
      height = ModFont->ascent + ModFont->descent;
      xp     = dx + i * dx - width/2 - ModFontWidth/2 + 1;
      yp     = ModFontHeight;
      XDrawImageString(display,pix,gc,xp,yp,Text,len);
   }
   XSetForeground(display,gc,Blue.pixel);
   for (i=0 ; i < 21 ; i+=2)
   {
      sprintf(Text,"%d",i);
      len    = strlen(Text);
      width  = XTextWidth(ModFont,Text,len);
      height = ModFont->ascent + ModFont->descent;
      xp     = Frame[k].x + Frame[k].w - width - 2;
      yp     = 2 * ModFontHeight + i * dx;
      XDrawImageString(display,pix,gc,xp,yp,Text,len);
   }
}


void ShowButtons(void)
{
   int i,k,l,x,y,d;

   XSetFont(display, gc, FixFont->fid);
   for (k=0 ; k < DIMBUTTON ; ++k)
   {
      if (k == 4 && SelAno->hide) continue;
      d = Button[k].h / 4;
      XSetForeground(display,gc,Button[k].b);
      XFillRectangle(display,Cow,gc,Button[k].x,Button[k].y,Button[k].w,Button[k].h);
      XSetForeground(display,gc,Button[k].f);
      XSetBackground(display,gc,Button[k].b);
      for (i=0 ; i < 3 ; ++i)
      if (Button[k].t[i][0])
      {
         l = strlen(Button[k].t[i]);
         x = Button[k].x + (Button[k].w - l * FixFontWidth) / 2;
         y = Button[k].y + i * d + FixFont->ascent + d/2;
         XDrawImageString(display,Cow,gc,x,y,Button[k].t[i],l);
      }
   }
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

void InitColors(void)
{
/*
   Visual *visual;
   int count,st;
   XStandardColormap *best_map_info;
   XStandardColormap bmi[8];

   visual = DefaultVisual(display,0);
   memset(bmi,0,8*sizeof(XStandardColormap));
   best_map_info = bmi;

   count = 0;
   st = XGetRGBColormaps(display,RootWindow(display,0),&best_map_info,&count,XA_RGB_DEFAULT_MAP);
   {
      printf("status = %d  count=%d\n",st,count);
      printf("colormap = %p\n",bmi[0].colormap);
      printf("mapping = %x/%x/%x\n",bmi[0].red_max,bmi[0].green_max,bmi[0].blue_max);
   }
*/

   // colormap = best_map_info->colormap;

   AllocateColorCells(OroStrip);
   AllocateColorCells(OroMarsStrip);
   AllocateColorCells(GibbStrip);
   AllocateColorCells(TStrip);

   XAllocNamedColor(display,colormap,"red"        ,&Red       ,&Dummy);
   XAllocNamedColor(display,colormap,"green"      ,&Green     ,&Dummy);
   XAllocNamedColor(display,colormap,"blue"       ,&Blue      ,&Dummy);
   XAllocNamedColor(display,colormap,"grey"       ,&Grey      ,&Dummy);
   XAllocNamedColor(display,colormap,"hot pink"   ,&LightRed  ,&Dummy);
   XAllocNamedColor(display,colormap,"dark red"   ,&DarkRed   ,&Dummy);
   XAllocNamedColor(display,colormap,"light blue" ,&LightBlue ,&Dummy);
   XAllocNamedColor(display,colormap,"dark blue"  ,&DarkBlue  ,&Dummy);
   XAllocNamedColor(display,colormap,"light green",&LightGreen,&Dummy);
   XAllocNamedColor(display,colormap,"dark green" ,&DarkGreen ,&Dummy);
   XAllocNamedColor(display,colormap,"yellow"     ,&Yellow    ,&Dummy);
   XAllocNamedColor(display,colormap,"cyan"       ,&Cyan      ,&Dummy);

   TextC = Yellow.pixel;
   HeadC = Cyan.pixel;
}

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

// Convert an RGB value to an X11 Pixel

unsigned long create_pixel(long red, long green, long blue)
{
   if (ScreenD == 24)       // 24 bit true color
   {
       return blue | green << 8 | red << 16;
   }
   else if (ScreenD == 16)  // 16 bit RGB 565
   {
      return blue >> 3 | (green >> 2) << 5 | (red >> 3) << 11;

   }
   else return 0;
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
   int k;
   fread(&k,sizeof(k),1,fpi);
   if (BigEndian) SwapIEEE32((char *)&k);
   return k;
}

short ReadSHORT(FILE *fpi)
{
   short k;
   fread(&k,sizeof(k),1,fpi);
   if (BigEndian) SwapIEEE16((char *)&k);
   return k;
}

struct BMIstruct ImageBMI;

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
      printf("BMI ClrImpo  = %d\n",ImageBMI.ClrImportant);
      printf("Pad Bytes    = %d\n",PadBytes);
      printf("Pad Width    = %d\n",PadWidth);
   }

   PicBytes = (bpp * ImageBMI.Width + PadBytes) * ImageBMI.Height;
   ImgBytes = 4 * PadWidth * ImageBMI.Height;
   ei->d = calloc(ImgBytes,1);
   BuffImageData = malloc(PicBytes);

   fseek(fp,OffBits,SEEK_SET);
   n = fread(BuffImageData,1,PicBytes,fp);
   fclose(fp);  

   if (Debug)
   {
      printf("Size Bytes   = %d\n",PicBytes);
      printf("Read Bytes   = %d\n",n);
   }
   
   ei->X = XCreateImage(display,CopyFromParent,ScreenD,ZPixmap,0,
                  ei->d,PadWidth,ImageBMI.Height,8,0);

   for (y = 0 ; y < ImageBMI.Height ; ++y)
   {
      for (x = 0 ; x < ImageBMI.Width  ; ++x)
      {
          i = bpp * x + (ImageBMI.Height - 1 - y) * (ImageBMI.Width*bpp+PadBytes);
          b = BuffImageData[i  ];
          g = BuffImageData[i+1];
          r = BuffImageData[i+2];
          XPutPixel(ei->X, x, y, create_pixel(r,g,b));
      }
   }
   free(BuffImageData);
   ei->w = ImageBMI.Width;
   ei->h = ImageBMI.Height;
}

void ReadLogo(int logo, char *filename)
{
   struct MapImageStruct MapI;

   ReadImage(&MapI,filename);

   Logo[logo].w = MapI.w;
   Logo[logo].h = MapI.h;
   Logo[logo].X = MapI.X;

   if (logo)
   {
      Logo[logo].x = Logo[logo-1].x + Logo[logo-1].w;
      Logo[logo].y = Logo[logo-1].y;
   }
}

void ShowCopyright(void)
{
   int x,y;
   x = CowSizeHints.min_width - 18.5 * ModFontWidth;
   y = CowSizeHints.min_height - ModFontHeight/2;
   XSetFont(display, gc, ModFont->fid);
   XSetForeground(display,gc,Grey.pixel);
   XSetBackground(display,gc,BlackPix);
   XDrawImageString(display,Cow,gc,x,y,"Image Credit: NASA",18);
}

int redaco;

int RedrawControlWindow(void)
{
   int i,l,x,y;
   double xrf,yrf;
   struct SelStruct *Sel;

   XWindowAttributes CurAtt;

   if (Debug) printf("Redraw %d\n",redaco++);
   XGetWindowAttributes(display,Cow,&CurAtt);
   WinXSize = CurAtt.width;
   WinYSize = CurAtt.height;
   
   XSetWindowBackground(display,Cow,WinBG);
   XClearWindow(display,Cow);

   // Paint some stars on the sky

   xrf = (double)WinXSize / pow(2.0,31),
   yrf = (double)WinYSize / pow(2.0,31);

   XSetForeground(display,gc,WhitePix);
   srandom(Seed);
   for (y = 0 ; y < ((WinXSize * WinYSize) >> 8)  ; ++y)
      XDrawPoint(display,Cow,gc,xrf*random(),yrf*random());

   XSetForeground(display,gc,BlackPix);
   XSetBackground(display,gc,WhitePix);

   for (Sel = &SelStart ; Sel ; Sel = Sel->Next) ShowSelection(Sel);

   if (CursorSel) CursorOn();

   ShowOrography();
   ShowFrame1();
   ShowButtons();
   for (l=0 ; l < 3 ; ++l)
   {
      XPutImage(display,Cow,gc,Logo[l].X,0,0,Logo[l].x,Logo[l].y,Logo[l].w,Logo[l].h);
   }
   if (Model == PLASIM) ShowCopyright();
   XSync(display,0);
   return 0;
}


void ChangeResolution(int NewRes)
{
   int i;
   struct SelStruct *Sel;

   for (i=RES_T21 , Sel = SelRes ; i <= RES_T42; ++i , Sel = Sel->Next)
   {
      if (i == NewRes) Sel->iv = 1;
      else             Sel->iv = 0;
   }
   Resolution = NewRes;
}


void UpdateResolution(void)
{
   int i;
   struct SelStruct *Sel;

   if (SelRes)
   for (i=RES_T21 , Sel = SelRes ; i <= RES_T42; ++i , Sel = Sel->Next)
   {
      if (Sel && Sel->iv == 1) Resolution = i;
   }
}


char *vcn[6] =
{
   "StaticGray",
   "GrayScale",
   "StaticColor",
   "PseudoColor",
   "TrueColor",
   "DirectColor"
};


void InitGUI(void)
{
   int argc = 1;
   int i,j,k;
   char *WinconTitle1 = {"MoSt - Model Suite (17) - University of Hamburg"};
   unsigned long valuemask = 0; /* ignore XGCvalues and use defaults */
   FILE *ftp;
   FILE *xpp;
   XVisualInfo vinfo;

   XGCValues values;
   XEvent Event;
   XWindowAttributes CurAtt;

   if ((display=XOpenDisplay(display_name)) == NULL)
   {
      fprintf(stderr,"%s: cannot connect to X server %s\n", 
              progname, XDisplayName(display_name));
      exit(1);
   }
   ScreenN  = DefaultScreen(display);
   ScreenD  = XDefaultDepth(display,ScreenN);
   ScreenW  = DisplayWidth (display,ScreenN);
   ScreenH  = DisplayHeight(display,ScreenN);
   if (ScreenHeight) ScreenH = ScreenHeight;

   if (Debug)
   {
      i=5;
      while (!XMatchVisualInfo(display,ScreenN,ScreenD,i--,&vinfo));
      printf("Found visual %s at depth %d\n",vcn[++i],ScreenD);
      printf("Red   mask = %8lx\n",vinfo.red_mask);
      printf("Green mask = %8lx\n",vinfo.green_mask);
      printf("Blue  mask = %8lx\n",vinfo.blue_mask);
   }

   BlackPix = BlackPixel(display,ScreenN);
   WhitePix = WhitePixel(display,ScreenN);
   SmallScreen = ScreenH < 768;

   wm_hints.initial_state = NormalState;
   wm_hints.input = True;
   wm_hints.flags = StateHint | InputHint;
   
   class_hints.res_name = progname;
   class_hints.res_class = "MOST";

   Delwin = XInternAtom(display,"WM_DELETE_WINDOW",0);


   LoadFonts();

   /* Setup control window */

   CowSizeHints.flags      = PPosition | PSize | PMinSize | PMaxSize;
   CowSizeHints.min_width  = MIN(WINDOW_WIDTH,ScreenW);
   CowSizeHints.min_height = MIN(740,ScreenH);
   CowSizeHints.max_width  = CowSizeHints.min_width;
   CowSizeHints.max_height = CowSizeHints.min_height;
   Cow = XCreateSimpleWindow(display,RootWindow(display,ScreenN), 
               0,0,
               CowSizeHints.max_width,CowSizeHints.max_height,
               4,BlackPix,WhitePix);
   XStringListToTextProperty(&WinconTitle1,1,&WinconName1);
   XSetWMProtocols(display,Cow,&Delwin,1);
   XSetWMProperties(display,Cow,&WinconName1,NULL,
		NULL,0,&CowSizeHints,&wm_hints,&class_hints);
   XSelectInput(display,Cow,ButtonPressMask | ButtonReleaseMask | 
                PointerMotionMask | KeyPressMask | ExposureMask);
   XMapWindow(display,Cow);


   XGetWindowAttributes(display,Cow,&CurAtt);
   WinXSize = CurAtt.width;
   WinYSize = CurAtt.height;

   /* Prepare GC */

   gc = XCreateGC(display, Cow, valuemask, &values);
   XSetFont(display, gc, FixFont->fid);
   colormap = XDefaultColormap(display,ScreenN);
   XSetForeground(display,gc,BlackPix);
   XSetBackground(display,gc,WhitePix);

   Mok = XGetModifierMapping(display);

   InitColors();

   /* Look for Num_Lock key */

   for (j=0 ; j < 8 ; ++j)
   {
      k = Mok->modifiermap[j * Mok->max_keypermod];
      if (XkbKeycodeToKeysym(display,k,0,0) == XK_Num_Lock   ) NumLockMask    = 1 << j;
      if (XkbKeycodeToKeysym(display,k,0,0) == XK_Mode_switch) ModeSwitchMask = 1 << j;
   }
   if (Debug) printf("NumLockMask    = %02X\n",NumLockMask   );
   if (Debug) printf("ModeSwitchMask = %02X\n",ModeSwitchMask);
   if (Debug)
   for (j=0 ; j < 8 ; ++j)
   {
      printf("Mok [%d] %d:",Mok->max_keypermod,j);
      for (i=0 ; i < Mok->max_keypermod ; ++i)
      {
	 k = Mok->modifiermap[i+j*Mok->max_keypermod];
         if (XkbKeycodeToKeysym(display,k,0,0) != NoSymbol)
         printf("  %4x %-16.16s   ",k,XKeysymToString(XkbKeycodeToKeysym(display,k,0,0)));
      }
      printf("\n");
   }
   XDisplayKeycodes(display,&EdiFirstKey,&EdiLastKey);
   EdiKeymap = XGetKeyboardMapping(display,EdiFirstKey,EdiLastKey-EdiFirstKey+1,&EdiSymsPerKey);

   /* Sun keymaps don't have lower case entries */

   for (j=EdiFirstKey ; j <= EdiLastKey ; ++j)
   {
      k = (j - EdiFirstKey) * EdiSymsPerKey;
      if (EdiKeymap[k] >= 'A' && EdiKeymap[k] < 0xE0 && EdiKeymap[k+1] == 0)
      {
         EdiKeymap[k+1] = EdiKeymap[k];
         EdiKeymap[k]  += 0x20;
      }
   }

   xpp = fopen("Beginner","r"); // Expert mode ?
   if (xpp)
   {
      Expert = 0;
      fclose(xpp);
   }

   xpp = fopen("lsg/src/lsgmod.f90","r"); // LSG there ?
   if (xpp)
   {
      LsgEnabled = 1;
      fclose(xpp);
   }

   xpp = fopen("puma/src/mpimod_multi.f90","r"); // Multirun module there ?
   if (xpp)
   {
      MultirunEnabled = 1;
      fclose(xpp);
   }

   // Read name of MPI execute command

   xpp = fopen("most_compiler_mpi","r"); // MPI installed ?
   if (xpp)
   {
      fgets(Buffer,LINEMAX,xpp);
      if (Buffer[strlen(Buffer)-1] == 10) Buffer[strlen(Buffer)-1] = 0;
      if (Buffer[strlen(Buffer)-1] == 13) Buffer[strlen(Buffer)-1] = 0;
      if (!strncmp(Buffer,"MPI_RUN=",8))  strcpy(mpirun,Buffer+8);
      fclose(xpp);
   }

   InitLogo();
   ReadLogo(0,"images/KC-Logo_RGB.bmp");
   ReadLogo(1,"images/puma.bmp");
   ReadLogo(2,"images/plasim.bmp");
   ReadImage(&MapHRE,"images/earth.bmp");
   ReadImage(&MapHRM,"images/mars.bmp");
   ReadImage(&MapLRK,"images/Kepler-16.bmp");
   ReadImage(&MapLRL,"images/habit360x180.bmp");
   opbox_y = 8 + 64 + 2 * FixFontHeight;
   InitSelections();

   InitNamelist();
   
   ChangeModel(PLASIM);
   NamelistSelector(PLASIM);
   ChangeModel(SAM);
   NamelistSelector(SAM);
   ChangeModel(PUMA);
   //AddPumaNamelist();
   NamelistSelector(PUMA);

   if (ReadSettings(cfg_file))
   {
      UpdateResolution();
      UpdateSelections(&SelStart);
   }

   InitFrames();
   CalcFrame(nlat);
   CreateModeSelector(1);
   InitButtons();
   MapLRE.w = Frame[1].w;
   MapLRE.h = Frame[1].h;
   MapLRK.w = Frame[1].w;
   MapLRK.h = Frame[1].h;
   MapLRM.w = Frame[1].w;
   MapLRM.h = Frame[1].h;
   MapLRL.w = Frame[0].w;
   MapLRL.h = Frame[0].h;
   MapLRM.l = -90;
   AzimuthalImage(&MapHRE,&MapLRE);
   AzimuthalImage(&MapHRM,&MapLRM);
   RedrawControlWindow();
}


int HitBox(struct SelStruct *Sel)
{
   return (
   (WinEvent.xbutton.button == Button1) &&
   (WinEvent.xbutton.x >= Sel->x) &&
   (WinEvent.xbutton.x <  Sel->x + Sel->w) &&
   (WinEvent.xbutton.y >= Sel->y) &&
   (WinEvent.xbutton.y <  Sel->y + Sel->h));
}


int HitButton(int k)
{
   return (
   (WinEvent.xbutton.button == Button1) &&
   (WinEvent.xbutton.x >= Button[k].x) &&
   (WinEvent.xbutton.x <  Button[k].x + Button[k].w) &&
   (WinEvent.xbutton.y >= Button[k].y) &&
   (WinEvent.xbutton.y <  Button[k].y + Button[k].h));
}


int HitFrame(int k)
{
   return (
   (WinEvent.xbutton.x >= Frame[k].x) &&
   (WinEvent.xbutton.x <  Frame[k].x + Frame[k].w) &&
   (WinEvent.xbutton.y >= Frame[k].y) &&
   (WinEvent.xbutton.y <  Frame[k].y + Frame[k].h));
}


int OnFrame(int i)
{
   if (HitFrame(i))
   {
      FrameNo = i;
      Frame[i].Action();
      return 1;
   }
   return 0;
}

void OnMouseClick(void)
{
   int i,j,cp;
   struct SelStruct *Sel;

   /* Check Action Buttons */

   for (i=0 ; i < DIMBUTTON ; ++i)
   {
      if (HitButton(i))
      {
         Button[i].Action();
         return;
      }
   }

   /* Check for model switch */

   for (i = PUMA , Sel = SelMod ; i < MODELS ; ++i , Sel = Sel->Next)
   if (HitBox(Sel) && Model != i)
   {
      if (Debug) printf("Change model from %d to %d\n",Model,i);
      ChangeModel(i);
      CalcFrame(Latitudes);
      return;
   }

   /* Continue with all other boxes */

   for ( ; Sel ; Sel = Sel->Next)
   {
      if (Sel->type == SEL_CHECK && HitBox(Sel) && !Sel->no)
      {
         cp = -1;
         for (i=0 ; i < PLANETS ; ++i) if (Sel == SelPlanet[i]) cp = i;
         if (cp >= 0 && Planet != cp) ChangePlanet(cp);
         else
         {
            Sel->iv = !Sel->iv;
            if (Sel->piv) *Sel->piv = Sel->iv;
         }
      }
      else if ((Sel->type == SEL_INT || Sel->type == SEL_REAL) && HitBox(Sel))
      {
         // printf("Target is <%s>\n",Sel->text);
         if (Sel != CursorSel) FinishLine();
         CursorSel = Sel;
         CursorCol = (WinEvent.xbutton.x - Sel->x) / FixFontWidth;
         if (CursorCol > Sel->edco-1) CursorCol = Sel->edco-1;
         CursorOn();
      }
   }

   /* Allow only one horizontal resolution */

   if (SelRes)
   {
      for (i = RES_T21 , Sel = SelRes ; i <= RES_T42 ; ++i , Sel = Sel->Next)
      if (HitBox(Sel))
      {
         ChangeResolution(i);
         return;
      }
   }
}


int EdiDecodeKey(void)
{
   int ShiftStatus;
   int ShiftIndex;
   int KeyCode;
   int KeyIndex;

   ShiftStatus =  WinEvent.xkey.state;
   KeyIndex    = (WinEvent.xkey.keycode-EdiFirstKey) * EdiSymsPerKey;
   KeyCode     = EdiKeymap[KeyIndex]; // Code with no modifiers
   if (ShiftStatus & 0x04)            // Control
   {
      if (KeyCode >= 0x40 && KeyCode < 0x80) KeyCode &= 0x1F;
   }
   else if (ShiftStatus & NumLockMask &&
            EdiKeymap[KeyIndex+1] >= XK_KP_Separator &&  // Numlock
	    EdiKeymap[KeyIndex+1] <= XK_KP_9)            // PC Keypad
      KeyCode = EdiKeymap[KeyIndex+1] - XK_KP_Space;
   else if (ShiftStatus & NumLockMask &&
            EdiKeymap[KeyIndex+2] >= XK_KP_Separator &&  // Numlock
	    EdiKeymap[KeyIndex+2] <= XK_KP_9)            // SUN Keypad
      KeyCode = EdiKeymap[KeyIndex+2] - XK_KP_Space;
   else
   {
      if (ShiftStatus & ModeSwitchMask) ShiftIndex = 2; // Alt Gr
      else ShiftIndex = ShiftStatus & 1;      // Normal & Shift
      if (ShiftIndex == 1 && KeyCode >= XK_F1 && KeyCode <= XK_F12) KeyCode += 12;
      else KeyCode = EdiKeymap[KeyIndex + ShiftIndex];
   }
   if (ShiftStatus & 0x02) // Caps Lock
   {
      if (KeyCode >= 'a'  && KeyCode <= 'z' ) KeyCode -= 0x20; // ASCII
      if (KeyCode >= 0xE0 && KeyCode <= 0xFD) KeyCode -= 0x20; // Latin-1
   }
   if (KeyCode == XK_KP_Left      ) KeyCode = XK_Left;
   if (KeyCode == XK_KP_Right     ) KeyCode = XK_Right;
   if (KeyCode == XK_KP_Up        ) KeyCode = XK_Up;
   if (KeyCode == XK_KP_Down      ) KeyCode = XK_Down;
   if (KeyCode == XK_KP_Home      ) KeyCode = XK_Home;
   if (KeyCode == XK_KP_End       ) KeyCode = XK_End  ;
   if (KeyCode == XK_KP_Page_Up   ) KeyCode = XK_Page_Up;
   if (KeyCode == XK_KP_Page_Down ) KeyCode = XK_Page_Down;
   if (KeyCode == XK_KP_Enter     ) KeyCode = XK_Return;
   if (KeyCode == XK_KP_Add       ) KeyCode = '+';
   if (KeyCode == XK_KP_Subtract  ) KeyCode = '-';
   if (KeyCode == XK_KP_Multiply  ) KeyCode = '*';
   if (KeyCode == XK_KP_Divide    ) KeyCode = '/';
   if (Debug)
   {
      if (KeyCode >= ' ' && KeyCode <= 255) printf("'%c' ",KeyCode);
      else                                  printf("--- ");
      printf("Key [%x] State [%2x]\n",WinEvent.xkey.keycode,WinEvent.xkey.state);
      for (ShiftIndex = 0 ; ShiftIndex < EdiSymsPerKey ; ++ShiftIndex)
      {
	 if (EdiKeymap[KeyIndex + ShiftIndex] && XKeysymToString(EdiKeymap[KeyIndex + ShiftIndex]))
         {
            printf("   %d: <%x> '%c' (%s)\n",
	       ShiftIndex,(unsigned int)EdiKeymap[KeyIndex + ShiftIndex],
	                  (int)EdiKeymap[KeyIndex + ShiftIndex],
	       XKeysymToString(EdiKeymap[KeyIndex + ShiftIndex]));
	 }
      }
      printf("\n");
   }
   return KeyCode;
}

void EditUpLine(void)
{
   FinishLine();
   if (CursorSel == NULL || CursorSel->Prev == NULL) return;
   if (CursorSel->Prev->type != SEL_INT && CursorSel->Prev->type != SEL_REAL) return;
   CursorSel = CursorSel->Prev;
   if (CursorCol > CursorSel->edco-1) CursorCol = CursorSel->edco-1;
}

void EditDownLine(void)
{
   FinishLine();
   if (CursorSel == NULL || CursorSel->Next == NULL) return;
   if (CursorSel->Next->type != SEL_INT && CursorSel->Next->type != SEL_REAL) return;
   CursorSel = CursorSel->Next;
   if (CursorCol > CursorSel->edco-1) CursorCol = CursorSel->edco-1;
}

void EditReturnLine(void)
{
   CursorCol = 0;
   EditDownLine();
}

void EditLeftChar(void)
{
   if (CursorCol > 0) --CursorCol;
}

void EditRightChar(void)
{
   if (CursorCol < CursorSel->edco-1) ++CursorCol;
}

void EmitChar(int ch)
{
   if (CursorSel)
   {
      CursorSel->teva[CursorCol] = ch;
      EditRightChar();
   }
} /* EmitChar */


void DeleteRightChar(void)
{
   int i;

   for (i=CursorCol ; i < CursorSel->edco-1 ; ++i)
      CursorSel->teva[i] = CursorSel->teva[i+1];
   CursorSel->teva[CursorSel->edco-1] = ' ';
}


void DeleteLeftChar(void)
{
   if (CursorCol == 0) return;
   EditLeftChar();
   DeleteRightChar();
}


void InsertBlank(void)
{
   int i;

   for (i = CursorSel->edco-1 ; i > CursorCol ; --i)
      CursorSel->teva[i] = CursorSel->teva[i-1];
   CursorSel->teva[CursorCol] = ' ';
}


void ProcChar(int Key)
{
   if (CursorSel)
   switch (Key)
   {
      case XK_Up       : EditUpLine();           break;
      case XK_Left     : EditLeftChar();         break;
      case XK_Right    : EditRightChar();        break;
      case XK_Down     : EditDownLine();         break;
      case XK_Return   : EditReturnLine();       break;
      case XK_BackSpace: DeleteLeftChar();       break;
      case XK_Delete   : DeleteRightChar();      break;
      case XK_Insert   : InsertBlank();          break;
   }
}


void LoopGUI(void)
{
   int w;

   while (1)
   {
      XNextEvent(display,&WinEvent);

      switch (WinEvent.type)
      {
         case ClientMessage: return;
         case MotionNotify : if (OnFrame(0)) ShowOrography();
                             if (Model == PUMA && OnFrame(1)) ShowModeSelector();
                             break;
         case ButtonRelease: Button1Down = 0; if (OnFrame(0)) ShowOrography(); break;
         case ButtonPress  : Button1Down = 1;
                             if (!OnFrame(0) && !OnFrame(1)) OnMouseClick();
         case Expose       : RedrawControlWindow(); break;
         case KeyPress     : EdiKeyCode = EdiDecodeKey();
                             if (EdiKeyCode >= ' ' && EdiKeyCode < 256) EmitChar(EdiKeyCode);
                             else                                       ProcChar(EdiKeyCode);
                             RedrawControlWindow();
                             break;

      }
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

int main(int argc, char *argv[])
{
   int ia;
   for (ia = 1 ; ia < argc ; ++ia)
   {
           if (!strcmp(argv[ia],"-d")) Debug = 1;
      else if (!strcmp(argv[ia],"-i")) cfg_file[0] = 0;
      else if (!strcmp(argv[ia],"-h")) ScreenHeight = atoi(argv[++ia]);
      else strncpy(cfg_file,argv[1],sizeof(cfg_file));
   }
   CurrentDate = time(NULL);
   BigEndian = CheckEndianess();
   InitGUI();
   LoopGUI();
   return 0;
}
