#include <stdlib.h>
#include <stdio.h>

#define NLEV 5

int DimX;
int DimY;
int DimT;
int head[10];

float *F;

int main(int argc, char *argv[])
{
   int i,j,jlev;

   FILE *ifp;
   FILE *ofp;

   ifp = fopen(argv[1],"r");
   ofp = fopen(argv[2],"w");

   for (jlev = 0 ; jlev < NLEV ; ++jlev)
   {
      fread(head,sizeof(int),10,ifp);
   
      DimX = head[4];
      DimY = head[5];
      DimT = DimX * DimY;
   
      F = (float *)malloc((DimT+2) * sizeof(float));
      fread(F,sizeof(float),DimT+2,ifp);
   
      for (i = 1 ; i < 9 ; ++i)
         fprintf(ofp,"%9d",head[i]);
      fprintf(ofp,"\n");
   
      for (j = 1 ; j <= DimT ; j += 8)
      {
         for (i = j ; i < j + 8 && i <= DimT ; ++i)
            fprintf(ofp,"%9.2f",F[i]);
         fprintf(ofp,"\n");
      }
   }
   return 0;
}
