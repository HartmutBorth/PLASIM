#include <stdlib.h>
#include <stdio.h>

int DimX;
int DimY;
int DimT;
int head[10];
int Count;

float *F;

int main(int argc, char *argv[])
{
   int i,j;

   FILE *ifp;
   FILE *ofp;

   ifp = fopen(argv[1],"r");
   ofp = fopen(argv[2],"w");

   fread(head,sizeof(int),10,ifp);

   DimX = head[5];
   DimY = head[6];
   DimT = DimX * DimY;

   F = (float *)malloc((DimT+2) * sizeof(float));

   while (!feof(ifp))
   {
      ++Count;
      fread(F,sizeof(float),DimT+2,ifp);
   
      for (i = 1 ; i < 9 ; ++i) fprintf(ofp,"%9d",head[i]);
      fprintf(ofp,"\n");
   
      for (j = 1 ; j <= DimT ; j += 8)
      {
         for (i = j ; i < j + 8 && i <= DimT ; ++i)
            fprintf(ofp," %8.3f",F[i]);
         fprintf(ofp,"\n");
      }
      fread(head,sizeof(int),10,ifp);
   }
   if (Count == 0) printf("\n*** Error - no service array found ***\n");
   else if (Count == 1)
      printf("\nConverted one array from binary service (srv) to ASCII service (sra)\n");
   else
      printf("\nConverted %d arrays from binary service (srv) to ASCII service (sra)\n",Count);
   return 0;
}
