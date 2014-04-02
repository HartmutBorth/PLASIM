#include <stdlib.h>
#include <stdio.h>

int fcw = 4;
int isz = 4;
int rsz = 4;
int ben = 0;

FILE *wp;

void check_fortran_integer(void)
{
   int b;
   int l;
   FILE *fp;
   fp = fopen("F90_INTEGER","r");
   fseek(fp,0,SEEK_END);
   l = ftell(fp);
        if (l == 16) isz = 8;
   else if (l == 20) fcw = 8;
   else if (l == 24)
   {
      isz = 8;
      fcw = 8;
   }
   rewind(fp);
   b = fgetc(fp);
   if (b == 0) ben = 1;
   fclose(fp);
}

void check_fortran_real(void)
{
   int b;
   int l;
   FILE *fp;
   fp = fopen("F90_REAL","r");
   fseek(fp,0,SEEK_END);
   l = ftell(fp);
   rsz = l - 2 * fcw;
   fclose(fp);
}

void print_info(void)
{
   wp = fopen("most_info.txt","w");
   fprintf(wp,"\nSystem info for <%s>\n",getenv("HOSTNAME"));
   fprintf(wp,"Architecture: %s\n",getenv("MOSTARCH"));
   if (ben) fprintf(wp,"Endian format             : big endian\n");
   else     fprintf(wp,"Endian format             : little endian\n");
   fprintf(wp,"FORTRAN control word size : %d bytes\n",fcw);
   fprintf(wp,"FORTRAN integer size      : %d bytes\n",isz);
   fprintf(wp,"FORTRAN real    size      : %d bytes\n",rsz);
   fprintf(wp,"C       int     size      : %ld bytes\n",sizeof(int));
   fprintf(wp,"C       float   size      : %ld bytes\n",sizeof(float));
   fprintf(wp,"C       long    size      : %ld bytes\n",sizeof(long));
   fprintf(wp,"C  long long    size      : %ld bytes\n",sizeof(long long));
   fclose(wp);
}

int main(int argc, char *argv[])
{
   check_fortran_integer();
   check_fortran_real();
   print_info();

   return 0;
}
