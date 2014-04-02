#include <stdlib.h>
#include <stdio.h>

FILE * fp;

int nstep;
int ostep;

int main(int argc, char *argv[])
{
   nstep = atoi(argv[2]);
   fp = fopen(argv[1],"r+");
   fseek(fp,28,0);
   fread(&ostep,4,1,fp);
   fseek(fp,28,0);
   fwrite(&nstep,4,1,fp);
   fclose(fp);
   printf("\nNSTEP changed from %d to %d in file <%s>\n",ostep,nstep,argv[1]);
   return 0;
}
