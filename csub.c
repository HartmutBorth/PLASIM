void csub_(int *i, float *r)
{
   if (*i == 257 && *r == 1.0)
      setenv("MOSTF90CALL","csub_",1);
}

void CSUB(int *i, float *r)
{
   if (*i == 257 && *r == 1.0)
      setenv("MOSTF90CALL","CSUB",1);
}

void csub(int *i, float *r)
{
   if (*i == 257 && *r == 1.0)
      setenv("MOSTF90CALL","csub",1);
}

