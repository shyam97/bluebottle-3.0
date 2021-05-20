#include "tracertryout.h"

void tracer_init(void)
{
  tracercheck = 7;
}

void tracer_print(void)
{
  cuda_dom_pull();
  printf("What's up, fam? Tracercheck=%d\n",tracercheck);
  printf("%d,%d\n",dom[rank].Gcc.s1, dom[rank].Gcc.s2);	
  tracercheck++;
}
