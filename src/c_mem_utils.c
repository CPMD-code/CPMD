#include <stddef.h>

void *cMemPointsToDbl ( double *p, int elemSize, int shift )
{
  return p + elemSize * shift;
}

size_t cGetMemAddrs ( void *p )
{
  size_t *pp = (size_t *) p;
  return pp;
}

