#ifndef _XMALLOC_H
#define _XMALLOC_H

#include <stdlib.h>

void *xmalloc(size_t size);
void *xcalloc(size_t n, size_t size);
void *xrealloc(void *ptr, size_t size);

#endif // _XMALLOC_H