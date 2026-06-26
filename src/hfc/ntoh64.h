#ifndef __NTOH64_H
#define __NTOH64_H

#include <stdint.h>

/* Host/network byte order for 64-bit GEB timestamps (inline avoids duplicate symbols). */
static inline uint64_t ntoh64(uint64_t input)
{
  uint64_t rval;
  uint8_t *data = (uint8_t *)&rval;
  data[0] = input >> 56;
  data[1] = input >> 48;
  data[2] = input >> 40;
  data[3] = input >> 32;
  data[4] = input >> 24;
  data[5] = input >> 16;
  data[6] = input >> 8;
  data[7] = input >> 0;
  return rval;
}

#endif
