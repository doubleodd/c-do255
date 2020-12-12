/*
 * This file is meant to be included, not compiled by itself.
 * Caller must have defined the following prior to inclusion:
 *
 *  - defined macro CURVE to the curve identifier (e.g. do255s)
 *
 * This file does the following:
 *
 *  - include "do255.h"
 *  - define the CN() macro to decorate names with the curve identifier
 */

#include <string.h>
#include "do255.h"

#define CN(x)            CNN(CURVE, x)
#define CNN(cname, x)    CNN_(cname, x)
#define CNN_(cname, x)   cname ## _ ## x

#if defined __GNUC__ || defined __clang__
#define UNUSED   __attribute__((unused))
#define FORCE_INLINE   __attribute__((always_inline))
#define NO_INLINE      __attribute__((noinline))
#else
#define UNUSED
#define FORCE_INLINE
#define NO_INLINE
#endif
