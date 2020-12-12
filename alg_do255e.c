#include <string.h>
#include "do255.h"
#include "do255_alg.h"

#define CN(x)    CNN(x)
#define CNN(x)   do255e_ ## x

#define DOM_ECDH     "do255e-ecdh:"
#define DOM_SIGN_K   "do255e-sign-k:"
#define DOM_SIGN_E   "do255e-sign-e:"

#include "alg.c"
