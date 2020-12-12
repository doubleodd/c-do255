#include <string.h>
#include "do255.h"
#include "do255_alg.h"

#define CN(x)    CNN(x)
#define CNN(x)   do255s_ ## x

#define DOM_ECDH     "do255s-ecdh:"
#define DOM_SIGN_K   "do255s-sign-k:"
#define DOM_SIGN_E   "do255s-sign-e:"

#include "alg.c"
