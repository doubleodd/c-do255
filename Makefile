CC = clang
CFLAGS = -Wall -Wextra -Wundef -Wshadow -O2 -march=skylake
LD = clang
LDFLAGS =
LIBS =

OBJ_TESTS = sha3.o test_do255.o
OBJ_DO255E_BMI2 = do255e_bmi2.o
OBJ_DO255E_W64 = do255e_w64.o
OBJ_DO255E_W32 = do255e_w32.o
OBJ_DO255S_BMI2 = do255s_bmi2.o
OBJ_DO255S_W64 = do255s_w64.o
OBJ_DO255S_W32 = do255s_w32.o

OBJ_ALG_DO255E = alg_do255e.o
OBJ_ALG_DO255S = alg_do255s.o

all: test_do255_bmi2 test_do255_w64 test_do255_w32

clean:
	-rm -f $(OBJ_DO255E_BMI2) $(OBJ_DO255S_BMI2) $(OBJ_DO255E_W64) $(OBJ_DO255S_W64) $(OBJ_DO255E_W32) $(OBJ_DO255S_W32) $(OBJ_ALG_DO255E) $(OBJ_ALG_DO255S) $(OBJ_TESTS) test_do255_bmi2 test_do255_w64 test_do255_w32

test_do255_bmi2: $(OBJ_DO255E_BMI2) $(OBJ_DO255S_BMI2) $(OBJ_ALG_DO255E) $(OBJ_ALG_DO255S) $(OBJ_TESTS)
	$(LD) $(LDFLAGS) -o test_do255_bmi2 $(OBJ_DO255E_BMI2) $(OBJ_DO255S_BMI2) $(OBJ_ALG_DO255E) $(OBJ_ALG_DO255S) $(OBJ_TESTS) $(LIBS)

test_do255_w32: $(OBJ_DO255E_W32) $(OBJ_DO255S_W32) $(OBJ_ALG_DO255E) $(OBJ_ALG_DO255S) $(OBJ_TESTS)
	$(LD) $(LDFLAGS) -o test_do255_w32 $(OBJ_DO255E_W32) $(OBJ_DO255S_W32) $(OBJ_ALG_DO255E) $(OBJ_ALG_DO255S) $(OBJ_TESTS) $(LIBS)

test_do255_w64: $(OBJ_DO255E_W64) $(OBJ_DO255S_W64) $(OBJ_ALG_DO255E) $(OBJ_ALG_DO255S) $(OBJ_TESTS)
	$(LD) $(LDFLAGS) -o test_do255_w64 $(OBJ_DO255E_W64) $(OBJ_DO255S_W64) $(OBJ_ALG_DO255E) $(OBJ_ALG_DO255S) $(OBJ_TESTS) $(LIBS)

alg_do255e.o: alg_do255e.c alg.c do255.h do255_alg.h sha3.h
	$(CC) $(CFLAGS) -c -o alg_do255e.o alg_do255e.c

alg_do255s.o: alg_do255s.c alg.c do255.h do255_alg.h sha3.h
	$(CC) $(CFLAGS) -c -o alg_do255s.o alg_do255s.c

do255e_bmi2.o: do255e_bmi2.c do255.h support.c gf_bmi2.c gf_do255e_bmi2.c sqrt_do255e_w64.c pcore_w64.c padd_do255e_w64.c icore_w64.c scalar_do255e_w64.c pmul_base_w64.c pmul_do255e_w64.c pvrfy_do255e_w64.c pmap_do255e_w64.c
	$(CC) $(CFLAGS) -c -o do255e_bmi2.o do255e_bmi2.c

do255s_bmi2.o: do255s_bmi2.c do255.h support.c gf_bmi2.c gf_do255s_bmi2.c sqrt_do255s_w64.c pcore_w64.c padd_do255s_w64.c icore_w64.c scalar_do255s_w64.c pmul_base_w64.c pmul_do255s_w64.c lagrange_do255s_w64.c pvrfy_do255s_w64.c pmap_do255s_w64.c
	$(CC) $(CFLAGS) -c -o do255s_bmi2.o do255s_bmi2.c

do255e_w64.o: do255e_w64.c do255.h support.c gf_w64.c gf_do255e_w64.c sqrt_do255e_w64.c pcore_w64.c padd_do255e_w64.c icore_w64.c scalar_do255e_w64.c pmul_base_w64.c pmul_do255e_w64.c pvrfy_do255e_w64.c pmap_do255e_w64.c
	$(CC) $(CFLAGS) -c -o do255e_w64.o do255e_w64.c

do255s_w64.o: do255s_w64.c do255.h support.c gf_w64.c gf_do255s_w64.c sqrt_do255s_w64.c pcore_w64.c padd_do255s_w64.c icore_w64.c scalar_do255s_w64.c pmul_base_w64.c pmul_do255s_w64.c lagrange_do255s_w64.c pvrfy_do255s_w64.c pmap_do255s_w64.c
	$(CC) $(CFLAGS) -c -o do255s_w64.o do255s_w64.c

do255e_w32.o: do255e_w32.c do255.h support.c gf_w32.c gf_do255e_w32.c sqrt_do255e_w32.c pcore_w32.c padd_do255e_w32.c icore_w32.c scalar_do255e_w32.c pmul_base_w32.c pmul_do255e_w32.c pvrfy_do255e_w32.c pmap_do255e_w32.c
	$(CC) $(CFLAGS) -c -o do255e_w32.o do255e_w32.c

do255s_w32.o: do255s_w32.c do255.h support.c gf_w32.c gf_do255s_w32.c sqrt_do255s_w32.c pcore_w32.c padd_do255s_w32.c icore_w32.c scalar_do255s_w32.c pmul_base_w32.c pmul_do255s_w32.c lagrange_do255s_w32.c pvrfy_do255s_w32.c pmap_do255s_w32.c
	$(CC) $(CFLAGS) -c -o do255s_w32.o do255s_w32.c

sha3.o: sha3.c sha3.h
	$(CC) $(CFLAGS) -c -o sha3.o sha3.c

test_do255.o: test_do255.c sha3.h do255.h do255_alg.h
	$(CC) $(CFLAGS) -c -o test_do255.o test_do255.c
